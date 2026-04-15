# =============================================================================
# 01_modular_topology.R — GABAergic Modular Topology Analysis
# =============================================================================
# Research Question (a):
#   How are GABAergic genes distributed across the modular architecture of
#   co-expression networks? Are they confined to specific modules (functional
#   confinement) or dispersed (systemic dysregulation)?
#
# Analyses:
#   1. Module enrichment (Fisher's exact test per module)
#   2. Shannon entropy of modular distribution (vs. permutation null)
#   3. Topological comparison: GABA genes vs. network background
#   4. All comparisons: AD vs. Normal, across 5 brain regions
#
# Inputs:  topology/*_nodes_summary.csv (Louvain communities + metrics)
# Outputs: results/tables/01_*, results/figures/01_*
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# 1. MODULE ENRICHMENT (Fisher's exact test)
# =============================================================================

#' Test enrichment of GABAergic genes in each Louvain module
#' @param topo data.table: topology file (node, degree, pagerank, ..., membership)
#' @param gaba_ids Character vector of Ensembl IDs for GABA genes
#' @param min_module_size Integer: skip modules with fewer total nodes
#' @return data.table with enrichment results per module
compute_module_enrichment <- function(topo, gaba_ids, min_module_size = 10) {
  
  # Total nodes in network
  N <- nrow(topo)
  # GABA genes present in this network
  gaba_in_net <- intersect(gaba_ids, topo$node)
  K <- length(gaba_in_net)
  
  if (K == 0) return(data.table())
  
  # Flag GABA membership
  topo[, is_gaba := node %in% gaba_in_net]
  
  # Get module sizes and GABA counts per module
  mod_stats <- topo[, .(
    module_size   = .N,
    n_gaba        = sum(is_gaba),
    gaba_genes    = paste(node[is_gaba], collapse = ";")
  ), by = membership]
  
  # Filter by minimum size

  mod_stats <- mod_stats[module_size >= min_module_size]
  
  # Fisher's exact test for each module
  mod_stats[, c("odds_ratio", "p_value") := {
    res <- lapply(seq_len(.N), function(i) {
      # 2x2 contingency: GABA vs non-GABA × in-module vs out-of-module
      a <- n_gaba[i]                          # GABA in module
      b <- module_size[i] - n_gaba[i]         # non-GABA in module
      c_val <- K - n_gaba[i]                  # GABA outside module
      d <- (N - module_size[i]) - c_val       # non-GABA outside module
      
      mat <- matrix(c(a, b, c_val, d), nrow = 2)
      ft  <- fisher.test(mat, alternative = "greater")
      list(or = ft$estimate, pv = ft$p.value)
    })
    list(
      sapply(res, `[[`, "or"),
      sapply(res, `[[`, "pv")
    )
  }]
  
  # FDR correction
  mod_stats[, p_adj := p.adjust(p_value, method = "BH")]
  mod_stats[, is_enriched := p_adj < 0.05]
  
  # Clean up
  topo[, is_gaba := NULL]
  
  mod_stats[order(p_value)]
}

# =============================================================================
# 2. SHANNON ENTROPY OF MODULAR DISTRIBUTION
# =============================================================================

#' Compute Shannon entropy of gene distribution across modules
#' @param topo data.table with membership column
#' @param gene_ids Character vector of gene IDs to evaluate
#' @return Numeric: Shannon entropy (bits)
compute_distribution_entropy <- function(topo, gene_ids) {
  genes_in_net <- intersect(gene_ids, topo$node)
  if (length(genes_in_net) < 2) return(NA_real_)
  
  # Frequency of genes across modules
  memberships <- topo[node %in% genes_in_net, membership]
  freq_table  <- table(memberships)
  props       <- as.numeric(freq_table) / sum(freq_table)
  
  # Shannon entropy: H = -sum(p * log2(p))
  -sum(props * log2(props))
}

#' Generate null distribution of entropy by permutation
#' @param topo data.table with membership column
#' @param n_genes Integer: number of genes in the query set
#' @param n_perm Integer: number of permutations
#' @return Numeric vector of null entropies
permute_entropy_null <- function(topo, n_genes, n_perm = 10000) {
  all_nodes <- topo$node
  
  vapply(seq_len(n_perm), function(i) {
    random_genes <- sample(all_nodes, n_genes, replace = FALSE)
    compute_distribution_entropy(topo, random_genes)
  }, numeric(1))
}

#' Test whether observed entropy differs from random expectation
#' @return data.table with: observed_entropy, mean_null, sd_null, z_score, p_value
test_entropy <- function(topo, gene_ids, n_perm = 10000) {
  obs   <- compute_distribution_entropy(topo, gene_ids)
  n_eff <- length(intersect(gene_ids, topo$node))
  
  if (is.na(obs) || n_eff < 2) {
    return(data.table(
      n_genes_in_net = n_eff, observed_entropy = obs,
      mean_null = NA, sd_null = NA, z_score = NA, p_empirical = NA
    ))
  }
  
  null_dist <- permute_entropy_null(topo, n_eff, n_perm)
  
  z <- (obs - mean(null_dist)) / sd(null_dist)
  # Two-sided p-value: is the entropy more extreme than expected?
  p_low  <- mean(null_dist <= obs)  # low entropy = confinement
  p_high <- mean(null_dist >= obs)  # high entropy = dispersion
  p_two  <- 2 * min(p_low, p_high)
  
  data.table(
    n_genes_in_net  = n_eff,
    observed_entropy = obs,
    mean_null       = mean(null_dist),
    sd_null         = sd(null_dist),
    z_score         = z,
    p_empirical     = p_two,
    direction       = ifelse(z < 0, "confined", "dispersed")
  )
}

# =============================================================================
# 3. TOPOLOGICAL COMPARISON: GABA vs. BACKGROUND
# =============================================================================

#' Compare topological metrics of GABA genes vs. all other genes
#' @param topo data.table with node, degree, pagerank, kcore columns
#' @param gaba_ids Character vector of GABA Ensembl IDs
#' @return data.table with Wilcoxon test results for each metric
compare_topology_gaba_vs_background <- function(topo, gaba_ids) {
  
  genes_in_net <- intersect(gaba_ids, topo$node)
  if (length(genes_in_net) < 3) return(data.table())
  
  topo[, group := fifelse(node %in% genes_in_net, "GABA", "Background")]
  
  metrics <- c("degree", "pagerank", "kcore")
  results <- rbindlist(lapply(metrics, function(m) {
    gaba_vals <- topo[group == "GABA", get(m)]
    bg_vals   <- topo[group == "Background", get(m)]
    
    wt <- wilcox.test(gaba_vals, bg_vals, alternative = "two.sided")
    
    data.table(
      metric         = m,
      n_gaba         = length(gaba_vals),
      n_background   = length(bg_vals),
      median_gaba    = median(gaba_vals, na.rm = TRUE),
      median_bg      = median(bg_vals, na.rm = TRUE),
      fold_change    = median(gaba_vals, na.rm = TRUE) / 
                        max(median(bg_vals, na.rm = TRUE), 1e-10),
      W_statistic    = wt$statistic,
      p_value        = wt$p.value
    )
  }))
  
  topo[, group := NULL]
  results
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

#' Bar plot: Number of GABA genes per module, colored by enrichment
plot_module_enrichment <- function(enrich_dt, region, condition, cfg) {
  
  if (nrow(enrich_dt) == 0 || all(enrich_dt$n_gaba == 0)) return(NULL)
  
  plot_data <- enrich_dt[n_gaba > 0]
  plot_data[, module_label := paste0("M", membership)]
  plot_data[, module_label := factor(module_label, 
    levels = plot_data[order(-n_gaba), module_label])]
  
  ggplot(plot_data, aes(x = module_label, y = n_gaba, fill = is_enriched)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
    geom_text(aes(label = n_gaba), vjust = -0.5, size = 3) +
    scale_fill_manual(
      values = c("TRUE" = "#D62828", "FALSE" = "grey70"),
      labels = c("TRUE" = "Enriched (FDR < 0.05)", "FALSE" = "Not enriched"),
      name   = NULL
    ) +
    labs(
      title    = sprintf("%s — %s", region, condition),
      subtitle = sprintf("GABAergic gene distribution across Louvain modules"),
      x = "Module", y = "Number of GABA genes"
    ) +
    theme_publication(cfg = cfg) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Dot plot: GABA vs. Background topological metrics across regions
plot_topology_comparison <- function(topo_comp_dt, cfg) {
  
  topo_comp_dt[, metric := factor(metric, levels = c("degree", "pagerank", "kcore"))]
  
  # Reshape for plotting: we need median_gaba and median_bg side by side
  plot_data <- melt(topo_comp_dt, 
    id.vars = c("region", "condition", "metric"),
    measure.vars = c("median_gaba", "median_bg"),
    variable.name = "group", value.name = "median_value"
  )
  plot_data[, group := fifelse(group == "median_gaba", "GABA", "Background")]
  plot_data[, sig_label := fifelse(
    topo_comp_dt$p_value[match(paste(region, condition, metric), 
      paste(topo_comp_dt$region, topo_comp_dt$condition, topo_comp_dt$metric))] < 0.05,
    "*", ""
  )]
  
  ggplot(plot_data, aes(x = region, y = median_value, color = group, shape = condition)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("GABA" = "#D62828", "Background" = "grey50")) +
    scale_shape_manual(values = c("AD" = 16, "control" = 17)) +
    labs(
      title = "GABAergic vs. Background: Topological Properties",
      x = NULL, y = "Median value", color = NULL, shape = NULL
    ) +
    theme_publication(cfg = cfg)
}

#' Heatmap: Entropy z-scores across regions and conditions
plot_entropy_heatmap <- function(entropy_dt, cfg) {
  
  ggplot(entropy_dt, aes(x = condition, y = region, fill = z_score)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("z=%.1f\np=%.3f", z_score, p_empirical)), 
              size = 3, color = "white") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
      name = "Entropy\nZ-score",
      limits = c(-4, 4), oob = squish
    ) +
    labs(
      title    = "Modular Distribution Entropy of GABAergic Genes",
      subtitle = "Z < 0: Confined to few modules | Z > 0: Dispersed across modules",
      x = NULL, y = NULL
    ) +
    theme_publication(cfg = cfg) +
    theme(panel.grid = element_blank())
}

#' GABA gene profile: subpathway composition per module (for enriched modules)
plot_subpathway_module_composition <- function(topo_all, gene_set, cfg) {
  
  # Get GABA genes with their module assignments
  gaba_nodes <- merge(
    topo_all[node %in% gene_set$ensembl_gene_id, 
             .(node, membership, region, condition)],
    gene_set[, .(ensembl_gene_id, subpathway)],
    by.x = "node", by.y = "ensembl_gene_id"
  )
  
  if (nrow(gaba_nodes) == 0) return(NULL)
  
  # Count subpathway composition per region × condition
  comp <- gaba_nodes[, .N, by = .(region, condition, subpathway)]
  
  ggplot(comp, aes(x = region, y = N, fill = subpathway)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~ condition) +
    scale_fill_manual(values = get_subpathway_colors(cfg), name = "Sub-pathway") +
    labs(
      title = "GABAergic Sub-pathway Representation in Networks",
      x = NULL, y = "Number of genes"
    ) +
    theme_publication(cfg = cfg) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# =============================================================================
# 5. MAIN RUNNER
# =============================================================================

#' Run the complete modular topology analysis
#' @param cfg Config list from load_config()
#' @param gene_set data.table from build_annotated_gene_set()
#' @return Invisible list of all results
run_modular_topology <- function(cfg, gene_set) {
  
  message("\n", strrep("=", 60))
  message("SCRIPT 01: MODULAR TOPOLOGY ANALYSIS")
  message(strrep("=", 60))
  
  gaba_ids <- gene_set[!is.na(ensembl_gene_id), ensembl_gene_id]
  min_mod  <- cfg$analysis$min_module_size
  n_perm   <- cfg$analysis$n_permutations
  
  # ---- A. Module Enrichment -------------------------------------------------
  message("\n--- A. Module Enrichment (Fisher's exact test) ---")
  
  enrichment_fn <- function(cohort, region, condition, cfg) {
    topo <- load_topology(cohort, region, condition, cfg)
    if (is.null(topo)) return(data.table())
    compute_module_enrichment(topo, gaba_ids, min_mod)
  }
  
  enrichment_all <- iterate_networks(cfg, enrichment_fn)
  fwrite(enrichment_all, file.path(cfg$paths$tables_dir, "01_module_enrichment.csv"))
  message("Enrichment results saved.")
  
  # ---- B. Entropy Analysis --------------------------------------------------
  message("\n--- B. Shannon Entropy Analysis ---")
  
  entropy_fn <- function(cohort, region, condition, cfg) {
    topo <- load_topology(cohort, region, condition, cfg)
    if (is.null(topo)) return(data.table())
    test_entropy(topo, gaba_ids, n_perm)
  }
  
  entropy_all <- iterate_networks(cfg, entropy_fn)
  fwrite(entropy_all, file.path(cfg$paths$tables_dir, "01_entropy_analysis.csv"))
  message("Entropy results saved.")
  
  # ---- C. Topological Comparison --------------------------------------------
  message("\n--- C. Topological Properties: GABA vs. Background ---")
  
  topo_comp_fn <- function(cohort, region, condition, cfg) {
    topo <- load_topology(cohort, region, condition, cfg)
    if (is.null(topo)) return(data.table())
    compare_topology_gaba_vs_background(topo, gaba_ids)
  }
  
  topo_comp_all <- iterate_networks(cfg, topo_comp_fn)
  fwrite(topo_comp_all, file.path(cfg$paths$tables_dir, "01_topology_comparison.csv"))
  message("Topology comparison saved.")
  
  # ---- D. Load all topologies for composition plot --------------------------
  message("\n--- D. Collecting topologies for composition plot ---")
  
  all_topos_fn <- function(cohort, region, condition, cfg) {
    topo <- load_topology(cohort, region, condition, cfg)
    if (is.null(topo)) return(data.table())
    topo
  }
  
  all_topos <- iterate_networks(cfg, all_topos_fn)
  
  # ---- E. Generate Plots ----------------------------------------------------
  message("\n--- E. Generating plots ---")
  
  # E1. Module enrichment bar plots (one per region × condition)
  for (reg in unique(enrichment_all$region)) {
    for (cond in unique(enrichment_all$condition)) {
      sub <- enrichment_all[region == reg & condition == cond]
      p <- plot_module_enrichment(sub, reg, cond, cfg)
      if (!is.null(p)) {
        save_plot(p, sprintf("01_enrichment_%s_%s.png", reg, cond), cfg,
                  subdir = "01_modular_topology", width = 8, height = 5)
      }
    }
  }
  
  # E2. Entropy heatmap
  if (nrow(entropy_all) > 0 && any(!is.na(entropy_all$z_score))) {
    p_entropy <- plot_entropy_heatmap(entropy_all, cfg)
    save_plot(p_entropy, "01_entropy_heatmap.png", cfg,
              subdir = "01_modular_topology", width = 6, height = 5)
  }
  
  # E3. Topology comparison dot plot
  if (nrow(topo_comp_all) > 0) {
    p_topo <- plot_topology_comparison(topo_comp_all, cfg)
    save_plot(p_topo, "01_topology_gaba_vs_background.png", cfg,
              subdir = "01_modular_topology", width = 12, height = 5)
  }
  
  # E4. Sub-pathway composition
  if (nrow(all_topos) > 0) {
    p_comp <- plot_subpathway_module_composition(all_topos, gene_set, cfg)
    if (!is.null(p_comp)) {
      save_plot(p_comp, "01_subpathway_composition.png", cfg,
                subdir = "01_modular_topology", width = 10, height = 6)
    }
  }
  
  # ---- F. Summary -----------------------------------------------------------
  message("\n--- Summary ---")
  
  if (nrow(enrichment_all) > 0) {
    n_enriched <- enrichment_all[is_enriched == TRUE, .N]
    n_total_mods <- nrow(enrichment_all)
    message(sprintf("  Enriched modules (FDR<0.05): %d / %d tested", 
                    n_enriched, n_total_mods))
  }
  
  if (nrow(entropy_all) > 0) {
    confined   <- entropy_all[direction == "confined" & p_empirical < 0.05]
    dispersed  <- entropy_all[direction == "dispersed" & p_empirical < 0.05]
    message(sprintf("  Significantly confined: %d networks", nrow(confined)))
    message(sprintf("  Significantly dispersed: %d networks", nrow(dispersed)))
  }
  
  message("\nScript 01 complete.\n")
  
  invisible(list(
    enrichment  = enrichment_all,
    entropy     = entropy_all,
    topo_comp   = topo_comp_all,
    all_topos   = all_topos
  ))
}
