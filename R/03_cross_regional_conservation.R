# =============================================================================
# 03_cross_regional_conservation.R — Cross-Regional Modular Conservation
# =============================================================================
# Question (c): Is the modular co-assignment of GABAergic genes conserved
# across brain regions in Normal, and does Aβ pathology disrupt this?
#
# Analyses:
#   1. Co-assignment matrix: do gene X and Y share a module in region R?
#   2. Cross-regional concordance (Jaccard of co-assignment vectors)
#   3. AD vs Normal comparison of conservation scores
#   4. Identification of region-specific vs. conserved GABA modules
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# -- Build co-assignment vector for GABA genes in one network -----------------
#' For each pair of GABA genes present in the network, return 1 if they share
#' the same Louvain module, 0 otherwise. Returns a named vector keyed by
#' "geneA__geneB" (alphabetically sorted pair).
build_coassignment <- function(cohort, region, condition, cfg, gaba_ids) {
  topo <- load_topology(cohort, region, condition, cfg)
  if (is.null(topo)) return(NULL)
  
  present <- intersect(gaba_ids, topo$node)
  if (length(present) < 2) return(NULL)
  
  # Get memberships
  memb <- topo[node %in% present, .(node, membership)]
  setkey(memb, node)
  
  # All pairs (upper triangle)
  pairs <- combn(sort(present), 2)
  pair_ids <- paste(pairs[1,], pairs[2,], sep = "__")
  
  co_assign <- as.integer(
    memb[pairs[1,], membership] == memb[pairs[2,], membership]
  )
  names(co_assign) <- pair_ids
  co_assign
}

# -- Jaccard between two co-assignment vectors --------------------------------
coassign_jaccard <- function(v1, v2) {
  shared <- intersect(names(v1), names(v2))
  if (length(shared) == 0) return(NA_real_)
  a <- v1[shared]; b <- v2[shared]
  # Jaccard: proportion of pairs where both agree they are co-assigned
  sum(a == 1 & b == 1) / max(sum(a == 1 | b == 1), 1)
}

# -- Compute all pairwise region concordances for one condition ---------------
compute_regional_concordance <- function(cfg, gaba_ids, condition) {
  # Build co-assignment for each region
  ca_list <- list()
  for (ri in cfg$regions) {
    ca <- build_coassignment(ri$cohort, ri$name, condition, cfg, gaba_ids)
    if (!is.null(ca)) ca_list[[ri$name]] <- ca
  }
  
  if (length(ca_list) < 2) return(data.table())
  
  regions <- names(ca_list)
  pairs <- combn(regions, 2)
  
  rbindlist(lapply(seq_len(ncol(pairs)), function(i) {
    r1 <- pairs[1, i]; r2 <- pairs[2, i]
    j <- coassign_jaccard(ca_list[[r1]], ca_list[[r2]])
    n_shared <- length(intersect(names(ca_list[[r1]]), names(ca_list[[r2]])))
    data.table(region1 = r1, region2 = r2, jaccard = j,
               n_shared_pairs = n_shared, condition = condition)
  }))
}

# -- Per-gene module stability across regions ---------------------------------
compute_gene_module_stability <- function(cfg, gaba_ids, gene_set) {
  # For each gene, collect its module membership across regions
  membership_list <- list()
  for (cond in cfg$conditions) {
    for (ri in cfg$regions) {
      topo <- load_topology(ri$cohort, ri$name, cond, cfg)
      if (is.null(topo)) next
      present <- topo[node %in% gaba_ids, .(node, membership)]
      if (nrow(present) > 0) {
        present[, `:=`(region = ri$name, condition = cond)]
        membership_list <- c(membership_list, list(present))
      }
    }
  }
  
  if (length(membership_list) == 0) return(data.table())
  memb_all <- rbindlist(membership_list)
  
  # For each gene × condition: how many DISTINCT modules across regions?
  stability <- memb_all[, .(
    n_regions   = .N,
    n_modules   = uniqueN(membership),
    module_list = paste(unique(membership), collapse = ",")
  ), by = .(node, condition)]
  
  # Stability score: 1 if always same module, decreases with more modules
  stability[, stability_score := 1 / n_modules]
  
  # Annotate
  stability <- merge(stability, gene_set[, .(ensembl_gene_id, symbol, subpathway)],
                     by.x = "node", by.y = "ensembl_gene_id", all.x = TRUE)
  stability
}

# -- Plots --------------------------------------------------------------------
plot_concordance_heatmap <- function(conc_dt, cfg) {
  if (nrow(conc_dt) == 0) return(NULL)
  
  # Make symmetric for heatmap
  d <- rbind(conc_dt[, .(r1 = region1, r2 = region2, jaccard, condition)],
             conc_dt[, .(r1 = region2, r2 = region1, jaccard, condition)])
  
  ggplot(d, aes(r1, r2, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
    facet_wrap(~ condition) +
    scale_fill_gradient(low = "white", high = "#2A9D8F", name = "Jaccard",
                        limits = c(0, 1)) +
    labs(title = "Cross-Regional Co-Assignment Concordance",
         subtitle = "Higher Jaccard = GABA genes share modules across regions",
         x = NULL, y = NULL) +
    theme_publication(cfg = cfg) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_concordance_delta <- function(conc_dt, cfg) {
  if (nrow(conc_dt) == 0) return(NULL)
  
  # Compute AD − Control delta per region pair
  wide <- dcast(conc_dt, region1 + region2 ~ condition, value.var = "jaccard")
  if (!"AD" %in% names(wide) || !"control" %in% names(wide)) return(NULL)
  wide[, delta := AD - control]
  wide[, pair := paste(region1, region2, sep = " – ")]
  
  ggplot(wide, aes(x = reorder(pair, delta), y = delta, fill = delta > 0)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "#D62828", "FALSE" = "#2166AC"),
                      labels = c("TRUE" = "More conserved in AD",
                                 "FALSE" = "Less conserved in AD"), name = NULL) +
    coord_flip() +
    labs(title = "Change in Cross-Regional Conservation (AD − Normal)",
         subtitle = "Negative = AD disrupts shared modular structure",
         x = NULL, y = "Δ Jaccard (AD − Normal)") +
    theme_publication(cfg = cfg)
}

plot_gene_stability <- function(stab_dt, cfg) {
  if (nrow(stab_dt) == 0) return(NULL)
  
  d <- stab_dt[n_regions >= 3]  # Only genes present in 3+ regions
  if (nrow(d) == 0) return(NULL)
  
  # Compare stability between conditions
  wide <- dcast(d, symbol + subpathway ~ condition, value.var = "stability_score")
  if (!"AD" %in% names(wide) || !"control" %in% names(wide)) return(NULL)
  wide[, delta_stability := AD - control]
  wide <- wide[!is.na(delta_stability)]
  
  ggplot(wide, aes(x = reorder(symbol, delta_stability), y = delta_stability,
                    fill = subpathway)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = get_subpathway_colors(cfg), name = "Sub-pathway") +
    coord_flip() +
    labs(title = "Gene-Level Modular Stability Change (AD − Normal)",
         subtitle = "Positive = more stable in AD | Negative = fragmented in AD",
         x = NULL, y = "Δ Stability Score") +
    theme_publication(cfg = cfg, base_size = 10) +
    theme(axis.text.y = element_text(size = 7))
}

# -- Main runner --------------------------------------------------------------
run_cross_regional_conservation <- function(cfg, gene_set) {
  message("\n", strrep("=", 60), "\nSCRIPT 03: CROSS-REGIONAL CONSERVATION\n", strrep("=", 60))
  
  gids <- gene_set$ensembl_gene_id
  sd <- "03_cross_regional"
  
  # A. Regional concordance
  message("\n--- A. Cross-regional concordance ---")
  conc_list <- list()
  for (cond in cfg$conditions) {
    message(sprintf("  Condition: %s", cond))
    conc <- compute_regional_concordance(cfg, gids, cond)
    if (nrow(conc) > 0) conc_list <- c(conc_list, list(conc))
  }
  conc_all <- rbindlist(conc_list, fill = TRUE)
  
  if (nrow(conc_all) > 0) {
    fwrite(conc_all, file.path(cfg$paths$tables_dir, "03_regional_concordance.csv"))
    message("  Concordance table saved.")
  }
  
  # B. Gene-level stability
  message("\n--- B. Gene-level modular stability ---")
  stab <- compute_gene_module_stability(cfg, gids, gene_set)
  if (nrow(stab) > 0) {
    fwrite(stab, file.path(cfg$paths$tables_dir, "03_gene_stability.csv"))
    message("  Stability table saved.")
  }
  
  # C. Plots
  message("\n--- C. Plots ---")
  if (nrow(conc_all) > 0) {
    p1 <- plot_concordance_heatmap(conc_all, cfg)
    if (!is.null(p1)) save_plot(p1, "03_concordance_heatmap.png", cfg, sd, 10, 6)
    
    p2 <- plot_concordance_delta(conc_all, cfg)
    if (!is.null(p2)) save_plot(p2, "03_concordance_delta.png", cfg, sd, 8, 6)
  }
  
  if (nrow(stab) > 0) {
    p3 <- plot_gene_stability(stab, cfg)
    if (!is.null(p3)) save_plot(p3, "03_gene_stability.png", cfg, sd, 10, 10)
  }
  
  # D. Summary
  message("\n--- Summary ---")
  if (nrow(conc_all) > 0) {
    for (cond in unique(conc_all$condition)) {
      mj <- conc_all[condition == cond, mean(jaccard, na.rm = TRUE)]
      message(sprintf("  %s: mean cross-regional Jaccard = %.3f", cond, mj))
    }
  }
  if (nrow(stab) > 0) {
    for (cond in unique(stab$condition)) {
      ms <- stab[condition == cond, mean(stability_score, na.rm = TRUE)]
      message(sprintf("  %s: mean gene stability = %.3f", cond, ms))
    }
  }
  
  message("\nScript 03 done.\n")
  invisible(list(concordance = conc_all, stability = stab))
}
