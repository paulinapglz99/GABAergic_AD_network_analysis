# =============================================================================
# 02_differential_connectivity.R — Differential Connectivity Analysis
# =============================================================================
# Question (b): Which GABAergic genes undergo the most drastic connectivity
# changes between Normal and AD? To which sub-pathways do they belong?
#
# Analyses:
#   1. Per-gene degree in AD vs Normal → delta_degree, ratio
#   2. Per-gene MI-weighted strength in AD vs Normal
#   3. Neighborhood rewiring: Jaccard index of neighbor sets
#   4. Sub-pathway annotation of top DC genes
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# -- Core: compute degree/strength for GABA genes in one network -------------
compute_gaba_connectivity <- function(cohort, region, condition, cfg, gaba_ids) {
  net <- load_network(cohort, region, condition, cfg)
  if (is.null(net)) return(data.table())
  
  # Filter edges touching at least one GABA gene
  gaba_edges <- net[gene1 %in% gaba_ids | gene2 %in% gaba_ids]
  if (nrow(gaba_edges) == 0) return(data.table())
  
  # For each GABA gene: degree (n edges) and strength (sum MI)
  e1 <- gaba_edges[gene1 %in% gaba_ids, .(degree = .N, strength = sum(MI)), by = .(gene = gene1)]
  e2 <- gaba_edges[gene2 %in% gaba_ids, .(degree = .N, strength = sum(MI)), by = .(gene = gene2)]
  conn <- rbind(e1, e2)[, .(degree = sum(degree), strength = sum(strength)), by = gene]
  conn
}

# -- Core: compute neighbor set for each GABA gene ---------------------------
get_neighbor_sets <- function(cohort, region, condition, cfg, gaba_ids) {
  net <- load_network(cohort, region, condition, cfg)
  if (is.null(net)) return(list())
  
  gaba_edges <- net[gene1 %in% gaba_ids | gene2 %in% gaba_ids]
  if (nrow(gaba_edges) == 0) return(list())
  
  # Build neighbor list
  neighbors <- list()
  for (g in gaba_ids) {
    n1 <- gaba_edges[gene1 == g, gene2]
    n2 <- gaba_edges[gene2 == g, gene1]
    nb <- unique(c(n1, n2))
    if (length(nb) > 0) neighbors[[g]] <- nb
  }
  neighbors
}

# -- Jaccard index between two sets -------------------------------------------
jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# -- Differential connectivity for one region ---------------------------------
compute_dc_one_region <- function(cohort, region, cfg, gaba_ids) {
  
  conn_ad   <- compute_gaba_connectivity(cohort, region, "AD", cfg, gaba_ids)
  conn_ctrl <- compute_gaba_connectivity(cohort, region, "control", cfg, gaba_ids)
  
  if (nrow(conn_ad) == 0 || nrow(conn_ctrl) == 0) return(data.table())
  
  # Merge AD and control
  dc <- merge(conn_ad, conn_ctrl, by = "gene", suffixes = c("_AD", "_ctrl"), all = TRUE)
  
  # Fill NAs with 0 (gene present in one condition but not the other)
  for (col in names(dc)[-1]) set(dc, which(is.na(dc[[col]])), col, 0)
  
  # Delta and ratio
  dc[, delta_degree   := degree_AD - degree_ctrl]
  dc[, delta_strength := strength_AD - strength_ctrl]
  dc[, ratio_degree   := fifelse(degree_ctrl > 0, degree_AD / degree_ctrl, NA_real_)]
  
  # Neighbor rewiring (Jaccard)
  nb_ad   <- get_neighbor_sets(cohort, region, "AD", cfg, gaba_ids)
  nb_ctrl <- get_neighbor_sets(cohort, region, "control", cfg, gaba_ids)
  
  dc[, jaccard_neighbors := sapply(gene, function(g) {
    jaccard(nb_ad[[g]], nb_ctrl[[g]])
  })]
  
  # Rewiring score: 1 - Jaccard (high = more rewiring)
  dc[, rewiring_score := 1 - jaccard_neighbors]
  
  dc
}

# -- Plots --------------------------------------------------------------------
plot_dc_volcano <- function(dc_all, gene_set, cfg) {
  # Merge with gene annotations
  d <- merge(dc_all, gene_set[, .(ensembl_gene_id, symbol, subpathway)],
             by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)
  d[is.na(symbol), symbol := substr(gene, 1, 15)]
  
  # Top movers: top 5 absolute delta_degree per region
  d[, rank_dc := frank(-abs(delta_degree)), by = region]
  d[, is_top := rank_dc <= 5]
  
  ggplot(d, aes(x = delta_degree, y = rewiring_score, color = subpathway)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_text(data = d[is_top == TRUE], aes(label = symbol),
              size = 2.5, vjust = -0.8, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ region, scales = "free_x", ncol = 3) +
    scale_color_manual(values = get_subpathway_colors(cfg), name = "Sub-pathway") +
    labs(
      title = "Differential Connectivity of GABAergic Genes (AD vs Normal)",
      subtitle = "X: change in degree | Y: rewiring score (1 - Jaccard)",
      x = "Δ Degree (AD − Normal)", y = "Rewiring Score"
    ) +
    theme_publication(cfg = cfg)
}

plot_dc_heatmap <- function(dc_all, gene_set, cfg) {
  d <- merge(dc_all, gene_set[, .(ensembl_gene_id, symbol, subpathway)],
             by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)
  d[is.na(symbol), symbol := substr(gene, 1, 15)]
  
  # Keep genes present in at least 3 regions
  gene_counts <- d[, .N, by = symbol]
  keep <- gene_counts[N >= 3, symbol]
  d <- d[symbol %in% keep]
  
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = region, y = reorder(symbol, delta_degree), fill = delta_degree)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, name = "Δ Degree") +
    labs(title = "Degree Change Heatmap: GABA Genes × Regions",
         subtitle = "Blue: lost connectivity in AD | Red: gained connectivity",
         x = NULL, y = NULL) +
    theme_publication(cfg = cfg, base_size = 10) +
    theme(axis.text.y = element_text(size = 7))
}

plot_dc_subpathway_summary <- function(dc_all, gene_set, cfg) {
  d <- merge(dc_all, gene_set[, .(ensembl_gene_id, subpathway)],
             by.x = "gene", by.y = "ensembl_gene_id")
  
  sp_summary <- d[, .(
    mean_delta = mean(delta_degree, na.rm = TRUE),
    mean_rewiring = mean(rewiring_score, na.rm = TRUE),
    n_genes = .N
  ), by = .(region, subpathway)]
  
  ggplot(sp_summary, aes(x = region, y = subpathway, size = n_genes, color = mean_delta)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "grey80", high = "#B2182B",
                          midpoint = 0, name = "Mean Δ Degree") +
    scale_size_continuous(range = c(2, 8), name = "Genes") +
    labs(title = "Sub-pathway Connectivity Changes Across Regions",
         x = NULL, y = NULL) +
    theme_publication(cfg = cfg)
}

# -- Main runner --------------------------------------------------------------
run_differential_connectivity <- function(cfg, gene_set) {
  message("\n", strrep("=", 60), "\nSCRIPT 02: DIFFERENTIAL CONNECTIVITY\n", strrep("=", 60))
  
  gids <- gene_set$ensembl_gene_id
  sd <- "02_differential_connectivity"
  
  message("\n--- Computing DC per region ---")
  dc_list <- list()
  for (ri in cfg$regions) {
    message(sprintf("  %s_%s", ri$cohort, ri$name))
    dc <- compute_dc_one_region(ri$cohort, ri$name, cfg, gids)
    if (nrow(dc) > 0) {
      dc[, `:=`(region = ri$name, cohort = ri$cohort)]
      dc_list <- c(dc_list, list(dc))
    }
  }
  dc_all <- rbindlist(dc_list, fill = TRUE)
  
  if (nrow(dc_all) == 0) {
    message("  No DC results. Check that both AD and control networks exist.")
    message("Script 02 done.\n")
    return(invisible(list()))
  }
  
  fwrite(dc_all, file.path(cfg$paths$tables_dir, "02_differential_connectivity.csv"))
  message("  DC table saved.")
  
  # Top DC genes table
  dc_ann <- merge(dc_all, gene_set[, .(ensembl_gene_id, symbol, subpathway)],
                  by.x = "gene", by.y = "ensembl_gene_id", all.x = TRUE)
  top_dc <- dc_ann[order(-abs(delta_degree)), head(.SD, 10), by = region]
  fwrite(top_dc, file.path(cfg$paths$tables_dir, "02_top_dc_genes.csv"))
  
  message("\n--- Plots ---")
  p1 <- plot_dc_volcano(dc_all, gene_set, cfg)
  save_plot(p1, "02_dc_volcano.png", cfg, sd, 14, 8)
  
  p2 <- plot_dc_heatmap(dc_all, gene_set, cfg)
  if (!is.null(p2)) save_plot(p2, "02_dc_heatmap.png", cfg, sd, 8, 10)
  
  p3 <- plot_dc_subpathway_summary(dc_all, gene_set, cfg)
  save_plot(p3, "02_dc_subpathway.png", cfg, sd, 10, 7)
  
  message("\n--- Summary ---")
  msg <- dc_ann[, .(mean_delta = mean(delta_degree, na.rm = TRUE),
                     mean_rewire = mean(rewiring_score, na.rm = TRUE)), by = region]
  for (i in seq_len(nrow(msg)))
    message(sprintf("  %s: mean Δdegree=%.1f, mean rewiring=%.2f", msg$region[i], msg$mean_delta[i], msg$mean_rewire[i]))
  
  message("\nScript 02 done.\n")
  invisible(list(dc = dc_all, top = top_dc))
}
