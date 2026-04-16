# =============================================================================
# 01_modular_topology.R — GABAergic Modular Topology Analysis
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

# -- Core analyses ------------------------------------------------------------
compute_module_enrichment <- function(topo, gaba_ids, min_mod = 10) {
  N <- nrow(topo); K <- length(intersect(gaba_ids, topo$node))
  if (K == 0) return(data.table())
  topo[, is_gaba := node %in% gaba_ids]
  ms <- topo[, .(mod_size = .N, n_gaba = sum(is_gaba),
                  gaba_genes = paste(node[is_gaba], collapse = ";")), by = membership]
  ms <- ms[mod_size >= min_mod]
  ms[, c("odds_ratio", "p_value") := {
    res <- lapply(seq_len(.N), function(i) {
      mat <- matrix(c(n_gaba[i], mod_size[i] - n_gaba[i],
                       K - n_gaba[i], (N - mod_size[i]) - (K - n_gaba[i])), nrow = 2)
      ft <- fisher.test(mat, alternative = "greater")
      list(unname(ft$estimate), ft$p.value)
    }); list(sapply(res, `[[`, 1), sapply(res, `[[`, 2))
  }]
  ms[, p_adj := p.adjust(p_value, "BH")][, is_enriched := p_adj < 0.05]
  topo[, is_gaba := NULL]
  ms[order(p_value)]
}

compute_entropy <- function(topo, ids) {
  m <- topo[node %in% intersect(ids, topo$node), membership]
  if (length(m) < 2) return(NA_real_)
  p <- as.numeric(table(m)) / length(m); -sum(p * log2(p))
}

test_entropy <- function(topo, ids, nperm = 10000) {
  obs <- compute_entropy(topo, ids); n <- length(intersect(ids, topo$node))
  if (is.na(obs) || n < 2) return(data.table(n_in_net = n, obs_H = obs,
    mean_null = NA, sd_null = NA, z = NA, p = NA, direction = NA_character_))
  nulls <- vapply(seq_len(nperm), function(i)
    compute_entropy(topo, sample(topo$node, n)), numeric(1))
  z <- (obs - mean(nulls)) / sd(nulls)
  data.table(n_in_net = n, obs_H = obs, mean_null = mean(nulls),
    sd_null = sd(nulls), z = z,
    p = 2 * min(mean(nulls <= obs), mean(nulls >= obs)),
    direction = fifelse(z < 0, "confined", "dispersed"))
}

compare_topo <- function(topo, gaba_ids) {
  gin <- intersect(gaba_ids, topo$node)
  if (length(gin) < 3) return(data.table())
  topo[, grp := fifelse(node %in% gin, "GABA", "BG")]
  res <- rbindlist(lapply(c("degree", "pagerank", "kcore"), function(m) {
    g <- topo[grp == "GABA", get(m)]; b <- topo[grp == "BG", get(m)]
    wt <- wilcox.test(g, b)
    data.table(metric = m, n_gaba = length(g), med_gaba = median(g),
               med_bg = median(b), W = unname(wt$statistic), p = wt$p.value)
  }))
  topo[, grp := NULL]; res
}

annotate_modules <- function(co, rg, cd, cfg, mods, top_n = 3) {
  fe <- load_enrichment(co, rg, cd, cfg)
  if (is.null(fe) || nrow(fe) == 0) return(data.table())
  target <- mods[n_gaba > 0, membership]
  fs <- fe[CommunityID %in% target]
  fs[, p.adjust := as.numeric(p.adjust)]
  fs <- fs[p.adjust < cfg$analysis$enrichment_padj_cutoff]
  if (nrow(fs) == 0) return(data.table())
  fs[, ontology := fcase(grepl("^GO:", ID), "GO", grepl("^hsa", ID), "KEGG", default = "Other")]
  fs[order(p.adjust), head(.SD, top_n), by = .(CommunityID, ontology),
     .SDcols = c("ID", "Description", "p.adjust", "Count")]
}

# -- Plots --------------------------------------------------------------------
plot_enrichment_bar <- function(dt, rg, cd, cfg) {
  d <- dt[n_gaba > 0]; if (nrow(d) == 0) return(NULL)
  d[, ml := factor(paste0("M", membership), levels = d[order(-n_gaba), paste0("M", membership)])]
  ggplot(d, aes(ml, n_gaba, fill = is_enriched)) +
    geom_col(width = .7, color = "grey30", linewidth = .3) +
    geom_text(aes(label = n_gaba), vjust = -.5, size = 3) +
    scale_fill_manual(values = c("TRUE" = "#D62828", "FALSE" = "grey70"),
      labels = c("TRUE" = "FDR<0.05", "FALSE" = "NS"), name = NULL) +
    labs(title = sprintf("%s — %s", rg, cd), x = "Module", y = "GABA genes") +
    theme_publication(cfg = cfg) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_entropy_hm <- function(dt, cfg) {
  ggplot(dt, aes(condition, region, fill = z)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("z=%.1f\np=%.3f", z, p)), size = 3, color = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B",
      midpoint = 0, name = "Z-score", limits = c(-4, 4), oob = squish) +
    labs(title = "Entropy of GABAergic Modular Distribution",
         subtitle = "Z<0: confined | Z>0: dispersed", x = NULL, y = NULL) +
    theme_publication(cfg = cfg) + theme(panel.grid = element_blank())
}

plot_topo_dot <- function(dt, cfg) {
  d <- melt(dt, id.vars = c("region", "condition", "metric"),
    measure.vars = c("med_gaba", "med_bg"), variable.name = "grp", value.name = "val")
  d[, grp := fifelse(grp == "med_gaba", "GABA", "Background")]
  ggplot(d, aes(region, val, color = grp, shape = condition)) +
    geom_point(size = 3, position = position_dodge(.5)) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c(GABA = "#D62828", Background = "grey50")) +
    scale_shape_manual(values = c(AD = 16, control = 17)) +
    labs(title = "GABA vs Background Topology", x = NULL, y = "Median", color = NULL, shape = NULL) +
    theme_publication(cfg = cfg)
}

plot_subpathway <- function(topos, gs, cfg) {
  d <- merge(topos[node %in% gs$ensembl_gene_id, .(node, region, condition)],
             gs[, .(ensembl_gene_id, subpathway)], by.x = "node", by.y = "ensembl_gene_id")
  if (nrow(d) == 0) return(NULL)
  ct <- d[, .N, by = .(region, condition, subpathway)]
  ggplot(ct, aes(region, N, fill = subpathway)) + geom_col(position = "stack", width = .7) +
    facet_wrap(~ condition) + scale_fill_manual(values = get_subpathway_colors(cfg)) +
    labs(title = "Sub-pathway Representation", x = NULL, y = "Genes") +
    theme_publication(cfg = cfg) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -- Main runner --------------------------------------------------------------
run_modular_topology <- function(cfg, gene_set) {
  message("\n", strrep("=", 60), "\nSCRIPT 01: MODULAR TOPOLOGY\n", strrep("=", 60))
  gids <- gene_set$ensembl_gene_id; mm <- cfg$analysis$min_module_size
  np <- cfg$analysis$n_permutations; tn <- cfg$analysis$top_n_terms; sd <- "01_modular_topology"

  message("\n--- A. Module Enrichment ---")
  enr <- iterate_networks(cfg, function(co, rg, cd, cfg)
    { t <- load_topology(co, rg, cd, cfg); if (is.null(t)) data.table() else compute_module_enrichment(t, gids, mm) })
  fwrite(enr, file.path(cfg$paths$tables_dir, "01_module_enrichment.csv"))

  message("\n--- B. Entropy ---")
  ent <- iterate_networks(cfg, function(co, rg, cd, cfg)
    { t <- load_topology(co, rg, cd, cfg); if (is.null(t)) data.table() else test_entropy(t, gids, np) })
  fwrite(ent, file.path(cfg$paths$tables_dir, "01_entropy.csv"))

  message("\n--- C. Topology GABA vs BG ---")
  tc <- iterate_networks(cfg, function(co, rg, cd, cfg)
    { t <- load_topology(co, rg, cd, cfg); if (is.null(t)) data.table() else compare_topo(t, gids) })
  fwrite(tc, file.path(cfg$paths$tables_dir, "01_topology_comparison.csv"))

  message("\n--- D. Functional Annotation ---")
  fa <- iterate_networks(cfg, function(co, rg, cd, cfg) {
    t <- load_topology(co, rg, cd, cfg); if (is.null(t)) return(data.table())
    me <- compute_module_enrichment(t, gids, mm)
    if (nrow(me) == 0) return(data.table())
    annotate_modules(co, rg, cd, cfg, me, tn) })
  if (nrow(fa) > 0) fwrite(fa, file.path(cfg$paths$tables_dir, "01_functional_annotation.csv"))

  message("\n--- E. Topologies ---")
  at <- iterate_networks(cfg, function(co, rg, cd, cfg)
    { t <- load_topology(co, rg, cd, cfg); if (is.null(t)) data.table() else t })

  message("\n--- F. Plots ---")
  for (r in unique(enr$region)) for (c in unique(enr$condition)) {
    p <- plot_enrichment_bar(enr[region == r & condition == c], r, c, cfg)
    if (!is.null(p)) save_plot(p, sprintf("01_enrich_%s_%s.png", r, c), cfg, sd, 8, 5)
  }
  if (nrow(ent) > 0 && any(!is.na(ent$z))) save_plot(plot_entropy_hm(ent, cfg), "01_entropy.png", cfg, sd, 6, 5)
  if (nrow(tc) > 0) save_plot(plot_topo_dot(tc, cfg), "01_topo_comparison.png", cfg, sd, 12, 5)
  if (nrow(at) > 0) { p <- plot_subpathway(at, gene_set, cfg); if (!is.null(p)) save_plot(p, "01_subpathway.png", cfg, sd, 10, 6) }

  message("\n--- G. Integrated table ---")
  if (nrow(enr) > 0 && nrow(fa) > 0) {
    fs <- fa[, .(top_terms = paste(Description, collapse = " | ")),
             by = .(region, condition, CommunityID)]
    setnames(fs, "CommunityID", "membership")
    intg <- merge(enr[n_gaba > 0, .(region, condition, membership, mod_size, n_gaba, odds_ratio, p_adj, is_enriched)],
                  fs, by = c("region", "condition", "membership"), all.x = TRUE)
    intg[is.na(top_terms), top_terms := "—"]
    fwrite(intg, file.path(cfg$paths$tables_dir, "01_integrated_summary.csv"))
  }

  message("\n--- Summary ---")
  if (nrow(enr) > 0) message(sprintf("  Enriched modules: %d / %d", sum(enr$is_enriched, na.rm = TRUE), nrow(enr)))
  if (nrow(ent) > 0 && any(!is.na(ent$z))) {
    message(sprintf("  Confined: %d | Dispersed: %d", nrow(ent[direction == "confined" & p < .05]),
                    nrow(ent[direction == "dispersed" & p < .05])))
  }
  message("\nScript 01 done.\n")
  invisible(list(enrichment = enr, entropy = ent, topo_comp = tc, func = fa, topos = at))
}
