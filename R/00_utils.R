# =============================================================================
# 00_utils.R — Shared utilities (NO biomaRt dependency)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(ggplot2)
  library(scales)
})

# -- Config -------------------------------------------------------------------
load_config <- function(config_path = "config/config.yaml") {
  stopifnot(file.exists(config_path))
  cfg <- yaml::read_yaml(config_path)
  dirs <- c(cfg$paths$results_dir, cfg$paths$figures_dir, cfg$paths$tables_dir,
            cfg$paths$gene_sets_dir)
  lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  cfg
}

# -- Gene set (reads symbol:ENSG mapping directly from config) ----------------
build_annotated_gene_set <- function(cfg) {
  gs <- cfg$gene_set
  dt_list <- lapply(names(gs), function(sp) {
    entries <- gs[[sp]]
    data.table(
      symbol          = names(entries),
      ensembl_gene_id = unlist(entries, use.names = FALSE),
      subpathway      = sp
    )
  })
  gs_table <- rbindlist(dt_list)
  message(sprintf("Gene set: %d genes, %d sub-pathways",
                  nrow(gs_table), uniqueN(gs_table$subpathway)))
  gs_table
}

# -- File path construction ---------------------------------------------------
get_file_path <- function(cohort, region, condition, type = "network", cfg) {
  base <- sprintf("%s_%s_counts_%s_topN200000", cohort, region, condition)
  switch(type,
    network    = file.path(cfg$paths$networks_dir, paste0(base, ".tsv")),
    topology   = file.path(cfg$paths$topology_dir, paste0(base, "_nodes_summary.csv")),
    enrichment = {
      primary  <- file.path(cfg$paths$enrichment_dir, paste0(base, "_enrichment.csv"))
      fallback <- file.path(cfg$paths$enrichment_dir,
                            paste0(base, "_nodes_summary_enrichment.csv"))
      if (file.exists(primary)) primary
      else if (file.exists(fallback)) { message("  (fallback enrichment name)"); fallback }
      else primary
    },
    stop("Unknown type: ", type)
  )
}

# -- Loaders ------------------------------------------------------------------
load_network <- function(cohort, region, condition, cfg) {
  fp <- get_file_path(cohort, region, condition, "network", cfg)
  if (!file.exists(fp)) { warning("Not found: ", fp); return(NULL) }
  fread(fp, sep = "\t", header = TRUE, col.names = c("gene1", "gene2", "MI"))
}

load_topology <- function(cohort, region, condition, cfg) {
  fp <- get_file_path(cohort, region, condition, "topology", cfg)
  if (!file.exists(fp)) { warning("Not found: ", fp); return(NULL) }
  fread(fp, header = TRUE)
}

load_enrichment <- function(cohort, region, condition, cfg) {
  fp <- get_file_path(cohort, region, condition, "enrichment", cfg)
  if (!file.exists(fp)) { warning("Enrichment not found: ", fp); return(NULL) }
  fread(fp, header = TRUE)
}

# -- Iterator -----------------------------------------------------------------
iterate_networks <- function(cfg, fn, ...) {
  results <- list()
  for (ri in cfg$regions) {
    for (cond in cfg$conditions) {
      message(sprintf("  %s_%s_%s", ri$cohort, ri$name, cond))
      res <- fn(ri$cohort, ri$name, cond, cfg, ...)
      if (!is.null(res) && nrow(res) > 0) {
        res[, `:=`(region = ri$name, condition = cond, cohort = ri$cohort)]
        results <- c(results, list(res))
      }
    }
  }
  rbindlist(results, fill = TRUE)
}

# -- Theme & helpers ----------------------------------------------------------
theme_publication <- function(base_size = 12, cfg = NULL) {
  ff <- if (!is.null(cfg)) cfg$plotting$font_family else "sans"
  theme_minimal(base_size = base_size, base_family = ff) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(color = "grey40", size = base_size),
      axis.title = element_text(face = "bold"), legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

get_region_colors    <- function(cfg) unlist(cfg$plotting$region_colors)
get_condition_colors <- function(cfg) unlist(cfg$plotting$condition_colors)
get_subpathway_colors <- function(cfg) unlist(cfg$plotting$subpathway_colors)

save_plot <- function(p, filename, cfg, subdir = NULL, width = NULL, height = NULL) {
  w <- width  %||% cfg$plotting$width
  h <- height %||% cfg$plotting$height
  od <- if (!is.null(subdir)) file.path(cfg$paths$figures_dir, subdir) else cfg$paths$figures_dir
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  fp <- file.path(od, filename)
  ggsave(fp, plot = p, width = w, height = h, dpi = cfg$plotting$dpi, bg = "white")
  message("  Saved: ", fp)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
message("00_utils.R loaded.")
