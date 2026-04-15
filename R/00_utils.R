# =============================================================================
# 00_utils.R — Shared utilities for GABAergic AD Network Analysis
# =============================================================================
# Functions: config loading, gene-set building, ID mapping, I/O, plot theme.
# Source this file before running any analysis script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(ggplot2)
  library(scales)
})

# -- Load configuration -------------------------------------------------------

#' Load and validate the project config.yaml
#' @param config_path Path to config.yaml
#' @return Named list with all config parameters
load_config <- function(config_path = "config/config.yaml") {
  stopifnot(file.exists(config_path))
  cfg <- yaml::read_yaml(config_path)
  
  # Ensure output directories exist
  dirs <- c(cfg$paths$results_dir, cfg$paths$figures_dir, cfg$paths$tables_dir,
            cfg$paths$gene_sets_dir)
  lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  
  cfg
}

# -- Gene set utilities -------------------------------------------------------

#' Flatten the nested gene_set config into a tidy data.table
#' @param cfg Config list (from load_config)
#' @return data.table with columns: symbol, subpathway
build_gene_set_table <- function(cfg) {
  gs <- cfg$gene_set
  dt_list <- lapply(names(gs), function(sp) {
    data.table(symbol = unlist(gs[[sp]]), subpathway = sp)
  })
  rbindlist(dt_list)
}

#' Map HGNC symbols to Ensembl IDs via biomaRt
#' @param symbols Character vector of HGNC symbols
#' @param dataset Ensembl dataset name
#' @return data.table with columns: ensembl_gene_id, hgnc_symbol
map_symbols_to_ensembl <- function(symbols, dataset = "hsapiens_gene_ensembl") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt is required for Ensembl ID mapping. ",
         "Install with: BiocManager::install('biomaRt')")
  }
  
  # Try to use a cached mapping file first
  cache_file <- "data/gene_sets/ensembl_mapping_cache.csv"
  if (file.exists(cache_file)) {
    message("Loading cached Ensembl mapping from: ", cache_file)
    cached <- fread(cache_file)
    # Check if all symbols are in cache
    missing <- setdiff(symbols, cached$hgnc_symbol)
    if (length(missing) == 0) return(cached[hgnc_symbol %in% symbols])
    message(length(missing), " symbols not in cache, querying biomaRt...")
  }
  
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  mapping <- as.data.table(biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "hgnc_symbol",
    values     = symbols,
    mart       = mart
  ))
  
  # Remove duplicates: keep one Ensembl ID per symbol (prefer ENSG on chr 1-22,X,Y)
  mapping <- mapping[hgnc_symbol != ""]
  mapping <- mapping[!duplicated(mapping, by = "hgnc_symbol")]
  
  # Cache for future runs
  fwrite(mapping, cache_file)
  message("Cached Ensembl mapping to: ", cache_file)
  
  mapping
}

#' Build the complete gene set with Ensembl IDs
#' @param cfg Config list
#' @return data.table with columns: ensembl_gene_id, symbol, subpathway
build_annotated_gene_set <- function(cfg) {
  gs_table <- build_gene_set_table(cfg)
  
  if (isTRUE(cfg$use_ensembl)) {
    mapping <- map_symbols_to_ensembl(gs_table$symbol, cfg$ensembl_dataset)
    gs_table <- merge(gs_table, mapping,
                      by.x = "symbol", by.y = "hgnc_symbol", all.x = TRUE)
    n_mapped <- sum(!is.na(gs_table$ensembl_gene_id))
    n_total  <- nrow(gs_table)
    message(sprintf("Gene set: %d / %d symbols mapped to Ensembl IDs", n_mapped, n_total))
    
    if (n_mapped < n_total) {
      unmapped <- gs_table[is.na(ensembl_gene_id), symbol]
      warning("Unmapped symbols: ", paste(unmapped, collapse = ", "))
    }
  }
  
  gs_table
}

# -- Network I/O --------------------------------------------------------------

#' Construct the expected filename for a network or topology file
#' @param cohort Character: "ROSMAP" or "Mayo"
#' @param region Character: brain region abbreviation
#' @param condition Character: "AD" or "control"
#' @param type Character: "network" or "topology"
#' @param cfg Config list
#' @return Character: full file path
get_file_path <- function(cohort, region, condition, type = "network", cfg) {
  base_name <- sprintf("%s_%s_counts_%s_topN200000", cohort, region, condition)
  if (type == "network") {
    file.path(cfg$paths$networks_dir, paste0(base_name, ".tsv"))
  } else if (type == "topology") {
    file.path(cfg$paths$topology_dir, paste0(base_name, "_nodes_summary.csv"))
  } else {
    stop("Unknown type: ", type)
  }
}

#' Load a single edge-list network
#' @return data.table with columns: gene1, gene2, MI
load_network <- function(cohort, region, condition, cfg) {
  fpath <- get_file_path(cohort, region, condition, "network", cfg)
  if (!file.exists(fpath)) {
    warning("Network file not found: ", fpath)
    return(NULL)
  }
  fread(fpath, sep = "\t", header = TRUE, col.names = c("gene1", "gene2", "MI"))
}

#' Load a single topology (node summary) file
#' @return data.table with columns: node, degree, pagerank, pagerank_norm, kcore, membership
load_topology <- function(cohort, region, condition, cfg) {
  fpath <- get_file_path(cohort, region, condition, "topology", cfg)
  if (!file.exists(fpath)) {
    warning("Topology file not found: ", fpath)
    return(NULL)
  }
  fread(fpath, header = TRUE)
}

#' Iterate over all region × condition combinations and apply a function
#' @param cfg Config list
#' @param fn Function taking (cohort, region, condition, cfg, ...) and returning a data.table
#' @param ... Additional arguments passed to fn
#' @return data.table with added columns: region, condition, cohort
iterate_networks <- function(cfg, fn, ...) {
  results <- list()
  for (reg_info in cfg$regions) {
    for (cond in cfg$conditions) {
      message(sprintf("Processing: %s_%s_%s", reg_info$cohort, reg_info$name, cond))
      res <- fn(reg_info$cohort, reg_info$name, cond, cfg, ...)
      if (!is.null(res) && nrow(res) > 0) {
        res[, `:=`(region = reg_info$name, condition = cond, cohort = reg_info$cohort)]
        results <- c(results, list(res))
      }
    }
  }
  rbindlist(results, fill = TRUE)
}

# -- Plotting theme -----------------------------------------------------------

#' Publication-quality ggplot2 theme
#' @param base_size Base font size
#' @param cfg Config list (optional, for font_family)
#' @return ggplot2 theme object
theme_publication <- function(base_size = 12, cfg = NULL) {
  font_fam <- if (!is.null(cfg)) cfg$plotting$font_family else "sans"
  
  theme_minimal(base_size = base_size, base_family = font_fam) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_text(color = "grey40", size = base_size, hjust = 0),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

#' Get region color palette from config
#' @param cfg Config list
#' @return Named character vector
get_region_colors <- function(cfg) {
  unlist(cfg$plotting$region_colors)
}

#' Get condition color palette from config
#' @param cfg Config list
#' @return Named character vector
get_condition_colors <- function(cfg) {
  unlist(cfg$plotting$condition_colors)
}

#' Get subpathway color palette from config
#' @param cfg Config list
#' @return Named character vector
get_subpathway_colors <- function(cfg) {
  unlist(cfg$plotting$subpathway_colors)
}

#' Save a ggplot with standard dimensions
#' @param p ggplot object
#' @param filename Filename (without path)
#' @param cfg Config list
#' @param subdir Subdirectory within figures_dir
save_plot <- function(p, filename, cfg, subdir = NULL, 
                      width = NULL, height = NULL) {
  w <- width  %||% cfg$plotting$width
  h <- height %||% cfg$plotting$height
  
  out_dir <- if (!is.null(subdir)) {
    file.path(cfg$paths$figures_dir, subdir)
  } else {
    cfg$paths$figures_dir
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  fpath <- file.path(out_dir, filename)
  ggsave(fpath, plot = p, width = w, height = h, dpi = cfg$plotting$dpi, bg = "white")
  message("Saved: ", fpath)
}

#' Null-coalescing operator (for R < 4.4 compatibility)
`%||%` <- function(x, y) if (is.null(x)) y else x

message("00_utils.R loaded successfully.")
