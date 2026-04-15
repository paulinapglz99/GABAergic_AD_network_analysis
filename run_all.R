# =============================================================================
# run_all.R — Controller Script for GABAergic AD Network Analysis
# =============================================================================
# This script orchestrates the full analysis pipeline.
# Run from the project root: source("run_all.R")
#
# Workflow:
#   1. Load config and shared utilities
#   2. Build annotated GABAergic gene set (with Ensembl mapping)
#   3. Run Script 01: Modular Topology
#   4. Run Script 02: Differential Connectivity
#   5. Run Script 03: Cross-Regional Conservation
#
# Each script can also be run independently after sourcing 00_utils.R.
# =============================================================================

# -- Setup --------------------------------------------------------------------
cat("\n")
cat(strrep("=", 70), "\n")
cat("  GABAergic System Reorganization in Alzheimer's Disease\n")
cat("  A Network Medicine Approach\n")
cat(strrep("=", 70), "\n\n")

t_start <- Sys.time()

# Source shared utilities
source("R/00_utils.R")

# Load configuration
config <- load_config("config/config.yaml")

# -- Build Gene Set -----------------------------------------------------------
message("\n--- Building annotated GABAergic gene set ---")
gene_set <- build_annotated_gene_set(config)

# Save gene set for reference
gs_path <- file.path(config$paths$gene_sets_dir, "gabaergic_gene_set.csv")
fwrite(gene_set, gs_path)
message("Gene set saved to: ", gs_path)
message(sprintf("Total genes: %d | Sub-pathways: %d", 
                nrow(gene_set), uniqueN(gene_set$subpathway)))

# -- Script 01: Modular Topology ----------------------------------------------
source("R/01_modular_topology.R")
results_01 <- run_modular_topology(config, gene_set)

# -- Script 02: Differential Connectivity -------------------------------------
source("R/02_differential_connectivity.R")
results_02 <- run_differential_connectivity(config, gene_set)

# -- Script 03: Cross-Regional Conservation -----------------------------------
source("R/03_cross_regional_conservation.R")
results_03 <- run_cross_regional_conservation(config, gene_set)

# -- Done ---------------------------------------------------------------------
t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "mins")

cat("\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Pipeline complete. Total runtime: %.1f minutes.\n", as.numeric(elapsed)))
cat(sprintf("  Results in: %s\n", config$paths$results_dir))
cat(strrep("=", 70), "\n\n")
