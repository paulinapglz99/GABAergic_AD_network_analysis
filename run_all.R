# =============================================================================
# run_all.R — Controller Script
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("  GABAergic System Reorganization in Alzheimer's Disease\n")
cat(strrep("=", 70), "\n\n")

t0 <- Sys.time()
source("R/00_utils.R")
config <- load_config("config/config.yaml")

message("\n--- Building gene set ---")
gene_set <- build_annotated_gene_set(config)
fwrite(gene_set, file.path(config$paths$gene_sets_dir, "gabaergic_gene_set.csv"))

source("R/01_modular_topology.R")
r01 <- run_modular_topology(config, gene_set)

source("R/02_differential_connectivity.R")
r02 <- run_differential_connectivity(config, gene_set)

source("R/03_cross_regional_conservation.R")
r03 <- run_cross_regional_conservation(config, gene_set)

cat("\n", strrep("=", 70), "\n")
cat(sprintf("  Done in %.1f min. Results: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), config$paths$results_dir))
cat(strrep("=", 70), "\n\n")
