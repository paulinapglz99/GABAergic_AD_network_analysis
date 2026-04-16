# GABAergic System Reorganization in Alzheimer's Disease

## Setup
```bash
# Copy data
cp ~/path/to/fomo_networks/*.tsv data/networks/
cp ~/path/to/results_topos_louvain/*_nodes_summary.csv data/topology/
cp ~/path/to/results_topos_louvain/results_comm/*_enrichment.csv data/enrichment/

# Run
Rscript run_all.R
```

## Requirements
R >= 4.1: `data.table`, `ggplot2`, `yaml`, `patchwork`, `scales`
No Bioconductor dependencies. Ensembl IDs are hardcoded in `config.yaml`.

## Scripts
| Script | Question | Key outputs |
|--------|----------|-------------|
| `01_modular_topology.R` | Are GABA genes confined or dispersed? | Enrichment, entropy, topology |
| `02_differential_connectivity.R` | Which GABA genes rewire in AD? | DC scores, Jaccard rewiring |
| `03_cross_regional_conservation.R` | Is modular structure conserved across regions? | Concordance, stability |

## Gene Set
40 genes, 7 sub-pathways. Curated from Carello-Collar 2023, Jiménez-Balado 2021, Krueger-Burg 2025, Marilovtseva 2025.
