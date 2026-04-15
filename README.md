# GABAergic System Reorganization in Alzheimer's Disease
## A Network Medicine Approach

### Overview
This project interrogates co-expression networks from 5 brain regions (DLPFC, HCN, PCC, TC, CRB) 
across two conditions (Aβ pathology vs. Normal aging) to characterize the reorganization of the 
GABAergic system in Alzheimer's Disease.

### Data Sources
- **ROSMAP**: DLPFC, HCN, PCC (Religious Orders Study / Memory and Aging Project)
- **Mayo**: TC, CRB (Mayo Clinic Brain Bank)
- Networks are mutual information-based co-expression edge-lists (top 200K edges)
- Community detection: Louvain algorithm

### Research Questions
1. **Modular Topology** (`01_modular_topology.R`): Are GABAergic genes confined to specific 
   functional modules or dispersed across the network architecture?
2. **Differential Connectivity** (`02_differential_connectivity.R`): Which GABAergic genes undergo 
   the most drastic connectivity rewiring between Normal and AD?
3. **Cross-Regional Conservation** (`03_cross_regional_conservation.R`): Is the modular organization 
   of GABAergic genes conserved across brain regions, and does AD disrupt this conservation?

### GABAergic Gene Set
Curated from five published reviews/meta-analyses, organized into functional sub-pathways:
- **Synthesis**: GAD1, GAD2
- **Transport**: SLC6A1, SLC6A11, SLC6A13, SLC32A1
- **GABA-A Receptors**: GABRA1-6, GABRB1-3, GABRG1-3, GABRD, GABRE, GABRP, GABRQ
- **GABA-B Receptors**: GABBR1, GABBR2
- **Catabolism**: ABAT, ALDH5A1
- **Interneuron Markers**: PVALB, SST, VIP, CCK, CALB1, CALB2, NPY
- **Synapse Organizers**: NLGN2, NLGN3, NRXN1, NRXN3, SLITRK3, GPHN, ARHGEF9

### Repository Structure
```
GABAergic_AD_network_analysis/
├── README.md
├── config/
│   └── config.yaml
├── R/
│   ├── 00_utils.R
│   ├── 01_modular_topology.R
│   ├── 02_differential_connectivity.R
│   └── 03_cross_regional_conservation.R
├── run_all.R
├── data/
│   ├── networks/          # .tsv edge-lists (user-provided)
│   ├── topology/          # _nodes_summary.csv (user-provided)
│   └── gene_sets/         # auto-generated from config
└── results/
    ├── figures/
    └── tables/
```

### Usage
```r
# From the project root:
source("run_all.R")

# Or run individual analyses:
source("R/00_utils.R")
source("R/01_modular_topology.R")
run_modular_topology(config)
```

### Requirements
R >= 4.1, with packages: data.table, igraph, ggplot2, yaml, patchwork, scales, biomaRt (optional)

### References
- Carello-Collar et al. (2023) Mol Psychiatry 28:5025-5036
- Jiménez-Balado & Eich (2021) Semin Cell Dev Biol 116:146-159
- Krueger-Burg (2025) Trends Neurosci 48(1):47-62
- Marilovtseva et al. (2025) Int J Mol Sci 26:9492
