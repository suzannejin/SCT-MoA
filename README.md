# Evaluating Measures of Association for Single-Cell Transcriptomics

This repository contains the code and data to reproduce the analyses presented in the paper *"Evaluating measures of association for single-cell transcriptomics."*

The project benchmarks **17 measures of association** (correlation, distance, and proportionality metrics) across **211 scRNA-seq datasets** to evaluate their utility for coexpression network inference, functional coherence, network overlap with known biology, cell clustering, reproducibility, and disease gene prediction.

---

## Codebase Structure

```
SCT-MoA/
│
├── main.nf                      # Nextflow pipeline entrypoint (DSL2)
├── nextflow.config              # Pipeline parameters and profiles
├── conf/                        # Nextflow profile configs
│   ├── base.config              #   Default resource limits (CPU, memory, time)
│   ├── crg.config               #   CRG HPC cluster settings
│   ├── test.config              #   Test profile
│   └── trace.config             #   Execution tracing
│
├── modules/                     # Nextflow process definitions
│   ├── coexpr.nf                #   Coexpression matrix generation & filtering
│   ├── evaluation.nf            #   EGAD (AUROC) & network overlap evaluation
│   └── network.nf               #   Reference network rewiring
│
├── bin/                         # Scripts called by Nextflow processes & other scripts I used
│   ├── functions.R              #   Shared utilities (gene mapping, species detection)
│   ├── coexpr/
│   │   ├── write-matrix.R       #     Compute coexpression matrix (dismay/propr)
│   │   └── filter-matrix.R      #     Filter genes by expression prevalence
│   ├── egad/
│   │   ├── calculate-auroc.R    #     AUROC for GO/Reactome annotations
│   │   ├── consolidate-auroc.R  #     Merge AUROC results across datasets
│   │   └── plot-auroc*.R        #     Visualization scripts
│   └── overlap/
│       ├── rewire-networks.R    #     Generate 1000 rewired null networks
│       ├── calculate-overlap.R  #     Compute overlap statistics
│       ├── consolidate-overlap.R#     Merge overlap results
│       └── plot-network-overlap*.R #  Visualization scripts
│
├── R/                           # Legacy code - from Skinnider
│   ├── functions.R              #   Shared R utilities
│   ├── theme.R                  #   ggplot2 theme for publication figures
│   ├── geo/                     #   GEO data preprocessing & filtering
│   ├── 10xgenomics.com/         #   10X Genomics data processing
│   ├── mousebrain.org/          #   mousebrain.org data processing
│   ├── coexpr/                  #   Standalone coexpression matrix scripts
│   ├── function/                #   Functional coherence analysis (AUROC)
│   ├── networks/                #   Network preprocessing, rewiring & overlap
│   ├── clustering/              #   Cell clustering (hierarchical & SNN/Louvain)
│   ├── reproducibility/         #   Cross-dataset reproducibility analysis
│   ├── disease/                 #   Disease gene prediction (Phenopedia)
│   ├── benchmark/               #   Performance benchmarking
│   └── tables/                  #   Supplementary table generation
│
├── data/
│   ├── geo/
│   │   ├── raw/                 #   Raw GEO expression files
│   │   ├── processed/           #   Standardized format
│   │   ├── filtered/            #   Gene-filtered matrices
│   │   └── one-dataset-per-publication.txt  # Consistent random sample
│   ├── one-per-publication/     #   Main pipeline input (filtered expression data)
│   ├── 10xgenomics.com/         #   10X Genomics processed/filtered data
│   ├── loom/                    #   mousebrain.org loom files (filtered)
│   ├── networks/
│   │   ├── HIPPIE/              #   Protein-protein interactions
│   │   ├── OmniPath/            #   Signalling networks
│   │   ├── Reactome/            #   Metabolic pathway co-membership
│   │   ├── STRING/              #   Text-mining gene co-occurrence
│   │   ├── idx                  #   Network index file
│   │   └── rewired/             #   Generated null models (gitignored)
│   ├── go/                      #   Gene Ontology annotations
│   └── orthologs/               #   Human-mouse ortholog mappings
│
├── results/                     # Pipeline output (gitignored)
├── fig/                         # Publication figures
└── jupyter/                     # Jupyter notebooks
```

---

## Pipeline Overview

The main analysis runs as a **Nextflow DSL2** pipeline inside a Docker/Singularity container (`suzannejin/dismay:v3.5`). It has three stages:

```
                     ┌───────────────────────────┐
                     │   Expression matrices     │
                     │ data/one-per-publication/ │
                     └────────────┬──────────────┘
                                  │
                                  ▼
                  ┌──────────────────────────────────┐
                  │  Step 1: COEXPRESSION MATRICES   │
                  │  GET_COEXPRESSION_MATRIX         │
                  │  → dismay: 17 association metrics│
                  │  FILTER_COEXPRESSION_MATRIX      │
                  │  → filter by gene prevalence     │
                  └─────────────┬────────────────────┘
                                │
                    ┌───────────┴───────────┐
                    ▼                       ▼
    ┌──────────────────────────┐  ┌─────────────────────────────┐
    │ Step 2: FUNCTIONAL       │  │ Step 3: NETWORK OVERLAP     │
    │ COHERENCE (EGAD)         │  │                             │
    │ → AUROC for GO terms     │  │ REWIRE_NETWORK              │
    │ → AUROC for Reactome     │  │ → 1000 rewired null models  │
    │                          │  │ CALCULATE_OVERLAP           │
    │                          │  │ → compare to HIPPIE,        │
    │                          │  │   OmniPath, Reactome,STRING │
    └──────────────────────────┘  └─────────────────────────────┘
```

---

## Running the Pipeline

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (DSL2 compatible)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

### Execution

```bash
# Run with Docker
nextflow run main.nf -profile docker

# Run with Singularity (HPC environments)
nextflow run main.nf -profile singularity

# Run on CRG cluster
nextflow run main.nf -profile singularity,crg

# Resume after interruption
nextflow run main.nf -profile docker -resume
```

Before this please run test, for example:
```bash
nextflow run main.nf -profile singularity,crg,test
```

### Configuration

Key parameters in `nextflow.config`:

| Parameter     | Default                                    | Description                              |
|---------------|--------------------------------------------|------------------------------------------|
| `input`       | `data/one-per-publication/*.txt.gz`        | Filtered expression matrices             |
| `methods`     | 5 SJIN variants                            | Association metrics to evaluate           |
| `network`     | `data/networks/*/{human,mouse}.txt.gz`     | Reference biological networks            |
| `eval_filt`   | `['main']`                                 | Filtering thresholds (expandable to `[50,60,70,80,90,95,'main']`) |
| `max_memory`  | 84 GB                                      | Resource ceiling                         |
| `max_cpus`    | 16                                         | Resource ceiling                         |

---

## Datasets

| Source            | Count | Location                     |
|-------------------|-------|------------------------------|
| GEO               | 162   | `data/geo/`                  |
| 10X Genomics      | 10    | `data/10xgenomics.com/`      |
| mousebrain.org    | 39    | `data/loom/`                 |

For the main analyses, one dataset per publication is randomly sampled (list in `data/geo/one-dataset-per-publication.txt`) to avoid over-representing studies with many datasets.

---

## Analyses

### Coexpression network generation

Coexpression matrices are computed using the [`dismay`](https://github.com/skinnider/dismay) R package, which implements 17 association metrics. Genes absent from ≥80% of cells are excluded for the main analysis.

### Functional coherence (EGAD)

Evaluates whether coexpression networks recover known functional relationships using AUROC scores for Gene Ontology and Reactome annotations via the [`EGAD`](https://bioconductor.org/packages/EGAD/) package.

### Network overlap

Compares coexpression networks against four types of biological reference networks. Each reference is rewired 1,000 times to establish a null distribution.

| Network   | Type                            | Source |
|-----------|---------------------------------|--------|
| HIPPIE    | Protein-protein interactions    | [HIPPIE](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/) |
| OmniPath  | Signalling networks             | [OmniPath](http://omnipathdb.org/) |
| Reactome  | Metabolic pathway co-membership | [Reactome](https://reactome.org/) |
| STRING    | Text-mining co-occurrence       | [STRING](https://string-db.org/) |


### Plotting & analysis

You can use the other scripts in `bin` to analyze/plot the results of the pipeline.

## Dependencies

All dependencies are pre-installed in the container image `suzannejin/dismay:v3.5`.
