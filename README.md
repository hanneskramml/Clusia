![GitHub Release](https://img.shields.io/github/v/release/hanneskramml/clusia)
[![DOI](https://zenodo.org/badge/804242299.svg)](https://doi.org/10.5281/zenodo.19212189)

# _Clusia_ code repository
Supplementary Code accompanying the manuscript **_Clusia_ genomes shed light on the evolution and diversity of CAM physiotypes (_Nature Communications_, 2026)**

## 📘 Overview
Pipelines and scripts to build and explore the _Clusia_ genome assemblies, annotations, and multiomics analyses. The `clusia_panomics_database` handles data processing and integration (from raw data to figures and tables).

## 📁 Repository Structure
```
Clusia/
 ├── assembly/                      # Genome assembly, scaffolding, and ploidy/subgenome identification scripts
 ├── annotation/                    # Pipelines for gene, repeat, and pseudogene annotation
 ├── comparative_genomics/          # Orthogroup sampling, alignments, features incl. cis-regulatory motifs
 ├── rnaseq/                        # RNA-seq/gene expression workflow
 ├── greenhouse_experiment/         # Environmental monitoring app written in Python/C (RaspberryPi/Arduino)
 ├── clusia_panomics_database/      # R pipeline for data processing and integration (figure and table generation)
 └── README.md                      # Project documentation
```

## 🚀 Getting Started
This section provides instructions on how to run the `clusia_panomics_database`, which generates figures and tables based on the raw data processed by the scripts and pipelines in the `assembly`, `annotation`, `comparative_genomics`, and `rnaseq` folders.

### Requirements
- A Windows/Linux/MacOS system with R installed (tested on R4.2.1)
- CRAN packages: BiocManager (1.30.25), circlize (0.4.16), curl (5.2.3), tidyverse (2.0.0), viridis (0.6.5)
- BioConductor packages: BiocVersion (3.16.0), Rsamtools (2.14.0), clusterProfiler (4.6.2), cogeqc (1.2.1), enrichplot (1.18.4), org.At.tair.db (3.16.0)
- GitHub packages: GENESPACE (1.2.3)

### Clone the repository
```bash
git clone https://github.com/hanneskramml/Clusia.git
cd Clusia
git checkout Manuscript1
```

### Get the rawdata
- Link to **figshare**: https://doi.org/10.6084/m9.figshare.27599406
- Download the following raw data and place it in the `data` folder:
- Dataset 2-6, Dataset 9-11

### Run analyses
Set paths in `ClusiaDB.R` accordingly. The entire pipeline can be executed by running one of the following scripts:
```bash
cd clusia_panomics_database
./run
```

Or directly in R:
```r
source("ClusiaDB.R")
```

### Demo
A dedicated compute capsule including all dependencies, Supplementary Data, and final results is published on **CodeOcean**: https://doi.org/10.24433/CO.2105665.v2
(runtime: approx. 40min)


## 📊 Reproducibility
- Scripts are modular and can be executed independently (unless stated otherwise).
- Data-processing steps are documented within each script.
- Figures and tables used in the manuscript are generated programmatically.

## 📝 Manuscript Integration
- This branch corresponds specifically to **_Clusia_ genomes shed light on the evolution and diversity of CAM physiotypes (_Nature Communications_, 2026)**.  
- Versioned analyses correspond to those used in the submitted/published manuscript.  
- Updates after peer review can be tracked via commits, tags and releases.

## 📄 License
This project is covered under the **GPL-3.0 License**.