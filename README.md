# TORCH biomarker analysis

R scripts for the downstream statistical analyses and figures accompanying the TORCH biomarker study:

> Application of genetic features and serial ctDNA in estimating response and prognosis of immunotherapy-based total neoadjuvant therapy for microsatellite-stable locally advanced rectal cancer

## Repository structure

```text
TORCH-biomarker-analysis/
├── Fig2/
├── Fig3/
├── Fig4/
├── Fig5/
├── Supplementary/
├── R/
├── scripts/
├── data/
└── results/
```

Figure folders contain the corresponding entry scripts. Shared functions and analysis scripts are stored in `R/` and `scripts/`.

## Usage

Install the required packages:

```r
source("scripts/00_install_packages.R")
```

Prepare the local input tables following `data/README.md`, then run the required figure script from the repository root:

```r
source("Fig2/Fig2.R")
source("Fig3/Fig3.R")
source("Fig4/Fig4.R")
source("Fig5/Fig5.R")
source("Supplementary/Supplementary.R")
```

Outputs are written to `results/`. The scripts should be checked against the final analysis dataset before release.

## Data availability

Patient-level data are not included in this repository. Genomic data are deposited in the Genome Sequence Archive for Human under controlled access (HRA013886; BioProject PRJCA048403).

## Software

The manuscript reports R 4.4.1. 

