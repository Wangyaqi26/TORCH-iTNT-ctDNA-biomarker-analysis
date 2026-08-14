# TORCH biomarker analysis

R code accompanying the manuscript:

> Application of genetic features and serial ctDNA in estimating response and prognosis of immunotherapy-based total neoadjuvant therapy (iTNT) for microsatellite-stable locally advanced rectal cancer: data from the TORCH trial

This repository contains the downstream statistical analysis and figure-generation code used for the post hoc translational biomarker analysis of the prospective TORCH trial (NCT04518280). It does **not** contain FASTQ processing, read alignment, germline calling, somatic variant calling, ATG-Seq implementation, or other proprietary upstream bioinformatics workflows.

## Important status note

`scripts/03_oncoprint.R` was rewritten from the plotting script supplied by the study team. The remaining scripts are clean, manuscript-aligned reproducibility implementations reconstructed from the Statistical analysis section, Results, figure legends, and supplementary-item descriptions. Before public release, the authors must run the scripts on the final analysis tables and verify that all numerical results and plots match the accepted manuscript. See `REVIEW_BEFORE_RELEASE.md`.

## Repository structure

```text
TORCH-biomarker-analysis/
├── R/
│   └── functions.R
├── scripts/
│   ├── 00_install_packages.R
│   ├── 01_prepare_analysis_data.R
│   ├── 02_cohort_summary.R
│   ├── 03_oncoprint.R
│   ├── 04_response_associations.R
│   ├── 05_logistic_loocv_roc.R
│   ├── 06_survival_analysis.R
│   ├── 07_swimmer_plot.R
│   ├── 08_immune_cytokine_analysis.R
│   └── 99_session_info.R
├── data/
│   ├── README.md
│   └── templates/
├── results/
├── CITATION.cff
└── REVIEW_BEFORE_RELEASE.md
```

## Analysis-to-manuscript map

| Script | Main purpose | Manuscript outputs |
| --- | --- | --- |
| `01_prepare_analysis_data.R` | Derive T1-T4 and latest available ctDNA variables | Table S2; model and survival inputs |
| `02_cohort_summary.R` | Cohort characteristics and arm comparisons | Table 1; Table S1 |
| `03_oncoprint.R` | Baseline genomic landscape | Figure 2A |
| `04_response_associations.R` | Fisher tests, BH correction, Welch tests, Cochran-Armitage trend tests, univariable logistic models | Figure 2B-D; Figure 3A-C; Tables S5-S7 |
| `05_logistic_loocv_roc.R` | Multivariable logistic regression, leave-one-out cross-validation, ROC curves and 1,000-iteration bootstrap AUC comparison | Figure 3D-I |
| `06_survival_analysis.R` | Kaplan-Meier analyses, univariable/multivariable Cox models and forest plots | Figures 4B-E and 5; Figures S4-S6; Tables S8-S11 |
| `07_swimmer_plot.R` | Patient follow-up and recurrence timeline | Figure 4A |
| `08_immune_cytokine_analysis.R` | Exploratory mIHC and cytokine comparisons | Figures S1-S3 |

## Input data

No individual-level patient data are included in this repository. Input tables should be created locally from the final locked analysis dataset using the schemas in `data/README.md`. Use de-identified study identifiers only.

The genomic dataset is deposited in the Genome Sequence Archive for Human under controlled access, accession **HRA013886**, BioProject **PRJCA048403**. Access is subject to the repository's Data Access Committee procedure and the applicable ethics approval, patient consent, and confidentiality requirements.

## Software

The manuscript reports analyses performed using R 4.4.1 and Python 3.10.14. These downstream scripts are written in R. Required packages are listed in `scripts/00_install_packages.R`.

From the repository root, install packages once:

```r
source("scripts/00_install_packages.R")
```

After placing the analysis tables in `data/`, run scripts in numerical order:

```r
source("scripts/01_prepare_analysis_data.R")
source("scripts/02_cohort_summary.R")
source("scripts/03_oncoprint.R")
source("scripts/04_response_associations.R")
source("scripts/05_logistic_loocv_roc.R")
source("scripts/06_survival_analysis.R")
source("scripts/07_swimmer_plot.R")
source("scripts/08_immune_cytokine_analysis.R")
source("scripts/99_session_info.R")
```

Outputs are written to `results/` and are ignored by Git by default.

## Reproducibility and privacy

- Run all scripts from the repository root.
- Use only de-identified patient identifiers.
- Do not commit clinical datasets, controlled-access genomic data, credentials, local absolute paths, or company-internal pipelines.
- The repository records analysis code only; access to controlled data does not imply permission to redistribute those data.
- The exact variable coding, complete-case populations, and statistical outputs must be checked against the final manuscript before a GitHub release is created.

## Repository organization reference

The one-analysis-per-script organization was informed by the public [SOG-Lab OAC IntraTumourHeterogeneity repository](https://github.com/SOG-Lab/OAC_IntraTumourHeterogeneity). No code was copied from that repository.

## Citation

After creating the public GitHub repository and Zenodo release, update `CITATION.cff` and cite both records in the manuscript. Suggested format:

> Wang Y, Xu Y, Lin Y, et al. Analysis code accompanying the TORCH biomarker study. GitHub. https://github.com/OWNER/TORCH-biomarker-analysis (2026).

> Wang Y, Xu Y, Lin Y, et al. Analysis code accompanying the TORCH biomarker study. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX (2026).

## Contact

Zhen Zhang, Fudan University Shanghai Cancer Center  
Email: zhen_zhang@fudan.edu.cn

