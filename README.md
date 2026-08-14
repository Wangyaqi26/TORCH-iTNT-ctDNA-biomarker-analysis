# Input data schemas

The repository intentionally excludes patient-level study data. Create the following local files in `data/` from the final locked analysis dataset. Do not commit them unless public sharing is permitted by the ethics approval, patient consent, confidentiality agreements, and Data Access Committee requirements.

Binary variables must be coded as `1` (present/positive/event) and `0` (absent/negative/censored). Missing values should be blank or `NA`.

## `clinical_data.csv`

One row per patient.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `arm` | `Consolidation` or `Induction` |
| `age` | Age at enrollment, years |
| `sex` | `Male` or `Female` |
| `stage` | Clinical stage, e.g. `II` or `III` |
| `cMRF` | Baseline mesorectal fascia involvement, `1/0` |
| `cEMVI` | Baseline extramural vascular invasion, `1/0` |
| `baseline_CEA` | Baseline carcinoembryonic antigen |
| `post_CEA` | Post-iTNT carcinoembryonic antigen |
| `mrTRG_T3` | Mid-iTNT mrTRG score, ordinal 1-5 |
| `mrTRG_T4` | Post-iTNT mrTRG score, ordinal 1-5 |
| `pTRG` | Postoperative pathological tumor regression grade; blank for non-surgical patients |
| `treatment_choice` | `W&W`, `Surgery`, `Refused surgery`, or `PD` |
| `CR` | Complete response, `1/0` |
| `TMB` | Tumor mutational burden in mut/Mb; confirm final model scale |
| `SMAD4_mut` | Baseline SMAD4 alteration, `1/0` |
| `FLT3_amp` | Baseline FLT3 amplification, `1/0` |
| `DFS_time` | Disease-free survival time |
| `DFS_event` | DFS event, `1/0` |
| `LRFS_time` | Local recurrence-free survival time |
| `LRFS_event` | LRFS event, `1/0` |
| `DMFS_time` | Distant metastasis-free survival time |
| `DMFS_event` | DMFS event, `1/0` |
| `OS_time` | Overall survival time |
| `OS_event` | Death, `1/0` |

Use one consistent time unit, preferably months, for all survival times. Confirm whether the origin is enrollment, treatment initiation, iTNT completion, or another prespecified date.

## `ctdna_long.csv`

One row per patient and plasma time point.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `timepoint` | `T1`, `T2`, `T3`, or `T4` |
| `ctDNA_status` | Tissue-informed ctDNA status, `1/0` |

The following variables are derived before figure generation:

- `ctDNA_T1` to `ctDNA_T4` from individual time points;
- `ctDNA_T23`, the latest available result at T2 or T3;
- `ctDNA_T234`, the latest available result at T2, T3, or T4.
- `ctDNA_clearance_T2` to `ctDNA_clearance_T4`, defined only among patients with baseline-positive ctDNA as `1` for subsequent negative status and `0` for persistent positivity.

The latest value is selected by time-point order and is not a pooled or aggregate result.

## `torch_population.csv` (optional; Table S1)

One row per patient in the overall evaluable TORCH population, including both patients included and not included in the biomarker sub-cohort.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `included_biomarker_cohort` | Included in the 63-patient biomarker cohort, `1/0` |
| `arm` | Treatment arm |
| `age` | Age at enrollment |
| `sex` | Sex |
| `stage` | Clinical stage |
| `cMRF` | Baseline mesorectal fascia involvement, `1/0` |
| `cEMVI` | Baseline extramural vascular invasion, `1/0` |

## `mutations_long.csv`

One row per patient-gene-alteration record from the final processed baseline tissue variant table.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `gene` | HGNC gene symbol |
| `alteration` | Standardized alteration class |
| `arm` | Treatment arm |
| `CR` | Complete response, `1/0` |

Accepted alteration labels in the plotting script include `Missense`, `Inframe_indel`, `CNV`, `Amplification`, `Deletion`, `Frameshift`, `Splicing`, `Stop_gained`, `Stop_lost`, and `Fusion`.

## `gene_order.csv`

Optional single-column file named `gene`, containing genes in the desired Oncoprint row order. When absent, the Oncoprint script selects the 30 most frequently altered genes.

## `swimmer_data.csv`

One row per patient.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `arm` | Treatment arm |
| `CR` | Complete response, `1/0` |
| `followup_months` | Total follow-up duration |
| `ongoing_followup` | Arrow at last follow-up, `1/0`; optional |
| `surgery_months` | Time to surgery; blank if no surgery |
| `local_recurrence_months` | Time to pelvic/local recurrence |
| `distant_metastasis_months` | Time to distant metastasis |
| `death_months` | Time to death |

The four wide event columns above are the simplest input option. For a
publication-style multilayer swimmer plot, the following two optional long-form
tables can be supplied instead.

## `swimmer_intervals.csv` (optional)

One row per patient interval.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `start_month` | Beginning of the interval |
| `end_month` | End of the interval |
| `interval_type` | Label such as `iTNT`, `W&W`, `Post-surgery follow-up`, or `Progression` |

## `swimmer_events.csv` (optional)

One row per patient event. When provided, this table takes precedence over the
wide event columns in `swimmer_data.csv`.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `event_month` | Time of event |
| `event_type` | `Surgery`, `Local recurrence`, `Distant metastasis`, `Death`, or another prespecified label |

## `immune_cytokine_long.csv`

One row per patient, marker, and time point.

| Column | Definition |
| --- | --- |
| `patient_id` | De-identified patient identifier |
| `assay` | `mIHC` or `Cytokine` |
| `marker` | Marker name, e.g. `CD8`, `CD8_PD1`, `CD4`, `CD4_FOXP3` |
| `timepoint` | `T1`/`Baseline` or `T4`/`Post` |
| `value` | Quantitative measurement |
| `CR` | Complete response, `1/0` |

Confirm the exact assay-specific normalization, units, and handling of values below the detection limit before public release.
