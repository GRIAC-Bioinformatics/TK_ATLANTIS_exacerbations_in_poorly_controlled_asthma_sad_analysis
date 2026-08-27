# SAD as a Predictor of Asthma Exacerbations — ATLANTIS Study

## Background

Small airway dysfunction (SAD) — measured via impulse oscillometry (IOS) — has been linked to poorer asthma outcomes. This project investigates whether SAD independently predicts **severe exacerbations** in adults with **poorly-controlled asthma** compared to those with **well-controlled asthma**, using data from the [ATLANTIS study](https://pubmed.ncbi.nlm.nih.gov/31229323/).

SAD is assessed using two definitions across three IOS parameters:

| Parameter        | Definition                                    | Standard threshold      | Severe threshold |
| ---------------- | --------------------------------------------- | ----------------------- | ---------------- |
| **R5-R20** | Frequency dependence of resistance (5–20 Hz) | Above ULN (1.645 x RSD) | Above 3 x RSD    |
| **AX**     | Area under the reactance curve                | Above ULN (1.645 x RSD) | Above 3 x RSD    |
| **X5**     | Peripheral airway reactance at 5 Hz           | Below LLN (1.645 x RSD) | Below 3 x RSD    |

The **severe SAD** thresholds are derived from sex- and height-adjusted reference equations from the ATLANTIS supplementary material, with abnormality defined at 3 residual standard deviations (RSD) rather than the standard 1.645 x RSD.

Asthma control is classified using ACQ-6 score:

- **Well-controlled**: ACQ-6 < 0.75
- **Partially controlled**: ACQ-6 0.75-1.5
- **Poorly controlled**: ACQ-6 > 1.5

---

## Repository Structure

```
.
├── R/
│   └── baseline_table_source.R              # Shared data preparation function (sourced by baseline scripts)
├── standard_SAD/
│   ├── TK_Baseline_table_asthma_control.R   # Baseline table stratified by asthma control group + SAD
│   ├── TK_SAD_univariate_model.R            # Univariate Cox + KM plots (SAD x asthma control)
│   ├── TK_SAD_multivariate_model.R          # Multivariate Cox (SAD x asthma control subgroups)
│   ├── TK_Multi_model_R520_AX_X5_whole.R   # Multivariate Cox in whole cohort (correcting for control)
│   ├── TK_Uni_mod_clinical_var.R            # Univariate Cox for individual clinical variables
│   └── Figure_1.R                           # Combines KM plots into Figure 1
├── severe_SAD/
│   ├── TK_severe_SAD_definition.R           # Computes severe SAD flags (3 x RSD) and saves lookup table
│   ├── TK_severe_SAD_baseline.R             # Baseline table stratified by severe SAD + asthma control
│   ├── Severe_SAD_univariate_model.R        # Univariate Cox + KM plots (severe SAD x asthma control)
│   └── Severe_SAD_multivariate_model.R      # Multivariate Cox (severe SAD x asthma control subgroups)
├── config.yaml                              # Local paths — NOT committed (see .gitignore)
├── config_template.yaml                     # Template for collaborators to fill in
├── .gitignore
└── README.md
```

> **Note:** Raw data (`atlantis_patient_data.csv`, `atlantis_exacerbations.Rdata`, `severe_SAD_table.csv`) are not included in this repository as they contain patient-level data subject to data sharing restrictions.

---

## Configuration

All file paths are managed through `config.yaml` in the project root. This file is listed in `.gitignore` and will never be committed.

To get started:

```bash
cp config_template.yaml config.yaml
# then open config.yaml and fill in your local paths
```

Paths are loaded in each script using:

```r
library(here)
library(yaml)
cfg <- yaml::read_yaml(here("config.yaml"))
```

---

## Scripts

### Shared

#### `R/baseline_table_source.R`

A **source script** (not run directly) that defines the `column_rename()` function. Centralises all shared data wrangling applied before creating baseline tables: ICS treatment classification, bronchial hyperresponsiveness (BHR) severity from PC20/PD20 cutoffs, ACQ-6 group assignment, standard SAD flag variables, and column renaming to publication-ready labels.

---

### Standard SAD analysis

#### `standard_SAD/TK_Baseline_table_asthma_control.R`

Creates **Table 1** — baseline characteristics stratified by asthma control group (well / partially / poorly controlled). Also produces separate tables cross-stratified by each standard SAD parameter (R5-R20, AX, X5) within each control group.

**Output:** `output_tables/asthma_control_baseline*.csv`

**Key packages:** `dplyr`, `tidyr`, `tableone`

#### `standard_SAD/TK_SAD_univariate_model.R`

Runs **univariate Cox proportional hazards regression** with each standard SAD parameter as the sole predictor of time to first exacerbation (follow-up capped at 420 days). For each SAD parameter x asthma control group combination, produces:

- A Kaplan-Meier cumulative incidence curve with risk table and inset zoom
- Hazard ratio with 95% CI and log-rank p-value annotated on the plot
- Saved `.rds` objects for well/poorly controlled + R5-R20 (used in Figure 1)

**Output:** `output_figures/sad_uni_model/uni_surv_*.png`

**Key packages:** `survival`, `survminer`, `dplyr`, `broom`, `ggplotify`

#### `standard_SAD/TK_SAD_multivariate_model.R`

Runs **multivariate Cox regression** per asthma control subgroup, adjusting each SAD parameter for: age, sex, smoking (ex / current), GINA step 4-5, prior exacerbations (>=1 last year), blood eosinophils, FEV1 % predicted, and RV/TLC.

**Output:** `output_tables/sad_multi_model/*_multivar_cox.csv`

**Key packages:** `survival`, `dplyr`, `broom`, `tidyr`

#### `standard_SAD/TK_Multi_model_R520_AX_X5_whole.R`

Runs **multivariate Cox regression in the whole cohort** (not split by asthma control group), adding poor asthma control (ACQ-6 > 1.5) as an explicit covariate alongside each standard SAD parameter. Allows assessment of the SAD-exacerbation association independent of disease control level.

**Output:** `output_tables/wholegroup_*_full_multivar_cox.csv`

**Key packages:** `survival`, `dplyr`, `broom`

#### `standard_SAD/TK_Uni_mod_clinical_var.R`

Runs **univariate Cox regression for each clinical variable individually** (sex, GINA step, smoking, FEV1, eosinophils, FeNO, prior exacerbations, all three SAD parameters) across the full cohort and each asthma control subgroup. Useful for variable screening before building the multivariate model.

**Output:** `output_tables/*_asthma_uni_clinical.csv`

**Key packages:** `survival`, `dplyr`, `broom`, `purrr`

#### `standard_SAD/Figure_1.R`

Combines the saved `.rds` KM plot objects (R5-R20, well- and poorly-controlled) into a two-panel **Figure 1** using `patchwork`, with shared legend and panel titles.

**Output:** `output_figures/sad_uni_model/uni_surv_figure1.png`

**Key packages:** `ggpubr`, `patchwork`

---

### Severe SAD analysis

The severe SAD pipeline mirrors the standard SAD analysis but uses stricter thresholds (3 x RSD). The severe SAD flags must be computed first by `TK_severe_SAD_definition.R`, which saves them to `severe_SAD_table.csv` — this file is then joined into all downstream severe SAD scripts.

#### `severe_SAD/TK_severe_SAD_definition.R`

Derives the **severe SAD classification** for each patient using sex- and height-adjusted reference equations from the ATLANTIS supplementary material. Abnormality is defined at 3 x RSD:

- `R520_higher_3rsd`: R5-R20 >= predicted + 3 x RSD
- `X5_lower_3rsd`: X5 <= predicted - 3 x RSD
- `AX_higher_3rsd`: AX >= predicted + 3 x RSD

Also includes a consistency check verifying the standard ULN computation against the existing `B_R520ABNR` flag in the source data.

**Output:** `severe_sad_table` path defined in `config.yaml`

**Key packages:** `dplyr`

#### `severe_SAD/TK_severe_SAD_baseline.R`

Creates baseline characteristic tables cross-stratified by **severe SAD status x asthma control group**. Sources `R/baseline_table_source.R` for column renaming and joins severe SAD flags from `severe_SAD_table.csv`.

**Output:** `output_tables/asthma_control_baseline_severe_*.csv`

**Key packages:** `dplyr`, `tidyr`, `tableone`

#### `severe_SAD/Severe_SAD_univariate_model.R`

Runs **univariate Cox regression and KM plots** for each severe SAD parameter x asthma control group combination. Structure mirrors `TK_SAD_univariate_model.R`.

**Output:** `output_severe_sad/uni_model_*_*.png`

**Key packages:** `survival`, `survminer`, `dplyr`

#### `severe_SAD/Severe_SAD_multivariate_model.R`

Runs **multivariate Cox regression** for each severe SAD parameter x asthma control group combination, adjusting for age, sex, smoking, GINA step 4-5, prior exacerbations, blood eosinophils, and FEV1 % predicted.

> Note: RV/TLC is not included as a covariate here, unlike in `TK_SAD_multivariate_model.R`.

**Output:** `output_severe_sad/*_multivar_cox.csv`

**Key packages:** `survival`, `dplyr`, `broom`

---

## Suggested Run Order

```
# --- Shared ---
1. R/baseline_table_source.R                         <- source only, defines column_rename()

# --- Standard SAD ---
2. standard_SAD/TK_Baseline_table_asthma_control.R
3. standard_SAD/TK_Uni_mod_clinical_var.R
4. standard_SAD/TK_SAD_univariate_model.R            <- saves .rds objects needed for Figure 1
5. standard_SAD/TK_SAD_multivariate_model.R
6. standard_SAD/TK_Multi_model_R520_AX_X5_whole.R
7. standard_SAD/Figure_1.R                           <- requires outputs from step 4

# --- Severe SAD ---
8.  severe_SAD/TK_severe_SAD_definition.R            <- must run first; creates severe_SAD_table.csv
9.  severe_SAD/TK_severe_SAD_baseline.R
10. severe_SAD/Severe_SAD_univariate_model.R
11. severe_SAD/Severe_SAD_multivariate_model.R
```

---

## Dependencies

Install all required packages with:

```r
install.packages(c(
  "survival", "survminer", "dplyr", "tidyr", "broom",
  "ggplot2", "ggpubr", "patchwork", "ggplotify",
  "tableone", "stringr", "purrr", "readr",
  "here", "yaml"
))
```

---

## Data

This analysis uses data from the **ATLANTIS** study. Data are not publicly available. Access requests should be directed to the ATLANTIS study group.

---

## Related Publications

The following ATLANTIS study publications are relevant to this analysis:

1. *Postma, Dirkje S., et al. "Exploring the relevance and extent of small airways dysfunction in asthma (ATLANTIS): baseline data from a prospective cohort study." *The Lancet Respiratory Medicine* 7.5 (2019): 402-416. DOI: [10.1016/S2213-2600(19)30049-9 ](https://doi.org/10.1016/S2213-2600(19)30049-9) - main ATLANTIS paper (study design and reference equations)*
2. *Galant, Stanley P., et al. "Assessment of the role of small airway dysfunction in relation to exacerbation risk in patients with well controlled asthma (ATLANTIS): an observational study." *The Lancet Respiratory Medicine* 13.11 (2025): 990-1000. DOI: [10.1016/S2213-2600(25)00283-8](https://doi.org/10.1016/S2213-2600(25)00283-8) - SAD and asthma control paper*

---

## Contact

**Tatiana Karp**
*t.karp@umcg.nl*

---

## Citation

*Galant, MD Stanley P., et al. "Small Airways Disease predicts exacerbation risk in well-controlled, but not in poorly-controlled asthma: a post-hoc analysis of the ATLANTIS study." *The Journal of Allergy and Clinical Immunology: In Practice (2026). DOI: [10.1016/j.jaip.2026.07.039](https://www.jaci-inpractice.org/article/S2213-2198(26)00636-7/abstract)
