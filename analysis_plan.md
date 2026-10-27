# Analysis plan

Drafted on 19 October 2026 and committed here before any methylation data
were downloaded or examined; the commit history of this file is the record.
The primary analysis below will not be changed after seeing results; anything
added later is reported as exploratory.

## Question

In adults, is epigenetic age acceleration higher in survivors of childhood
cancer than in community controls with no cancer history?

## Design

Secondary analysis of public St. Jude Lifetime Cohort (SJLIFE) DNA methylation
data on NCBI GEO (SuperSeries GSE197678). Cross-sectional: one blood sample per
person, Illumina MethylationEPIC array.

| Item | Decision |
|---|---|
| Controls | GSE197676: 282 community controls with no cancer history |
| Survivors | GSE197674 (SJLIFE1, 2,138 adult survivors): about 300 matched to the controls on age and sex, selected before any clock is computed, with a fixed random seed |
| Outcome | Epigenetic age acceleration in years: the residual from regressing clock age on chronological age across the whole sample |
| Primary clock | Levine (PhenoAge) |
| Secondary clocks | Horvath, Hannum |
| Main comparison | Survivor versus control |
| Covariate | Sex (age is already removed by the residual) |
| Primary analysis | Linear regression: acceleration (PhenoAge) ~ group + sex. Report the difference in years with its 95% confidence interval and p-value; significance level 0.05 |
| Secondary analyses | The same model for Horvath and Hannum. Within survivors only: acceleration ~ chest radiotherapy + sex, if the data record chest radiotherapy |
| Exploratory | An elastic-net age predictor trained on the controls (glmnet), compared with the Horvath clock |

## Data dictionary

Exact GEO field names are filled in after `scripts/01_geo_inventory.R` has run.

| Column | Meaning | Type | Unit or codes | GEO field |
|---|---|---|---|---|
| `sample_id` | GEO sample accession (GSM) | text | | `!Sample_geo_accession` |
| `group` | Survivor or control | category | survivor, control | |
| `age` | Age when blood was drawn | continuous | years | |
| `sex` | Sex as recorded | category | female, male | |
| `diagnosis` | Childhood cancer type (survivors only) | category | as recorded | |
| `chest_rt` | Chest radiotherapy (survivors only) | category | yes, no | |
| `age_pheno`, `age_horvath`, `age_hannum` | Clock age | continuous | years | computed |
| `eaa_pheno`, `eaa_horvath`, `eaa_hannum` | Age acceleration | continuous | years | computed |

Missing values are `NA` in R.

## Data handling

- The data are de-identified. No attempt will be made to link a sample to a person.
- Downloads are kept unchanged in `data/raw/`, with a second copy elsewhere.
  Scripts read from `data/raw/` and write only to `data/processed/`.
- `data/` is never committed. Only code, summaries and figures go to GitHub.

## Reference

Qin N, et al. Genome-wide association studies identify novel genetic loci for
epigenetic age acceleration among survivors of childhood cancer. *Genome
Medicine*. 2022. https://doi.org/10.1186/s13073-022-01038-6
