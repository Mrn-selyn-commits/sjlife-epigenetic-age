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

GEO field names confirmed by `scripts/01_geo_inventory.R` (27 October 2026).

| Column | Meaning | Type | Unit or codes | GEO field |
|---|---|---|---|---|
| `sample_id` | GEO sample accession (GSM) | text | | `!Sample_geo_accession` |
| `group` | Survivor or control | category | survivor, control | series: GSE197674 survivors, GSE197676 controls |
| `age` | Age when blood was drawn | continuous | years | `age` |
| `sex` | Sex as recorded | category | Female, Male | `gender` |
| `diagnosis` | Childhood cancer type (survivors only) | category | as recorded | not in GEO |
| `chest_rt` | Chest radiotherapy (survivors only) | category | 1 = yes, 0 = no | `chestrt` |
| `age_pheno`, `age_horvath`, `age_hannum` | Clock age | continuous | years | computed |
| `eaa_pheno`, `eaa_horvath`, `eaa_hannum` | Age acceleration | continuous | years | computed |

Missing values are `NA` in R.

## Amendments

Each amendment is dated and was made before any methylation values were
downloaded or examined. The original text above is otherwise unchanged.

**27 October 2026, after reading sample characteristics only (script 01):**

1. *Survivor sample.* "About 300 survivors" is made exact: 1:1 nearest-age
   matching within sex, without replacement, one survivor per control, caliper
   2 years, seed 20261027, survivors aged 18 or over at sampling
   (`scripts/02_match_survivors.R`). Controls without a survivor within 2 years
   stay in the analysis.
2. *Diagnosis.* GEO does not record the cancer diagnosis, so it is dropped from
   the data dictionary. It was never part of the analysis.
3. *Sample type.* GEO describes the samples as peripheral blood mononuclear
   cells, while the publication describes blood-derived DNA. This will be
   reported as a limitation; the analysis is unchanged.
4. *Methylation files.* The series matrices hold sample information only. The
   methylation values are in each series' supplementary `processed.txt.gz`
   file (controls about 1 GB, SJLIFE1 about 7 GB, compressed).

## Data handling

- The data are de-identified. No attempt will be made to link a sample to a person.
- Downloads are kept unchanged in `data/raw/`, with a second copy elsewhere.
  Scripts read from `data/raw/` and write only to `data/processed/`.
- `data/` is never committed. Only code, summaries and figures go to GitHub.

## Reference

Qin N, et al. Genome-wide association studies identify novel genetic loci for
epigenetic age acceleration among survivors of childhood cancer. *Genome
Medicine*. 2022. https://doi.org/10.1186/s13073-022-01038-6
