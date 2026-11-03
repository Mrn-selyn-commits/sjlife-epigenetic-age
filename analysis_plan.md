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

Each amendment is dated. Amendments 1 to 4 were made before any methylation
values were downloaded. Amendments 5 to 8 were made after the controls' clocks
had been computed (script 03), but before any survivor methylation values were
opened. Amendment 9 was made after the primary results had
been seen. The original text above is otherwise unchanged.

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

**3 November 2026, after computing clocks for the controls only (script 03), before any survivor methylation data were opened:**

5. *Clock offsets in controls.* In these peripheral blood mononuclear cell samples,
   all three clocks track chronological age closely (r = 0.92 Horvath, 0.90
   Hannum, 0.88 PhenoAge), but Hannum and PhenoAge read on average about 10 and
   20 years below chronological age (Horvath +1.4 years). Both clocks were
   developed in whole blood. The offsets are removed by the pre-specified
   definition of age acceleration (residual of clock age on chronological age),
   so the analysis is unchanged; the point is reported as a limitation.
6. *GEO labelling.* GEO describes the control samples as "blood from survivor of
   childhood cancer", while their IDs (SJNORM…) and the series identify them as
   community controls. This is treated as a labelling error in GEO.
7. *Implementation only.* Clocks are computed in batches of 30 samples, and the
   survivors' columns are extracted from the GEO file before computing, to stay
   within cloud memory limits. Each clock is a per-sample weighted sum, so this
   does not change any estimate.
8. *Secondary analysis detail.* The chest radiotherapy comparison is reported for
   all three clocks, PhenoAge first.

**3 November 2026, after the primary results (script 05) had been seen:**

9. *Corrections and post hoc analyses.* No pre-specified analysis is changed.
   - The introduction to this section previously said every amendment was made
     before any methylation values were downloaded or examined. That was true
     only of amendments 1 to 4; the wording above now states when each was made.
   - Amendment 7 called each clock "a per-sample weighted sum". More exactly,
     each clock is computed per sample from a weighted sum of CpGs plus an
     intercept (Horvath's is then transformed back to years). In methylclock,
     clock CpGs absent from the data are left out of the sum rather than
     imputed, which is a likely contributor to the offsets in amendment 5, and
     missing values at measured CpGs are imputed from other samples in the same
     batch. Whether any such values occurred has not yet been checked.
   - The reference below is corrected: the Genome Medicine 2022 paper is by
     Dong, Song, Qin et al.; Qin et al. (JNCI 2021) is added. This project is a
     reproduction and extension of that work, with the same controls.
   - Following external review, the analyses in `scripts/06_sensitivity.R` were
     added. They are exploratory and post hoc: the group difference estimated
     without residuals and with the age line fitted in controls only; a
     group-by-age interaction; effect sizes in SD units; Holm correction of the
     secondary tests; the smallest difference the design could detect; and the
     chest radiotherapy model adjusted for age and all other treatment flags.
   - Known limitation added: controls (GSE197676) and survivors (GSE197674) are
     separate GEO series and may have been processed in different batches, so
     group is confounded with series. This design cannot remove that.

## Data handling

- The data are de-identified. No attempt will be made to link a sample to a person.
- Downloads are kept unchanged in `data/raw/`, with a second copy elsewhere.
  Scripts read from `data/raw/` and write only to `data/processed/`.
- `data/` is never committed. Only code, summaries and figures go to GitHub.

## References

1. Qin N, Li Z, Song N, et al. Epigenetic age acceleration and chronic health
   conditions among adult survivors of childhood cancer. *Journal of the
   National Cancer Institute*. 2021;113(5):597-605.
   https://doi.org/10.1093/jnci/djaa147
2. Dong Q, Song N, Qin N, et al. Genome-wide association studies identify novel
   genetic loci for epigenetic age acceleration among survivors of childhood
   cancer. *Genome Medicine*. 2022;14:32.
   https://doi.org/10.1186/s13073-022-01038-6
