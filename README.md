# sjlife-epigenetic-age

Epigenetic age acceleration in adult survivors of childhood cancer: a
reanalysis of public data.

**Question.** Do adult survivors of childhood cancer age faster biologically
than people who never had cancer, and is that linked to their treatment?

**Data.** St. Jude Lifetime Cohort (SJLIFE) blood DNA methylation, Illumina
MethylationEPIC, from NCBI GEO SuperSeries
[GSE197678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197678).
The data are not stored in this repository.

**Plan.** The analysis plan was fixed before the data were opened:
[analysis_plan.md](analysis_plan.md).

## How to run

Everything runs in R (4.3 or later), on a laptop or in a free cloud notebook
(Google Colab with the R runtime, or Kaggle). Run the scripts from the
repository folder, in order:

| Script | What it does | Data used |
|---|---|---|
| `scripts/00_setup.R` | Installs the packages; saves their versions | small |
| `scripts/01_geo_inventory.R` | Lists the GEO files and their sizes; reads the sample characteristics without downloading methylation values | a few MB |
| `scripts/02_match_survivors.R` | Chooses one adult survivor per control, matched on age and sex, before any methylation data are opened | none (uses script 01 output) |
| `scripts/03_controls_clocks.R` | Computes the Horvath, Hannum and PhenoAge clocks for the 282 controls and checks them against chronological age. Run it in Colab with `notebooks/03_controls_clocks.ipynb` | about 1 GB |
| `scripts/04_survivors_clocks.R` | Computes the same clocks for the 281 matched survivors. Run it in Colab with `notebooks/04_survivors_clocks.ipynb` | about 7 GB |
| `scripts/05_compare.R` | The pre-specified comparison: age acceleration in survivors vs controls, and by chest radiotherapy. Runs on any computer | none (uses outputs of 03 and 04) |

Later scripts (clocks on the controls, matching, the main comparison) are
added as the project moves on.

## Status

- [x] Analysis plan fixed
- [x] GEO inventory and sample characteristics
- [x] Survivors matched to controls (281 pairs)
- [x] Clocks on the 282 controls
- [ ] Clocks on the 281 matched survivors
- [ ] Survivors versus controls
- [ ] Write-up

## Author

Dr. Mirembe Joselyn Nayiga, MBChB
