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

Later scripts (clocks on the controls, matching, the main comparison) are
added as the project moves on.

## Status

- [x] Analysis plan fixed
- [ ] GEO inventory and sample characteristics
- [ ] Clocks on the 282 controls
- [ ] Survivors versus controls
- [ ] Write-up

## Author

Dr. Mirembe Joselyn Nayiga, MBChB
