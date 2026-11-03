# sjlife-epigenetic-age

Epigenetic age acceleration in adult survivors of childhood cancer: a
pre-registered reanalysis of public data.

**Question.** Do adult survivors of childhood cancer have a higher epigenetic
age relative to their chronological age than people who never had cancer, and
is this associated with chest radiotherapy?

**Prior work.** Qin et al. (*JNCI* 2021) found higher PhenoAge acceleration in
St. Jude Lifetime Cohort (SJLIFE) survivors than in community controls,
associated with specific treatments. Dong et al. (*Genome Medicine* 2022)
studied its genetic basis and deposited the methylation data used here. This
project is a **reproduction and extension** of that work. It adds a
pre-registered plan, individual 1:1 matching on age and sex, three clocks
computed in R with methylclock instead of the online calculator, and the chest
radiotherapy comparison for all three clocks. It uses the same controls, so it
is not an independent replication.

**Data.** Public, de-identified DNA methylation data (Illumina
MethylationEPIC, processed average beta values) from SJLIFE on NCBI GEO,
SuperSeries
[GSE197678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197678).
No data are stored in this repository; the scripts download them from GEO. No
controlled-access data are used.

**Plan.** The [analysis plan](analysis_plan.md) was committed before any
methylation values were downloaded. Survivors were matched to controls on age
and sex before any survivor methylation values were opened. Every later change
is a dated amendment.

**Why matching.** There are two reasons. By design, 1:1 matching gives two
groups with the same age and sex distribution; the adult survivors available
were younger than the controls (median 32.2 vs 35.0 years). In practice,
computing clocks for 281 survivors instead of 2,068 kept the 7 GB survivor
file within free cloud computing limits. Using all adult survivors with age
adjustment would give more statistical power and is a possible extension.

## Results

Pre-specified analyses ([results/05_results.txt](results/05_results.txt)).
Age acceleration is the residual of clock age on chronological age; models
are adjusted for sex.

| Comparison | Clock | n | Difference, years (95% CI) | p |
|---|---|---|---|---|
| Survivors vs controls (primary) | PhenoAge | 563 | +3.79 (2.68 to 4.91) | < 0.001 |
| Survivors vs controls | Horvath | 563 | +1.63 (0.97 to 2.29) | < 0.001 |
| Survivors vs controls | Hannum | 563 | +2.16 (1.39 to 2.93) | < 0.001 |
| Chest radiotherapy vs none, survivors | PhenoAge | 281 | +4.95 (3.03 to 6.86) | < 0.001 |
| Chest radiotherapy vs none, survivors | Horvath | 281 | +2.51 (1.43 to 3.60) | < 0.001 |
| Chest radiotherapy vs none, survivors | Hannum | 281 | +2.72 (1.39 to 4.05) | < 0.001 |

Exploratory, post hoc checks (added after the results above were seen) are in
`results/06_sensitivity.txt`.

**Clock validation.** In both groups the clocks tracked chronological age
closely (r = 0.88 to 0.92; Dong et al. reported 0.87 to 0.95 with the online
calculator). In these mononuclear-cell samples, Hannum and PhenoAge read on
average about 10 and 20 years below chronological age; the residual
definition removes this constant offset
([figure](results/03_controls_clock_vs_age.png),
[summary](results/03_controls_clock_summary.txt)). Clocks were computed from
the depositors' processed beta values without further normalisation.

## Limitations

- Cross-sectional: one sample per person shows a difference in epigenetic age,
  not a measured rate of ageing.
- Controls and survivors are separate GEO series and may have been processed
  in different batches, so group is confounded with series.
- No adjustment for blood cell composition in the pre-specified analysis.
- The pre-specified chest radiotherapy model is adjusted for sex only;
  treatments overlap, so it is an association, not a causal estimate.
- Treatment is recorded as yes/no; cancer diagnosis is not available in GEO.
- One US cohort; results may not apply to survivors elsewhere, including
  Africa.

## How to run

Everything runs in R (4.3 or later), on a laptop or in Google Colab (R
runtime). Run the scripts from the repository folder, in order. Colab deletes
its files when a session ends, so download each script's outputs before
moving on.

| Script | What it does | Data used |
|---|---|---|
| `scripts/00_setup.R` | Installs the packages; saves their versions to `results/session_info.txt` | small |
| `scripts/01_geo_inventory.R` | Lists the GEO files and their sizes; reads sample characteristics without downloading methylation values | a few MB |
| `scripts/02_match_survivors.R` | Chooses one adult survivor per control, matched on age and sex, before any methylation data are opened | none (uses script 01 output) |
| `scripts/03_controls_clocks.R` | Horvath, Hannum and PhenoAge clocks for the 282 controls, checked against chronological age. Run in Colab with `notebooks/03_controls_clocks.ipynb` | about 1 GB |
| `scripts/04_survivors_clocks.R` | The same clocks for the 281 matched survivors. Run in Colab with `notebooks/04_survivors_clocks.ipynb` | about 7 GB |
| `scripts/05_compare.R` | The pre-specified comparison | none (uses outputs of 03 and 04) |
| `scripts/06_sensitivity.R` | Exploratory, post hoc checks | none (uses outputs of 03 and 04) |

## Status

- [x] Analysis plan fixed
- [x] GEO inventory and sample characteristics
- [x] Survivors matched to controls (281 pairs)
- [x] Clocks on the 282 controls
- [x] Clocks on the 281 matched survivors
- [x] Survivors versus controls
- [x] Write-up
- [ ] Exploratory: cell-composition-adjusted acceleration, batch and sex checks

## References

1. Qin N, Li Z, Song N, et al. Epigenetic age acceleration and chronic health
   conditions among adult survivors of childhood cancer. *J Natl Cancer Inst*.
   2021;113(5):597-605. https://doi.org/10.1093/jnci/djaa147
2. Dong Q, Song N, Qin N, et al. Genome-wide association studies identify novel
   genetic loci for epigenetic age acceleration among survivors of childhood
   cancer. *Genome Med*. 2022;14:32.
   https://doi.org/10.1186/s13073-022-01038-6

## How to cite

Nayiga Joselyn Mirembe (2026). *sjlife-epigenetic-age: epigenetic age
acceleration in adult survivors of childhood cancer, a pre-registered
reanalysis of public data*. GitHub repository,
https://github.com/Mrn-selyn-commits/sjlife-epigenetic-age

## Licence

Code: MIT licence (see [LICENSE](LICENSE)). The data belong to their
depositors and are used under NCBI GEO's terms.

AI tools (Claude, Anthropic) were used responsibly in developing this project.

## Author

Dr. Nayiga Joselyn Mirembe, MBChB
