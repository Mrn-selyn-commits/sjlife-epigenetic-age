# 00_setup.R
# Installs the packages the project needs. Run once per machine
# (in Google Colab or Kaggle, run it at the start of every session).
#
# Script 01 needs only base R. The packages below are for the later steps.

cran_pkgs <- c("ggplot2", "data.table", "BiocManager")
missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

bioc_pkgs <- c("GEOquery", "methylclock")
missing <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)

# Record the versions used, so the analysis can be rerun exactly
dir.create("results", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("Setup complete. Versions saved to results/session_info.txt")
