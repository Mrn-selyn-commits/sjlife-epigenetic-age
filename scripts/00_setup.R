# 00_setup.R
# Installs the packages the project needs. Run once per machine
# (in Google Colab or Kaggle, run it at the start of every session).
#
# On Linux cloud machines (Colab, Kaggle) it first installs the system
# libraries that the Bioconductor packages need, and uses pre-built CRAN
# packages so installation takes minutes rather than hours.

is_linux_root <- Sys.info()[["sysname"]] == "Linux" &&
  Sys.info()[["user"]] == "root" && nzchar(Sys.which("apt-get"))

if (is_linux_root) {
  message("Installing system libraries...")
  system(paste("apt-get update -qq && apt-get install -y -qq",
               "libxml2-dev libssl-dev libcurl4-openssl-dev libgsl-dev",
               "libglpk-dev libhdf5-dev libpng-dev libjpeg-dev"))
  os <- readLines("/etc/os-release")
  codename <- sub("^VERSION_CODENAME=", "", grep("^VERSION_CODENAME=", os, value = TRUE))
  if (length(codename) == 1 && nzchar(codename)) {
    options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)),
            HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
              paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
    message("Using pre-built CRAN packages for Ubuntu ", codename)
  }
}

cran_pkgs <- c("ggplot2", "data.table", "R.utils", "BiocManager")
missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

bioc_pkgs <- c("methylclock")
missing <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)

ok <- vapply(c(cran_pkgs, bioc_pkgs), requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Not installed: ", paste(names(ok)[!ok], collapse = ", "),
                   ". Scroll up to the first line that says ERROR and send it.")

# Record the versions used, so the analysis can be rerun exactly
dir.create("results", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("Setup complete. Versions saved to results/session_info.txt")
