# 04_survivors_clocks.R
# Phase 3 of the analysis plan, part 1: compute the Horvath, Hannum and
# Levine (PhenoAge) clocks for the 281 matched survivors chosen by script 02
# (results/02_selected_survivors.txt), from SJLIFE1 (GSE197674).
#
# Run in Google Colab (R runtime). The GEO file is about 7 GB compressed:
# the script downloads it once, then copies out only the 281 selected
# samples' columns in a single pass, so memory and disk stay manageable.
# Takes roughly 30-60 minutes, mostly downloading and unpacking.
#
# Outputs (small, commit these):
#   results/04_survivors_file_check.txt
#   results/04_survivors_clock_summary.txt
# Output (sample-level, keep out of GitHub; needed by script 05):
#   data/processed/04_survivors_clock_ages.csv

suppressPackageStartupMessages({ library(data.table); library(methylclock) })
source("scripts/geo_helpers.R")
source("scripts/link_columns.R")
options(timeout = 14400)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

base      <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197674/"
proc_name <- "GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_processed.txt.gz"
meta_name <- "GSE197674_SJLIFE1_IlluminaEPIC_2138samples_GEO_03012022_metadata.txt.gz"
proc_file <- file.path("data/raw", proc_name)
sub_file  <- "data/raw/04_selected_survivors_values.txt"

# ---------------------------------------------- 1. who was selected -------
selected <- readLines("results/02_selected_survivors.txt")
selected <- selected[nzchar(selected)]
cat("Survivors selected by script 02:", length(selected), "\n")

# ---------------------------------------------- 2. sample information ------
matrix_file <- fetch_geo(paste0(base, "matrix/GSE197674_series_matrix.txt.gz"),
                         "data/raw/GSE197674_series_matrix.txt.gz")
fields <- read_sample_rows(matrix_file)
keys   <- vapply(fields, `[`, "", 1)
gsm    <- fields[[match("!Sample_geo_accession", keys)]][-1]
info   <- series_characteristics(fields, gsm)          # age, gender, treatment flags
names(info)[names(info) == "gender"] <- "sex"
cat("Samples in series matrix:", nrow(info), "\n")
if (!all(selected %in% info$sample_id)) stop("Some selected survivors are not in GSE197674.")

# ---------------------------------------------- 3. download ---------------
if (!file.exists(sub_file)) {
  if (!file.exists(proc_file)) cat("Downloading", proc_name, "(about 7 GB; 10-40 minutes)...\n")
  fetch_geo(paste0(base, "suppl/", proc_name), proc_file)
}
meta_file <- fetch_geo(paste0(base, "suppl/", meta_name), file.path("data/raw", meta_name))
meta <- tryCatch(fread(meta_file), error = function(e) NULL)

# ---------------------------------------------- 4. find the columns --------
header <- NULL
if (!file.exists(sub_file)) {
  con <- gzfile(proc_file, "r"); header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]; close(con)
}
if (!is.null(header)) {
  header <- gsub('^"|"$', "", header)
  pval_idx <- grep("pval|detection", header, ignore.case = TRUE)
  beta_idx <- setdiff(seq_along(header)[-1], pval_idx)
  link <- link_columns(header[beta_idx], fields, gsm, meta)
  col_gsm <- link$gsm
  want <- beta_idx[col_gsm %in% selected]
  found <- unique(col_gsm[col_gsm %in% selected])
  check <- c(sprintf("File: %s", proc_name),
             sprintf("Columns: %d (value columns: %d; detection p-value columns: %d)",
                     length(header), length(beta_idx), length(pval_idx)),
             sprintf("First value column names: %s", paste(head(header[beta_idx], 5), collapse = ", ")),
             "", "Linking columns to samples:", link$report, sprintf("Method used: %s", link$method),
             sprintf("Selected survivors found in the file: %d of %d", length(found), length(selected)))
  writeLines(check, "results/04_survivors_file_check.txt")
  cat(check, sep = "\n")
  if (link$method == "none" || length(found) < length(selected))
    stop("Not all selected survivors could be found. Send me the text printed above.")

  # ---------------------------------------------- 5. extract once ----------
  cat("Copying out the", length(want), "selected columns (one pass through the file; 10-30 minutes)...\n")
  tmp <- paste0(sub_file, ".part")
  cmd <- sprintf("gzip -dc %s | cut -f %s > %s", shQuote(proc_file),
                 paste(c(1, want), collapse = ","), shQuote(tmp))
  if (system(cmd) != 0 || !file.exists(tmp)) stop("Extracting the columns failed. Run this cell again.")
  file.rename(tmp, sub_file)   # only a complete extract gets the final name
  unlink(proc_file)   # free about 7 GB of disk; the extract is all we need
}

# ---------------------------------------------- 6. compute the clocks ------
sub_header <- names(fread(sub_file, nrows = 0))
peek <- fread(sub_file, nrows = 2000)
vals <- unlist(peek[, 2:min(6, ncol(peek))])
is_beta <- min(vals, na.rm = TRUE) >= 0 && max(vals, na.rm = TRUE) <= 1
rm(peek)
sub_link <- link_columns(sub_header[-1], fields, gsm, meta)
ids <- sub_link$gsm
if (anyNA(ids) || !all(selected %in% ids)) stop("Lost track of sample IDs after extraction.")
cat(sprintf("Computing clocks for %d survivors (%d CpGs)...\n", length(ids), nrow(fread(sub_file, select = 1))))
ages <- clocks_in_batches(sub_file, sub_header[1], sub_header[-1], ids, is_beta = is_beta)

# ---------------------------------------------- 7. save --------------------
out <- merge(info[info$sample_id %in% selected, ], ages, by = "sample_id")
out$age <- as.numeric(out$age)
fwrite(out, "data/processed/04_survivors_clock_ages.csv")

stats <- rbindlist(lapply(c(PhenoAge = "age_pheno", Horvath = "age_horvath", Hannum = "age_hannum"), function(col)
  data.table(n = sum(!is.na(out[[col]])), r = cor(out$age, out[[col]], use = "complete.obs"))), idcol = "clock")
summary_lines <- c(
  sprintf("Matched survivors with clock estimates: %d of %d", nrow(out), length(selected)),
  "",
  sprintf("%-9s n = %d | correlation with age r = %.2f", stats$clock, stats$n, stats$r),
  "",
  "Comparison with controls is done in script 05, as set out in the analysis plan."
)
writeLines(summary_lines, "results/04_survivors_clock_summary.txt")
cat(summary_lines, sep = "\n")
cat("\nDone. Download results/04_* (for GitHub) and data/processed/04_survivors_clock_ages.csv (keep on your computer).\n")
