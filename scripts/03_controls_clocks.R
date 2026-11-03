# 03_controls_clocks.R
# Phase 2 of the analysis plan: compute the Horvath, Hannum and Levine
# (PhenoAge) clocks for the 282 community controls (GSE197676) and check that
# clock age tracks chronological age.
#
# Run in Google Colab (R runtime) or any machine with ~8 GB free RAM and
# ~5 GB free disk. It downloads one file of about 1 GB from GEO.
# Needs packages data.table, R.utils, ggplot2 and methylclock
# (installed by scripts/00_setup.R).
#
# Outputs (small, commit these):
#   results/03_controls_file_check.txt     what the methylation file contains
#   results/03_controls_clock_summary.txt  how well each clock tracks age
#   results/03_controls_clock_vs_age.png   the Phase 2 figure
# Outputs (sample-level, keep out of GitHub):
#   data/processed/03_controls_clock_ages.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(methylclock)
})
options(timeout = 7200)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

gse   <- "GSE197676"
base  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197676/"
proc_name <- "GSE197676_SJLIFE_controls_IlluminaEPIC_282samples_GEO_03012022_processed.txt.gz"
proc_url  <- paste0(base, "suppl/", proc_name)
proc_file <- file.path("data/raw", proc_name)
clocks <- c("Horvath", "Hannum", "Levine")

# ---------------------------------------------- 1. sample information ----
# Download a GEO file once into data/raw, retrying if GEO refuses briefly
# (it sometimes answers "403 Forbidden" when busy or after many requests).
fetch_geo <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  for (attempt in 1:5) {
    ok <- tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl",
                    headers = c("User-Agent" = "Mozilla/5.0 (research script; R)"))
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok && nzchar(Sys.which("curl")))
      ok <- system2("curl", c("-sSfL", "-A", shQuote("Mozilla/5.0"), "-o", shQuote(dest), shQuote(url))) == 0
    if (ok && file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
    unlink(dest)
    message(sprintf("GEO did not respond (attempt %d of 5); waiting %d seconds...", attempt, 30 * attempt))
    Sys.sleep(30 * attempt)
  }
  stop("GEO refused the download after 5 attempts: ", url,
       "\nWait 15-30 minutes and run this cell again.")
}

# Read only the sample rows of the series matrix (a few KB): GSM, title, age, sex.
read_sample_rows <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  out <- character()
  repeat {
    chunk <- readLines(con, n = 500, warn = FALSE)
    if (!length(chunk)) break
    stop_at <- grep("^!series_matrix_table_begin", chunk)
    if (length(stop_at)) { out <- c(out, chunk[seq_len(stop_at[1] - 1)]); break }
    out <- c(out, chunk)
  }
  out <- out[startsWith(out, "!Sample_")]
  lapply(out, function(r) gsub('^"|"$', "", strsplit(r, "\t", fixed = TRUE)[[1]]))
}
matrix_file <- fetch_geo(paste0(base, "matrix/GSE197676_series_matrix.txt.gz"),
                         "data/raw/GSE197676_series_matrix.txt.gz")
fields <- read_sample_rows(matrix_file)
keys   <- vapply(fields, `[`, "", 1)
gsm    <- fields[[match("!Sample_geo_accession", keys)]][-1]
char_value <- function(name) {
  for (i in which(keys == "!Sample_characteristics_ch1")) {
    v <- fields[[i]][-1]
    if (grepl(paste0("^", name, ":"), v[1], ignore.case = TRUE))
      return(trimws(sub("^[^:]*:", "", v)))
  }
  rep(NA_character_, length(gsm))
}
samples <- data.frame(sample_id = gsm,
                      age = as.numeric(char_value("age")),
                      sex = char_value("gender"))
cat("Samples in series matrix:", nrow(samples), "\n")

# ---------------------------------------------- 2. download the file -----
if (!file.exists(proc_file)) cat("Downloading", proc_name, "(about 1 GB)...\n")
fetch_geo(proc_url, proc_file)

# Unzip once (about 3-4 GB) so it is not decompressed again at every read
proc_txt <- sub("\\.gz$", "", proc_file)
if (!file.exists(proc_txt)) R.utils::gunzip(proc_file, destname = proc_txt, remove = FALSE)

# ---------------------------------------------- 3. what is in the file ---
peek   <- fread(proc_txt, nrows = 2000)
header <- names(peek)
id_col   <- header[1]
pval_col <- grep("pval|detection", header, ignore.case = TRUE, value = TRUE)
beta_col <- setdiff(header[-1], pval_col)
first_vals <- unlist(peek[, ..beta_col][, 1:min(5, length(beta_col))])
check <- c(
  sprintf("File: %s", proc_name),
  sprintf("Columns: %d (ID: '%s'; detection p-value columns: %d; value columns: %d)",
          length(header), id_col, length(pval_col), length(beta_col)),
  sprintf("First value column names: %s", paste(head(beta_col, 5), collapse = ", ")),
  sprintf("First row IDs: %s", paste(head(peek[[1]], 5), collapse = ", ")),
  sprintf("Value range in first 2000 rows (5 samples): %.3f to %.3f",
          min(first_vals, na.rm = TRUE), max(first_vals, na.rm = TRUE))
)
writeLines(check, "results/03_controls_file_check.txt")
cat(check, sep = "\n")
is_beta <- min(first_vals, na.rm = TRUE) >= 0 && max(first_vals, na.rm = TRUE) <= 1

# ---------------------------------------------- 4. link columns to GSM ---
source("scripts/link_columns.R")
meta_name <- "GSE197676_SJLIFE_controls_IlluminaEPIC_282samples_GEO_03012022_metadata.txt.gz"
meta_file <- file.path("data/raw", meta_name)
fetch_geo(paste0(base, "suppl/", meta_name), meta_file)
meta <- tryCatch(fread(meta_file), error = function(e) NULL)
link <- link_columns(beta_col, fields, gsm, meta)
col_gsm <- link$gsm

# D. The controls only need age and sex. If no exact GSM link works (or the
#    fingerprint link leaves some samples out), take them straight from the
#    metadata file, which then covers every column.
use_meta_direct <- FALSE
if ((link$method == "none" || startsWith(link$method, "C")) && !is.null(meta)) {
  m <- as.data.frame(meta)
  key <- names(m)[which.max(vapply(m, function(v) mean(norm_name(beta_col) %in% norm_name(v)), numeric(1)))]
  bare <- sub("^characteristics[ _:.-]*", "", names(m), ignore.case = TRUE)
  agec <- names(m)[grep("^age", bare, ignore.case = TRUE)[1]]
  sexc <- names(m)[grep("^(sex|gender)", bare, ignore.case = TRUE)[1]]
  if (!is.na(agec) && !is.na(sexc) && mean(norm_name(beta_col) %in% norm_name(m[[key]])) >= 0.9) {
    r <- match(norm_name(beta_col), norm_name(m[[key]]))
    samples <- data.frame(sample_id = beta_col, age = as.numeric(m[[agec]][r]), sex = m[[sexc]][r])
    col_gsm <- beta_col
    use_meta_direct <- TRUE
    link$method <- sprintf("D: age and sex taken from the metadata file ('%s', '%s'), keyed by column name", agec, sexc)
  }
}
link_report <- c("", "Linking columns to samples:", link$report,
                 sprintf("Method used: %s", link$method),
                 if (!is.null(meta)) c("", sprintf("Metadata file columns: %s", paste(names(meta), collapse = ", ")),
                                       sprintf("Metadata first row: %s", paste(unlist(meta[1]), collapse = " | "))),
                 sprintf("GEO record fields: %s", paste(unique(vapply(fields, `[`, "", 1)), collapse = ", ")),
                 sprintf("GEO first sample: %s", paste(vapply(fields, `[`, "", 2), collapse = " | ")))
write(link_report, "results/03_controls_file_check.txt", append = TRUE)
cat(link_report, sep = "\n")
if (link$method == "none" || mean(!is.na(col_gsm)) < 0.9) {
  stop("Could not link the columns to samples. Send me the text printed above ",
       "(or the file results/03_controls_file_check.txt).")
}

# ---------------------------------------------- 5-6. read and compute ----
# The full matrix (about 866,000 CpGs x 282 samples) can exhaust Colab's
# memory, so samples are read and scored in batches of 30. Each clock is a
# per-sample weighted sum, so batching does not change the results.
rm(peek); invisible(gc())
run_clocks <- function(m) {
  out <- tryCatch(DNAmAge(m, clocks = clocks, cell.count = FALSE),
                  error = function(e) DNAmAge(m, clocks = clocks))
  out <- as.data.frame(out)
  id_name <- intersect(c("id", "ID", "sample", "Sample", "sample_id"), names(out))[1]
  if (is.na(id_name)) out$id <- colnames(m) else names(out)[names(out) == id_name] <- "id"
  out
}
keep <- which(!is.na(col_gsm))
batches <- split(keep, ceiling(seq_along(keep) / 30))
results <- vector("list", length(batches))
for (i in seq_along(batches)) {
  idx <- batches[[i]]
  cat(sprintf("Batch %d of %d: reading %d samples and computing clocks...\n", i, length(batches), length(idx)))
  dt <- fread(proc_txt, select = c(id_col, beta_col[idx]), showProgress = FALSE)
  m <- as.matrix(dt[, -1, with = FALSE]); rownames(m) <- dt[[1]]; rm(dt)
  colnames(m) <- col_gsm[idx]
  if (!is_beta) m <- 2^m / (2^m + 1)
  results[[i]] <- run_clocks(m)
  if (i == 1) cat(sprintf("  (%d CpGs per sample%s)\n", nrow(m), if (!is_beta) "; M-values converted to beta" else ""))
  rm(m); invisible(gc())
}
clock_out <- as.data.frame(rbindlist(results, fill = TRUE))
pick <- function(pattern) {
  hit <- grep(pattern, names(clock_out), ignore.case = TRUE, value = TRUE)
  if (length(hit)) as.numeric(clock_out[[hit[1]]]) else NA_real_
}
ages <- data.frame(sample_id   = clock_out$id,
                   age_horvath = pick("^horvath$"),
                   age_hannum  = pick("^hannum$"),
                   age_pheno   = pick("^levine$|pheno"))
ages <- merge(samples, ages, by = "sample_id")
fwrite(ages, "data/processed/03_controls_clock_ages.csv")

# ---------------------------------------------- 7. how well do they fit --
long <- rbindlist(lapply(c(Horvath = "age_horvath", Hannum = "age_hannum", PhenoAge = "age_pheno"),
  function(col) data.table(age = ages$age, clock_age = ages[[col]])), idcol = "clock")
long$clock <- factor(long$clock, c("PhenoAge", "Horvath", "Hannum"))
stats <- long[!is.na(clock_age), .(n = .N,
                                   r = cor(age, clock_age),
                                   median_abs_error = median(abs(clock_age - age)),
                                   mean_difference = mean(clock_age - age)), by = clock]
summary_lines <- c(
  sprintf("Controls with clock estimates: %d of %d", nrow(ages), nrow(samples)),
  "",
  sprintf("%-9s n = %d | correlation with age r = %.2f | median absolute error %.1f years | mean clock minus age %+.1f years",
          stats$clock, stats$n, stats$r, stats$median_abs_error, stats$mean_difference),
  "",
  "Expected: r roughly 0.7-0.95 for all three clocks. A constant offset (mean",
  "difference) is normal and is removed when age acceleration is calculated."
)
writeLines(summary_lines, "results/03_controls_clock_summary.txt")
cat(summary_lines, sep = "\n")

p <- ggplot(long, aes(age, clock_age)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.55, size = 1.4, colour = "#5E4391") +
  geom_smooth(method = "lm", se = FALSE, colour = "#2E2B30", linewidth = 0.7) +
  facet_wrap(~clock) +
  labs(x = "Chronological age (years)", y = "Clock age (years)",
       title = "Epigenetic clocks in 282 community controls (SJLIFE, GSE197676)",
       subtitle = "Dashed line: clock age = chronological age") +
  theme_minimal(base_size = 12)
ggsave("results/03_controls_clock_vs_age.png", p, width = 10, height = 4, dpi = 200)
cat("\nDone. Download the three files in results/ and upload them to GitHub.\n")
