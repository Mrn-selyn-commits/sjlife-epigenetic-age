# 01_geo_inventory.R
# Purpose: find out what the SJLIFE methylation series on GEO actually contain
#          BEFORE downloading anything large.
#
# For each series it:
#   1. lists every file (series matrix and supplementary files) with its size,
#   2. reads only the header of each series matrix to get the sample
#      characteristics (age, sex, diagnosis, treatment ...), without
#      downloading the methylation values.
#
# Needs only base R and an internet connection. Uses a few MB of data.
# Outputs (small, safe to commit):
#   results/01_file_inventory.csv
#   results/01_<series>_characteristics_summary.txt
# Outputs (sample-level, kept out of GitHub):
#   data/processed/01_<series>_samples.csv

series <- c(
  controls = "GSE197676",  # community controls, no cancer history
  sjlife1  = "GSE197674",  # adult survivors
  sjlife2  = "GSE197675"   # younger survivors
)
superseries <- "GSE197678"  # links the three; checked for extra files only

options(timeout = 600)
dir.create("results", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- helpers ---

geo_dir <- function(gse, sub) {
  stub <- sub("\\d{1,3}$", "nnn", gse)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/", stub, gse, sub)
}

# File names listed in a GEO directory page (empty if the folder does not exist)
list_dir <- function(url) {
  page <- tryCatch(readLines(url, warn = FALSE), error = function(e) character())
  hrefs <- regmatches(page, regexpr('href="[^"?/][^"]*"', page))
  files <- gsub('^href="|"$', "", hrefs)
  files[!grepl("^https?:|/$", files)]
}

# Size in MB from the HTTP headers, without downloading the file
size_mb <- function(url) {
  h <- tryCatch(curlGetHeaders(url), error = function(e) character())
  starts <- grep("^HTTP/", h)                      # one block per redirect
  if (length(starts)) h <- h[tail(starts, 1):length(h)]
  len <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  if (!length(len)) return(NA_real_)
  signif(as.numeric(sub("^[^:]+:\\s*", "", tail(len, 1))) / 1024^2, 3)
}

# Sample rows of a series matrix header, read until the data table starts.
# `src` can be a URL or a local .gz file.
read_matrix_header <- function(src) {
  con <- if (grepl("^https?:", src)) gzcon(url(src, "rb")) else gzfile(src, "r")
  on.exit(close(con))
  out <- character()
  repeat {
    chunk <- readLines(con, n = 500, warn = FALSE)
    if (!length(chunk)) break
    stop_at <- grep("^!series_matrix_table_begin", chunk)
    if (length(stop_at)) {
      out <- c(out, chunk[seq_len(stop_at[1] - 1)])
      break
    }
    out <- c(out, chunk)
  }
  out[startsWith(out, "!Sample_")]
}

# Turn the "!Sample_..." rows into one row per sample.
# Each "characteristics" row holds "key: value"; keys become columns.
parse_samples <- function(rows) {
  split_row <- function(r) gsub('^"|"$', "", strsplit(r, "\t", fixed = TRUE)[[1]])
  fields <- lapply(rows, split_row)
  keys <- vapply(fields, `[`, "", 1)
  n <- length(fields[[1]]) - 1
  pick <- function(key) {
    i <- match(key, keys)
    if (is.na(i)) rep(NA_character_, n) else fields[[i]][-1]
  }
  df <- data.frame(sample_id = pick("!Sample_geo_accession"),
                   title     = pick("!Sample_title"),
                   source    = pick("!Sample_source_name_ch1"),
                   stringsAsFactors = FALSE)
  for (i in which(keys == "!Sample_characteristics_ch1")) {
    vals <- fields[[i]][-1]
    key <- sub(":.*$", "", vals[nzchar(vals)][1])
    col <- make.names(trimws(key))
    has <- grepl(":", vals)
    df[[col]] <- ifelse(has, trimws(sub("^[^:]*:", "", vals)), NA)
  }
  df
}

# A short, human-readable summary of each characteristic column
summarise_columns <- function(df) {
  lines <- c(sprintf("Samples: %d", nrow(df)), "")
  for (col in setdiff(names(df), c("sample_id", "title"))) {
    x <- df[[col]]
    num <- suppressWarnings(as.numeric(x))
    lines <- c(lines, sprintf("== %s (missing: %d)", col, sum(is.na(x) | x == "")))
    if (mean(!is.na(num)) > 0.9) {
      q <- quantile(num, c(0, .25, .5, .75, 1), na.rm = TRUE)
      lines <- c(lines, sprintf("  numeric: min %.1f | Q1 %.1f | median %.1f | Q3 %.1f | max %.1f",
                                q[1], q[2], q[3], q[4], q[5]))
    } else {
      tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
      shown <- head(tab, 15)
      lines <- c(lines, sprintf("  %s: %d", names(shown), as.integer(shown)))
      if (length(tab) > 15) lines <- c(lines, sprintf("  ... %d more values", length(tab) - 15))
    }
    lines <- c(lines, "")
  }
  lines
}

# ------------------------------------------------------- 1. file inventory ---

inventory <- list()
for (gse in c(series, superseries)) {
  for (sub in c("matrix", "suppl")) {
    base <- geo_dir(gse, sub)
    for (f in list_dir(base)) {
      inventory[[length(inventory) + 1]] <- data.frame(
        series = gse, folder = sub, file = f,
        size_mb = size_mb(paste0(base, f)), url = paste0(base, f))
    }
  }
}
inventory <- do.call(rbind, inventory)
if (is.null(inventory)) stop("No files found. Check the internet connection and try again.")
write.csv(inventory, "results/01_file_inventory.csv", row.names = FALSE)
cat("\nFILES ON GEO\n")
print(inventory[, c("series", "folder", "file", "size_mb")], row.names = FALSE)

# ------------------------------------------- 2. sample characteristics ------

for (set in names(series)) {
  gse <- series[[set]]
  mats <- inventory[inventory$series == gse & inventory$folder == "matrix" &
                    grepl("series_matrix\\.txt\\.gz$", inventory$file), ]
  if (!nrow(mats)) { message(gse, ": no series matrix found"); next }
  samples <- do.call(rbind, lapply(mats$url, function(u) parse_samples(read_matrix_header(u))))
  write.csv(samples, sprintf("data/processed/01_%s_samples.csv", gse), row.names = FALSE)
  summary_lines <- c(sprintf("%s (%s)", gse, set), summarise_columns(samples))
  writeLines(summary_lines, sprintf("results/01_%s_characteristics_summary.txt", gse))
  cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
}

cat("\nDone. Commit the files in results/ and send them to me for the next step.\n")
