# geo_helpers.R
# Shared helpers for downloading GEO files and computing the clocks.

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

# Read only the sample rows of a series matrix (a few KB).
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

# Run methylclock on a CpG x sample beta matrix; returns id + one column per clock.
run_clocks <- function(m, clocks = c("Horvath", "Hannum", "Levine")) {
  out <- tryCatch(methylclock::DNAmAge(m, clocks = clocks, cell.count = FALSE),
                  error = function(e) methylclock::DNAmAge(m, clocks = clocks))
  out <- as.data.frame(out)
  id_name <- intersect(c("id", "ID", "sample", "Sample", "sample_id"), names(out))[1]
  if (is.na(id_name)) out$id <- colnames(m) else names(out)[names(out) == id_name] <- "id"
  pick <- function(pattern) {
    hit <- grep(pattern, names(out), ignore.case = TRUE, value = TRUE)
    if (length(hit)) as.numeric(out[[hit[1]]]) else NA_real_
  }
  data.frame(sample_id   = out$id,
             age_horvath = pick("^horvath$"),
             age_hannum  = pick("^hannum$"),
             age_pheno   = pick("^levine$|pheno"))
}

# Clocks for selected columns of a tab-separated methylation file, in batches
# so memory stays low. `cols` are column names; `ids` the sample IDs to use.
clocks_in_batches <- function(path, id_col, cols, ids, is_beta = TRUE, size = 30) {
  batches <- split(seq_along(cols), ceiling(seq_along(cols) / size))
  res <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    idx <- batches[[i]]
    cat(sprintf("Batch %d of %d: %d samples...\n", i, length(batches), length(idx)))
    dt <- data.table::fread(path, select = c(id_col, cols[idx]), showProgress = FALSE)
    m <- as.matrix(dt[, -1, with = FALSE]); rownames(m) <- dt[[1]]; rm(dt)
    colnames(m) <- ids[idx]
    if (!is_beta) m <- 2^m / (2^m + 1)
    res[[i]] <- run_clocks(m)
    rm(m); invisible(gc())
  }
  do.call(rbind, res)
}
