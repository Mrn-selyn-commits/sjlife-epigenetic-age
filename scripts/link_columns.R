# link_columns.R
# Links the sample columns of a GEO "processed" methylation file to GEO sample
# accessions (GSM). Used by the clock scripts. Tries, in order:
#   A. the column name appears in the sample's GEO record (title, description...)
#      - exactly, or after removing suffixes such as "AVG_Beta" and punctuation;
#   B. the GEO supplementary metadata file names the columns, and its rows are
#      tied to GSMs through a shared field;
#   C. the metadata row is tied to a GSM by a fingerprint of its recorded
#      characteristics (age, sex, and treatment flags where present). Only
#      fingerprints that are unique on both sides are used.
# Returns a list: gsm (one per column, NA if unlinked), method, report lines.

norm_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("[._ -]*(red|grn)?\\.idat(\\.gz)?$", "", x)          # 2034..._R03C01_Red.idat
  x <- sub("[._ -]*(avg[._ -]*)?beta$", "", x)
  gsub("[^a-z0-9]", "", x)
}

norm_value <- function(x) {
  x <- trimws(tolower(as.character(x)))
  num <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(num), sprintf("%.4f", num), x)
}

# Characteristics from series-matrix rows: list of "key: value" vectors -> data.frame
series_characteristics <- function(fields, gsm) {
  keys <- vapply(fields, `[`, "", 1)
  out <- data.frame(sample_id = gsm, stringsAsFactors = FALSE)
  for (i in which(keys == "!Sample_characteristics_ch1")) {
    v <- fields[[i]][-1]
    k <- sub(":.*$", "", v[grepl(":", v)][1])
    if (is.na(k)) next
    out[[norm_name(k)]] <- ifelse(grepl(":", v), trimws(sub("^[^:]*:", "", v)), NA)
  }
  out
}

link_columns <- function(cols, fields, gsm, meta = NULL) {
  report <- character()
  # A. direct match against every value in the GEO sample rows
  exact <- norm <- character()
  for (f in fields) for (j in seq_along(f)[-1]) {
    exact[trimws(f[j])] <- gsm[j - 1]
    norm[norm_name(f[j])] <- gsm[j - 1]
    # IDs written inside square brackets, e.g. "... [SJNORM044990_G1]"
    br <- regmatches(f[j], gregexpr("\\[([^]]+)\\]", f[j]))[[1]]
    for (b in gsub("^\\[|\\]$", "", br)) {
      exact[trimws(b)] <- gsm[j - 1]
      norm[norm_name(b)] <- gsm[j - 1]
    }
  }
  exact[gsm] <- gsm
  a <- unname(exact[trimws(cols)])
  a[is.na(a)] <- unname(norm[norm_name(cols)])[is.na(a)]
  report <- c(report, sprintf("A. column names found in GEO sample records: %d of %d", sum(!is.na(a)), length(cols)))
  if (mean(!is.na(a)) >= 0.9) return(list(gsm = a, method = "A: column names in GEO sample records", report = report))

  if (is.null(meta) || !nrow(meta)) {
    report <- c(report, "B/C. no metadata file available")
    return(list(gsm = a, method = "none", report = report))
  }
  meta <- as.data.frame(meta)
  # which metadata column holds the processed-file column names?
  hit <- vapply(meta, function(v) mean(norm_name(cols) %in% norm_name(v)), numeric(1))
  key_col <- names(hit)[which.max(hit)]
  report <- c(report, sprintf("B. metadata column naming the samples: '%s' (%.0f%% of columns found)", key_col, 100 * max(hit)))
  if (max(hit) < 0.9) return(list(gsm = a, method = "none", report = report))
  meta_row <- match(norm_name(cols), norm_name(meta[[key_col]]))

  # B. a metadata field equal to a GEO record value (title, GSM, ...)
  b_map <- rep(NA_character_, nrow(meta))
  for (cn in names(meta)) {
    g <- unname(exact[trimws(as.character(meta[[cn]]))])
    if (mean(!is.na(g)) >= 0.9 && !anyDuplicated(na.omit(g))) { b_map <- g; break }
  }
  b <- b_map[meta_row]
  report <- c(report, sprintf("B. linked through a shared metadata field: %d of %d", sum(!is.na(b)), length(cols)))
  if (mean(!is.na(b)) >= 0.9) return(list(gsm = b, method = "B: metadata file field shared with GEO records", report = report))

  # C. fingerprint of shared characteristics
  ch <- series_characteristics(fields, gsm)
  mnames <- setNames(names(meta), norm_name(sub("^characteristics[ _:.-]*", "", names(meta), ignore.case = TRUE)))
  alias <- c(gender = "sex", sex = "gender")
  shared <- character()
  for (k in setdiff(names(ch), "sample_id")) {
    mk <- if (k %in% names(mnames)) k else if (!is.na(alias[k]) && alias[k] %in% names(mnames)) alias[k] else NA
    if (!is.na(mk)) shared[k] <- mnames[[mk]]
  }
  report <- c(report, sprintf("C. characteristics shared with the metadata file: %s",
                              if (length(shared)) paste(names(shared), collapse = ", ") else "none"))
  if (length(shared) < 2) return(list(gsm = a, method = "none", report = report))
  fp_geo  <- do.call(paste, c(lapply(names(shared), function(k) norm_value(ch[[k]])), sep = "|"))
  fp_meta <- do.call(paste, c(lapply(shared, function(k) norm_value(meta[[k]])), sep = "|"))
  uniq_geo  <- fp_geo[!fp_geo %in% fp_geo[duplicated(fp_geo)]]
  uniq_meta <- !fp_meta %in% fp_meta[duplicated(fp_meta)]
  c_map <- ifelse(uniq_meta, gsm[match(fp_meta, ifelse(fp_geo %in% uniq_geo, fp_geo, NA))], NA)
  cc <- c_map[meta_row]
  report <- c(report, sprintf("C. linked by unique fingerprint: %d of %d (others had a duplicate fingerprint)", sum(!is.na(cc)), length(cols)))
  method <- if (mean(!is.na(cc)) >= 0.9) "C: unique fingerprint of age, sex and other recorded characteristics" else "none"
  list(gsm = cc, method = method, report = report)
}
