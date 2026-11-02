# 02_match_survivors.R
# Chooses the survivors to compare with the 282 controls, using only age and
# sex, BEFORE any methylation data are opened (see analysis_plan.md).
#
# Method: 1:1 nearest-age matching within sex, without replacement.
#   - survivors only from SJLIFE1 (GSE197674), aged 18 or over at sampling
#   - controls are taken in a random order (fixed seed) so no one is favoured
#   - a survivor is accepted only if within `caliper` years of the control
#
# Needs (from script 01):
#   data/processed/01_GSE197676_samples.csv   controls
#   data/processed/01_GSE197674_samples.csv   SJLIFE1 survivors
# Outputs:
#   data/processed/02_matched_survivors.csv    survivors with treatment data
#   results/02_selected_survivors.txt          GSM accessions chosen (commit this)
#   results/02_matching_summary.txt            balance table (commit this)

seed    <- 20261027
caliper <- 2      # years

ctrl <- read.csv("data/processed/01_GSE197676_samples.csv", stringsAsFactors = FALSE)
surv <- read.csv("data/processed/01_GSE197674_samples.csv", stringsAsFactors = FALSE)
for (d in list(ctrl, surv)) {
  missing <- setdiff(c("sample_id", "gender", "age"), names(d))
  if (length(missing)) stop("Column(s) not found: ", paste(missing, collapse = ", "),
                            ". Rerun scripts/01_geo_inventory.R first.")
}
ctrl$age <- as.numeric(ctrl$age)
surv$age <- as.numeric(surv$age)
pool <- surv[!is.na(surv$age) & surv$age >= 18, ]

# ------------------------------------------------------------- matching ---
set.seed(seed)
taken <- rep(FALSE, nrow(pool))
pairs <- vector("list", nrow(ctrl))
for (i in sample(nrow(ctrl))) {
  cand <- which(!taken & pool$gender == ctrl$gender[i])
  if (!length(cand)) next
  gap  <- abs(pool$age[cand] - ctrl$age[i])
  best <- cand[gap == min(gap)]
  if (length(best) > 1) best <- best[sample.int(length(best), 1)]
  if (abs(pool$age[best] - ctrl$age[i]) > caliper) next
  taken[best] <- TRUE
  pairs[[i]] <- data.frame(control_id = ctrl$sample_id[i],
                           survivor_id = pool$sample_id[best],
                           age_gap = pool$age[best] - ctrl$age[i])
}
pairs <- do.call(rbind, pairs)
matched <- pool[match(pairs$survivor_id, pool$sample_id), ]
matched$matched_control <- pairs$control_id

dir.create("results", showWarnings = FALSE)
write.csv(matched, "data/processed/02_matched_survivors.csv", row.names = FALSE)
writeLines(sort(matched$sample_id), "results/02_selected_survivors.txt")

# -------------------------------------------------------------- balance ---
describe <- function(d, label) {
  sprintf("%-22s n = %4d | age mean %.1f (SD %.1f), median %.1f, range %.1f-%.1f | female %d (%.0f%%)",
          label, nrow(d), mean(d$age), sd(d$age), median(d$age), min(d$age), max(d$age),
          sum(d$gender == "Female"), 100 * mean(d$gender == "Female"))
}
smd <- (mean(matched$age) - mean(ctrl$age)) / sqrt((var(matched$age) + var(ctrl$age)) / 2)
unmatched <- ctrl[!ctrl$sample_id %in% pairs$control_id, ]

treat_cols <- intersect(c("chestrt", "abdominal_pelvic_rt", "brainrt", "alkylating.agent.classic",
                          "anthracyclines", "epipodophyllotoxins", "platinum", "corticosteroids",
                          "vincristine"), names(matched))
treat_lines <- vapply(treat_cols, function(col) {
  x <- as.numeric(matched[[col]])
  sprintf("  %-26s %4d of %d (%.0f%%)", col, sum(x == 1, na.rm = TRUE), sum(!is.na(x)),
          100 * mean(x == 1, na.rm = TRUE))
}, "")

summary_lines <- c(
  sprintf("Matching run: seed %d, caliper %g years, 1:1 within sex, survivors aged >= 18", seed, caliper),
  "",
  describe(ctrl, "Controls (all)"),
  describe(pool, "Survivor pool (18+)"),
  describe(matched, "Matched survivors"),
  "",
  sprintf("Standardised mean difference in age, matched survivors vs controls: %.3f (below 0.1 is good balance)", smd),
  sprintf("Mean absolute age gap within pairs: %.2f years", mean(abs(pairs$age_gap))),
  sprintf("Controls with no survivor within %g years: %d%s", caliper, nrow(unmatched),
          if (nrow(unmatched)) sprintf(" (ages %s)", paste(sprintf("%.1f", sort(unmatched$age)), collapse = ", ")) else ""),
  "",
  "Treatment exposure among matched survivors (1 = received):",
  treat_lines
)
writeLines(summary_lines, "results/02_matching_summary.txt")
cat(summary_lines, sep = "\n")
cat("\nDone. Commit results/02_selected_survivors.txt and results/02_matching_summary.txt",
    "BEFORE running any clock script: the commit time is the record that the",
    "survivors were chosen before their methylation data were seen.\n", sep = "\n")
