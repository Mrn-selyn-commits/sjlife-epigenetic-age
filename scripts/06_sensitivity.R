# 06_sensitivity.R
# EXPLORATORY / POST HOC. Added on 3 November 2026, after the primary results
# (script 05) had been seen, in response to external review. None of these
# analyses replaces the pre-specified results in results/05_results.txt.
#
# Runs on any computer with R. Needs the same two local files as script 05:
#   data/processed/03_controls_clock_ages.csv
#   data/processed/04_survivors_clock_ages.csv
#
# What it checks:
#   A. Group difference estimated directly: clock age ~ age + group + sex
#      (no residuals), and with the age line fitted in controls only
#   B. Whether the group difference changes with age: clock age ~ age * group + sex
#   C. Effect sizes in standard-deviation units (SD of age acceleration in controls)
#   D. Holm correction across the five secondary tests
#   E. Smallest difference the design could detect (80% power, alpha 0.05)
#   F. Chest radiotherapy adjusted for age, sex and every other treatment flag in GEO
#
# Outputs (commit these):
#   results/06_sensitivity.txt, results/06_sensitivity_table.csv

suppressPackageStartupMessages(library(data.table))
dir.create("results", showWarnings = FALSE)

ctrl <- fread("data/processed/03_controls_clock_ages.csv")[, group := "control"]
surv <- fread("data/processed/04_survivors_clock_ages.csv")[, group := "survivor"]
d <- rbindlist(list(ctrl, surv), fill = TRUE)
d[, sex := factor(tolower(substr(trimws(sex), 1, 1)), levels = c("f", "m"), labels = c("female", "male"))]
d[, group := factor(group, levels = c("control", "survivor"))]
d[, age := as.numeric(age)]
clocks <- c(PhenoAge = "age_pheno", Horvath = "age_horvath", Hannum = "age_hannum")

est <- function(fit, term) {
  co <- summary(fit)$coefficients; ci <- confint(fit)
  list(est = co[term, 1], lo = ci[term, 1], hi = ci[term, 2], p = co[term, 4], n = nobs(fit))
}
rows <- list()
add <- function(section, clock, model, e, sd_units = NA_real_) {
  rows[[length(rows) + 1]] <<- data.table(section = section, clock = clock, model = model, n = e$n,
                                          estimate = e$est, ci_low = e$lo, ci_high = e$hi,
                                          p_value = e$p, estimate_sd_units = sd_units)
}

for (k in names(clocks)) {
  col <- clocks[[k]]
  # Pre-specified residual (as in script 05), kept for the SD and power calculations
  d[, eaa := residuals(lm(as.formula(paste(col, "~ age")), data = d, na.action = na.exclude))]
  sd_ctrl <- sd(d[group == "control", eaa], na.rm = TRUE)

  # A1. direct model, no residuals
  f1 <- lm(as.formula(paste(col, "~ age + group + sex")), data = d)
  e1 <- est(f1, "groupsurvivor")
  add("A. Direct model", k, "clock ~ age + group + sex", e1, e1$est / sd_ctrl)

  # A2. age line fitted in controls only, applied to everyone
  ref <- lm(as.formula(paste(col, "~ age")), data = d[group == "control"])
  d[, eaa_ctrl := get(col) - predict(ref, newdata = d)]
  f2 <- lm(eaa_ctrl ~ group + sex, data = d)
  e2 <- est(f2, "groupsurvivor")
  add("A. Controls-only age line", k, "EAA(controls' line) ~ group + sex", e2, e2$est / sd_ctrl)

  # B. does the gap change with age? (difference in slope, years of clock age per year of age)
  f3 <- lm(as.formula(paste(col, "~ age * group + sex")), data = d)
  add("B. Group x age", k, "clock ~ age * group + sex (interaction term)", est(f3, "age:groupsurvivor"))

  # C/E. pre-specified estimate in SD units, and the smallest detectable difference
  f0 <- lm(eaa ~ group + sex, data = d)
  e0 <- est(f0, "groupsurvivor")
  add("C. Pre-specified, SD units", k, "EAA ~ group + sex (as script 05)", e0, e0$est / sd_ctrl)
  mdd <- power.t.test(n = min(table(d$group[!is.na(d$eaa)])), sd = sd_ctrl, sig.level = 0.05, power = 0.8)$delta
  rows[[length(rows) + 1]] <- data.table(section = "E. Detectable difference", clock = k,
                                         model = "80% power, alpha 0.05, two-sided, SD of EAA in controls",
                                         n = nobs(f0), estimate = mdd, ci_low = NA, ci_high = NA,
                                         p_value = NA, estimate_sd_units = mdd / sd_ctrl)
}

# F. chest radiotherapy adjusted for age, sex and the other treatment flags
s <- d[group == "survivor"]
flag_cols <- names(s)[vapply(s, function(v) {
  x <- suppressWarnings(as.numeric(v)); ok <- !is.na(x)
  is.numeric(v) && any(ok) && all(x[ok] %in% c(0, 1))
}, logical(1))]
flag_cols <- setdiff(flag_cols, c("eaa", "eaa_ctrl"))
rt_col <- grep("^chest", flag_cols, value = TRUE)[1]
others <- setdiff(flag_cols, rt_col)
if (!is.na(rt_col)) {
  for (k in names(clocks)) {
    col <- clocks[[k]]
    d[, eaa := residuals(lm(as.formula(paste(col, "~ age")), data = d, na.action = na.exclude))]
    s <- d[group == "survivor"]
    form <- paste("eaa ~", paste(c(rt_col, "sex", "age", others), collapse = " + "))
    f <- lm(as.formula(form), data = s)
    add("F. Chest RT, fully adjusted", k, paste("EAA ~", rt_col, "+ sex + age + all other treatment flags"), est(f, rt_col))
  }
}

res <- rbindlist(rows, fill = TRUE)

# D. Holm correction across the five pre-specified secondary tests (exact p-values from script 05)
sec <- tryCatch(fread("results/05_results_table.csv")[grepl("^Secondary", analysis)], error = function(e) NULL)
holm_lines <- character()
if (!is.null(sec) && nrow(sec)) {
  sec[, p_holm := p.adjust(p_value, method = "holm")]
  holm_lines <- c("D. Holm correction across the five secondary tests (from script 05)",
                  sprintf("  %-46s %-9s p = %.2g -> Holm p = %.2g", sec$analysis, sec$clock, sec$p_value, sec$p_holm), "")
}

fwrite(res, "results/06_sensitivity_table.csv")
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
lines <- c("EXPLORATORY / POST HOC analyses, added 3 November 2026 after the primary results were seen.",
           "They do not replace the pre-specified results in results/05_results.txt.", "")
for (sec_name in unique(res$section)) {
  lines <- c(lines, sec_name)
  r <- res[section == sec_name]
  if (startsWith(sec_name, "E.")) {
    lines <- c(lines, sprintf("  %-9s smallest detectable difference %.2f years (%.2f SD of EAA in controls), n per group %d",
                              r$clock, r$estimate, r$estimate_sd_units, floor(r$n / 2)))
  } else if (startsWith(sec_name, "B.")) {
    lines <- c(lines, sprintf("  %-9s difference in slope %+.3f years per year of age (95%% CI %+.3f to %+.3f) | %s",
                              r$clock, r$estimate, r$ci_low, r$ci_high, fmt_p(r$p_value)))
  } else {
    lines <- c(lines, sprintf("  %-9s %+.2f years (95%% CI %+.2f to %+.2f) | %s%s | n = %d",
                              r$clock, r$estimate, r$ci_low, r$ci_high, fmt_p(r$p_value),
                              ifelse(is.na(r$estimate_sd_units), "", sprintf(" | %.2f SD", r$estimate_sd_units)), r$n))
  }
  lines <- c(lines, "")
}
lines <- c(lines, holm_lines)
if (!is.na(rt_col)) lines <- c(lines, sprintf("Treatment flags used in F: %s (+ sex, age)", paste(c(rt_col, others), collapse = ", ")))
writeLines(lines, "results/06_sensitivity.txt")
cat(lines, sep = "\n")
cat("\nDone. Upload results/06_sensitivity.txt and results/06_sensitivity_table.csv to GitHub.\n")
