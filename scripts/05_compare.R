# 05_compare.R
# Phase 3 of the analysis plan, part 2: the pre-specified comparison.
# Runs on any computer with R (no large downloads). Needs two files made by
# scripts 03 and 04 (keep them in data/processed/, not on GitHub):
#   data/processed/03_controls_clock_ages.csv
#   data/processed/04_survivors_clock_ages.csv
#
# Analysis plan (analysis_plan.md):
#   Outcome   epigenetic age acceleration (EAA) = residual from regressing
#             clock age on chronological age across the whole sample
#   Primary   lm(EAA_PhenoAge ~ group + sex): survivors vs controls
#   Secondary the same for Horvath and Hannum; within survivors,
#             lm(EAA ~ chest radiotherapy + sex)
#
# Outputs (commit these):
#   results/05_results.txt, results/05_results_table.csv
#   results/05_fig1_clock_vs_age.png, results/05_fig2_eaa_by_group.png,
#   results/05_fig3_eaa_by_chest_rt.png

for (p in c("data.table", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("results", showWarnings = FALSE)

ctrl <- fread("data/processed/03_controls_clock_ages.csv")[, group := "control"]
surv <- fread("data/processed/04_survivors_clock_ages.csv")[, group := "survivor"]
d <- rbindlist(list(ctrl, surv), fill = TRUE)
d[, sex := factor(tolower(substr(trimws(sex), 1, 1)), levels = c("f", "m"), labels = c("female", "male"))]
d[, group := factor(group, levels = c("control", "survivor"))]
d[, age := as.numeric(age)]

clocks <- c(PhenoAge = "age_pheno", Horvath = "age_horvath", Hannum = "age_hannum")
for (k in names(clocks)) {
  col <- clocks[[k]]; eaa <- paste0("eaa_", k)
  fit <- lm(as.formula(paste(col, "~ age")), data = d, na.action = na.exclude)
  d[, (eaa) := as.numeric(residuals(fit))]
}

tidy <- function(fit, term, label, analysis) {
  co <- summary(fit)$coefficients; ci <- confint(fit)
  data.table(analysis = analysis, clock = label, n = nobs(fit),
             difference_years = co[term, 1], ci_low = ci[term, 1], ci_high = ci[term, 2],
             p_value = co[term, 4])
}
rows <- list()
for (k in names(clocks)) {
  f <- lm(as.formula(paste0("eaa_", k, " ~ group + sex")), data = d)
  rows[[length(rows) + 1]] <- tidy(f, "groupsurvivor", k,
                                   if (k == "PhenoAge") "Primary: survivors vs controls" else "Secondary: survivors vs controls")
}
s <- d[group == "survivor"]
has_rt <- "chestrt" %in% names(s)
if (has_rt) {
  s[, chest_rt := factor(as.integer(as.character(chestrt)), levels = c(0, 1), labels = c("no", "yes"))]
  for (k in names(clocks)) {
    f <- lm(as.formula(paste0("eaa_", k, " ~ chest_rt + sex")), data = s)
    rows[[length(rows) + 1]] <- tidy(f, "chest_rtyes", k, "Secondary: chest radiotherapy (survivors only)")
  }
}
res <- rbindlist(rows)
fwrite(res, "results/05_results_table.csv")

fmt <- function(r) sprintf("  %-9s n = %3d | difference %+.2f years (95%% CI %+.2f to %+.2f) | p %s",
                           r$clock, r$n, r$difference_years, r$ci_low, r$ci_high,
                           ifelse(r$p_value < 0.001, "< 0.001", sprintf("= %.3f", r$p_value)))
lines <- c(sprintf("Controls: %d | Matched survivors: %d", sum(d$group == "control"), sum(d$group == "survivor")),
           "Epigenetic age acceleration = residual of clock age on chronological age (whole sample).",
           "Positive difference = survivors (or those with chest radiotherapy) biologically older.", "")
for (a in unique(res$analysis)) {
  lines <- c(lines, a)
  for (i in which(res$analysis == a)) lines <- c(lines, fmt(res[i]))
  lines <- c(lines, "")
}
if (has_rt) lines <- c(lines, sprintf("Chest radiotherapy among survivors: %d yes, %d no", sum(s$chest_rt == "yes", na.rm = TRUE), sum(s$chest_rt == "no", na.rm = TRUE)))
lines <- c(lines, "Significance level 0.05 (analysis plan). Analyses not in the plan are labelled exploratory elsewhere.")
writeLines(lines, "results/05_results.txt")
cat(lines, sep = "\n")

# ---- figures (lilac for survivors, slate for controls)
pal <- c(control = "#3B3049", survivor = "#7A5BB0")
p1 <- ggplot(d, aes(age, age_pheno, colour = group)) +
  geom_point(alpha = 0.55, size = 1.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Chronological age (years)", y = "PhenoAge (years)",
       title = "PhenoAge against chronological age") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("results/05_fig1_clock_vs_age.png", p1, width = 6, height = 5, dpi = 200)

long <- melt(d, id.vars = c("sample_id", "group"), measure.vars = paste0("eaa_", names(clocks)),
             variable.name = "clock", value.name = "eaa")
long[, clock := factor(sub("eaa_", "", clock), levels = names(clocks))]
p2 <- ggplot(long, aes(group, eaa, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.35) +
  facet_wrap(~clock) + scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Epigenetic age acceleration (years)", title = "Age acceleration: survivors vs controls") +
  theme_minimal(base_size = 12)
ggsave("results/05_fig2_eaa_by_group.png", p2, width = 9, height = 4.5, dpi = 200)

if (has_rt) {
  p3 <- ggplot(s[!is.na(chest_rt)], aes(chest_rt, eaa_PhenoAge)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_boxplot(width = 0.5, fill = "#CDBDE3", outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.4) +
    labs(x = "Chest radiotherapy", y = "PhenoAge acceleration (years)", title = "Survivors: chest radiotherapy") +
    theme_minimal(base_size = 12)
  ggsave("results/05_fig3_eaa_by_chest_rt.png", p3, width = 4.5, height = 4.5, dpi = 200)
}
cat("\nDone. Upload the results/05_* files to GitHub.\n")
