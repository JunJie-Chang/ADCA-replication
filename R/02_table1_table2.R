# =============================================================================
# 02_table1_table2.R
# Phase 1: Reproduce Table 1 (demographics) and Table 2 (choice frequencies)
#
# Requires: output/df_wide.rds and output/df_long.rds from 01_data_prep.R
# =============================================================================

library(dplyr)
library(tidyr)
library(knitr)

# Suppress R CMD CHECK notes for dplyr NSE column names used in this script.
utils::globalVariables(c(
  "treatment", "variable", "nonhyp", "hyp",
  "n", "choice", "scenario", "NONE", "diff"
))

# ── 0. Load data ──────────────────────────────────────────────────────────────

project_dir <- dirname(
  dirname(rstudioapi::getActiveDocumentContext()$path)
)
output_dir <- file.path(project_dir, "output")

df_wide <- readRDS(file.path(output_dir, "df_wide.rds"))
df_long <- readRDS(file.path(output_dir, "df_long.rds"))

# =============================================================================
# TABLE 1: Summary Statistics of Selected Participant Demographics
# =============================================================================
# Reproduces Table 1 in Lusk & Schroeder (2004, p. 474).
# Null of equal means across treatments cannot be rejected for any variable.

mean_sd <- function(data, var) {
  data |>
    group_by(.data[["treatment"]]) |>
    summarise(
      mean = round(mean(.data[[var]], na.rm = TRUE), 2),
      sd   = round(sd(.data[[var]],   na.rm = TRUE), 2),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from  = "treatment",
      values_from = c("mean", "sd")
    ) |>
    mutate(variable = var) |>
    select("variable", "mean_nonhyp", "sd_nonhyp", "mean_hyp", "sd_hyp")
}

table1_vars <- c("gender", "age", "education", "student", "income")

table1 <- bind_rows(lapply(table1_vars, mean_sd, data = df_wide))

n_counts <- df_wide |>
  count(.data[["treatment"]]) |>
  pivot_wider(names_from = "treatment", values_from = "n")

n_row <- tibble(
  variable    = "N",
  mean_nonhyp = n_counts[["nonhyp"]],
  mean_hyp    = n_counts[["hyp"]],
  sd_nonhyp   = NA_real_,
  sd_hyp      = NA_real_
)

var_labels <- c(
  gender    = "Gender (1=female; 0=male)",
  age       = "Age (years)",
  education = "Education (1=HS; 5=bachelor's; 8=PhD)",
  student   = "Student (1=yes; 0=no)",
  income    = "Household income level",
  N         = "N"
)

table1 <- bind_rows(table1, n_row) |>
  mutate(variable = var_labels[variable])

cat("=== TABLE 1 ===\n\n")
print(kable(
  table1,
  format    = "simple",
  digits    = 2,
  col.names = c("Variable", "Mean", "(SD)", "Mean", "(SD)"),
  caption   = "Table 1. Summary Statistics by Treatment"
))
cat("\n  Left columns: Nonhypothetical | Right columns: Hypothetical\n")

cat("\n--- t-tests: Ho = equal means across treatments ---\n")
for (v in table1_vars) {
  x  <- df_wide |> filter(.data[["treatment"]] == "nonhyp") |> pull(v)
  y  <- df_wide |> filter(.data[["treatment"]] == "hyp")    |> pull(v)
  pv <- t.test(x, y)$p.value
  verdict <- if (pv < 0.05) "* SIGNIFICANT" else "(cannot reject equality)"
  cat(sprintf("  %-12s p = %.3f  %s\n", v, pv, verdict))
}

write.csv(
  table1,
  file.path(output_dir, "table1_demographics.csv"),
  row.names = FALSE
)
cat("\nSaved: output/table1_demographics.csv\n")

# =============================================================================
# TABLE 2: Hypothetical and Nonhypothetical Choices in 17 Scenarios
# =============================================================================
# Reproduces Table 2 in Lusk & Schroeder (2004, p. 475).
# Cells show % of respondents choosing each steak type per scenario.

alt_order <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")

table2 <- df_long |>
  group_by(.data[["treatment"]], .data[["scenario"]], .data[["choice"]]) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(.data[["treatment"]], .data[["scenario"]]) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  select(-"n") |>
  pivot_wider(
    names_from  = "choice",
    values_from = "pct",
    values_fill = 0
  ) |>
  select("treatment", "scenario", all_of(alt_order)) |>
  arrange(.data[["treatment"]], .data[["scenario"]])

col_labels <- c(
  "Scenario", "Generic", "GT", "Natural", "USDA Choice", "CAB", "None"
)

cat("\n=== TABLE 2: Nonhypothetical Treatment ===\n\n")
print(kable(
  table2 |>
    filter(.data[["treatment"]] == "nonhyp") |>
    select(-"treatment"),
  format    = "simple",
  col.names = col_labels,
  caption   = "Table 2 (Nonhyp). % choosing each steak by scenario"
))

cat("\n=== TABLE 2: Hypothetical Treatment ===\n\n")
print(kable(
  table2 |>
    filter(.data[["treatment"]] == "hyp") |>
    select(-"treatment"),
  format    = "simple",
  col.names = col_labels,
  caption   = "Table 2 (Hyp). % choosing each steak by scenario"
))

# Verification: paper states NONE gap is smallest at scenario 4 (8.31%)
# and largest at scenario 7 (32.84%).
cat("\n--- NONE choice rate comparison (scenarios 1-16) ---\n")
none_check <- table2 |>
  filter(.data[["scenario"]] <= 16) |>
  select("treatment", "scenario", "NONE") |>
  pivot_wider(names_from = "treatment", values_from = "NONE") |>
  mutate(diff = .data[["nonhyp"]] - .data[["hyp"]])

cat(sprintf(
  "  Smallest NONE gap: scenario %d (diff = %.1f%%)\n",
  none_check$scenario[which.min(abs(none_check$diff))],
  min(abs(none_check$diff))
))
cat(sprintf(
  "  Largest NONE gap:  scenario %d (diff = %.1f%%)\n",
  none_check$scenario[which.max(none_check$diff)],
  max(none_check$diff)
))
cat("  Paper: smallest = scenario 4 (8.31%), largest = scenario 7 (32.84%)\n")

write.csv(
  table2,
  file.path(output_dir, "table2_choice_frequencies.csv"),
  row.names = FALSE
)
cat("\nSaved: output/table2_choice_frequencies.csv\n")
