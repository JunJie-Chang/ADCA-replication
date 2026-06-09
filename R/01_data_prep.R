# =============================================================================
# 01_data_prep.R
# Phase 1: Read, clean, and reshape harrison.xls
#
# NOTE: Run from RStudio. The project folder contains Unicode characters
# that Rscript in terminal cannot encode.
#
# Output objects saved to output/:
#   df_wide      one row per respondent (demographics + scenario choices)
#   df_long      one row per respondent x scenario
#   df_mlogit    one row per respondent x scenario x alternative
#   price_matrix 17 scenarios x 5 steak price matrix
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)

# ── 0. Paths ──────────────────────────────────────────────────────────────────

project_dir <- dirname(
  dirname(rstudioapi::getActiveDocumentContext()$path)
)
data_path  <- file.path(project_dir, "harrison.xls")
output_dir <- file.path(project_dir, "output")

# ── 1. Column names ───────────────────────────────────────────────────────────

demo_names <- c(
  "subject",     "gender",     "age",         "household", "child",
  "education",   "income",     "student",      "occupation","race",
  "gbeef",       "steak_buy",  "other_beef",
  "gbeef_freq",  "steak_freq",
  "color",       "brand",      "quality_gd",   "ex_fat",    "in_fat",
  "price_imp",   "safety",     "juiciness",    "flavor",    "tenderness",
  "consistency", "doneness",   "production",   "quality_gd2",
  "slaughter",   "food_safe",  "quality_pref", "cooked",
  "chance1",     "chance2",    "chance3",      "chance4",   "chance5"
)
sc_names      <- paste0("sc", 1:17)
all_col_names <- c(demo_names, sc_names, "blank", "session")

# ── 2. Read raw data ──────────────────────────────────────────────────────────
# skip = 3: row 1 (variable labels), row 2 (Q-numbers), row 3 (blank)

raw <- read_excel(
  data_path,
  sheet     = "trt 10&11",
  col_names = all_col_names,
  skip      = 3,
  col_types = "text"
)

# ── 3. Filter valid respondent rows ───────────────────────────────────────────
# Keep only rows where subject is a numeric ID in 1001-1199.
# Drops: treatment label rows, blank rows, and summary rows at the bottom.
#
# Subject 1060 (Treatment 10) has all 17 scenarios missing — excluded.
# This reconciles the 68 subjects in the file vs. 67 used in the paper.

df_wide <- raw |>
  mutate(subject = suppressWarnings(as.numeric(subject))) |>
  filter(!is.na(subject), subject >= 1001, subject <= 1199) |>
  filter(subject != 1060) |>
  mutate(
    treatment = if_else(subject < 1100, "nonhyp", "hyp"),
    treatment = factor(treatment, levels = c("nonhyp", "hyp"))
  )

cat("nonhyp:", sum(df_wide$treatment == "nonhyp"), "(paper: 67)\n")
cat("hyp:   ", sum(df_wide$treatment == "hyp"),    "(paper: 37)\n")

# ── 4. Convert numeric columns ────────────────────────────────────────────────

numeric_cols <- c(
  "gender", "age", "household", "child", "education", "income",
  "student", "race", "gbeef", "steak_buy", "other_beef",
  "gbeef_freq", "steak_freq", "color", "brand", "quality_gd",
  "ex_fat", "in_fat", "price_imp", "safety", "juiciness", "flavor",
  "tenderness", "consistency", "doneness", "production",
  "quality_gd2", "slaughter", "food_safe", "quality_pref", "cooked",
  "chance1", "chance2", "chance3", "chance4", "chance5",
  sc_names, "session"
)

df_wide <- df_wide |>
  mutate(across(all_of(numeric_cols), as.numeric))

# ── 5. Price matrix (from paper Appendix) ─────────────────────────────────────
# 17 scenarios x 5 steak alternatives; prices in USD for 12-oz ribeye steaks

price_matrix <- tribble(
  ~scenario, ~GEN,  ~GT,   ~NAT,  ~CHO,  ~CAB,
  1,          3.38,  7.88,  6.75,  6.75,  6.75,
  2,          4.50,  9.00,  6.75,  5.63,  9.00,
  3,          6.75,  9.00,  5.63,  9.00,  6.75,
  4,          5.63,  9.00,  9.00,  6.75,  5.63,
  5,          6.75,  6.75,  6.75,  7.88,  7.88,
  6,          4.50,  7.88,  7.88,  9.00,  9.00,
  7,          4.50,  6.75,  5.63,  6.75,  7.88,
  8,          5.63,  5.63,  5.63,  7.88,  9.00,
  9,          5.63,  7.88,  5.63,  7.88,  9.00,
  10,         5.63,  6.75,  7.88,  5.63,  9.00,
  11,         6.75,  7.88,  9.00,  5.63,  7.88,
  12,         4.50,  5.63,  5.63,  7.88,  7.88,
  13,         3.38,  9.00,  7.88,  7.88,  9.00,
  14,         3.38,  5.63,  5.63,  5.63,  5.63,
  15,         3.38,  6.75,  9.00,  9.00,  9.00,
  16,         6.75,  5.63,  7.88,  6.75,  6.75,
  17,         5.63,  5.63,  5.63,  5.63,  5.63
)

# ── 6. Reshape wide → long (one row per respondent x scenario) ───────────────

df_long <- df_wide |>
  select(subject, treatment, session, all_of(sc_names)) |>
  pivot_longer(
    cols      = all_of(sc_names),
    names_to  = "sc_var",
    values_to = "choice_code"
  ) |>
  mutate(
    scenario = as.integer(sub("sc", "", sc_var)),
    choice   = factor(
      choice_code,
      levels = 1:6,
      labels = c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")
    )
  ) |>
  select(-sc_var, -choice_code) |>
  left_join(price_matrix, by = "scenario") |>
  arrange(subject, scenario)

# ── 7. Build mlogit dataset (one row per respondent x scenario x alt) ─────────
# Estimation uses scenarios 1-16 only.
# Scenario 17 (all prices equal at $5.63) is the binding scenario and excluded.

alt_labels <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")

df_mlogit <- df_long |>
  filter(scenario <= 16) |>
  crossing(alt = factor(alt_labels, levels = alt_labels)) |>
  mutate(
    chosen    = (as.character(choice) == as.character(alt)),
    price_alt = case_when(
      alt == "GEN"  ~ GEN,
      alt == "GT"   ~ GT,
      alt == "NAT"  ~ NAT,
      alt == "CHO"  ~ CHO,
      alt == "CAB"  ~ CAB,
      alt == "NONE" ~ NA_real_
    ),
    chid = paste(subject, scenario, sep = "_")
  ) |>
  select(chid, subject, treatment, session, scenario, alt, chosen, price_alt) |>
  arrange(subject, scenario, alt)

cat(
  "nonhyp obs:",
  sum(df_mlogit$treatment == "nonhyp" & df_mlogit$alt == "GEN"),
  "(paper: 1072)\n"
)
cat(
  "hyp obs:   ",
  sum(df_mlogit$treatment == "hyp" & df_mlogit$alt == "GEN"),
  "(paper: 592)\n"
)

# ── 8. Save ───────────────────────────────────────────────────────────────────

saveRDS(df_wide,      file.path(output_dir, "df_wide.rds"))
saveRDS(df_long,      file.path(output_dir, "df_long.rds"))
saveRDS(df_mlogit,    file.path(output_dir, "df_mlogit.rds"))
saveRDS(price_matrix, file.path(output_dir, "price_matrix.rds"))

cat("Saved: df_wide, df_long, df_mlogit, price_matrix -> output/\n")
