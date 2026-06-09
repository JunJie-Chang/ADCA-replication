# =============================================================================
# _common.R — shared setup sourced at the top of every chapter.
# Provides: packages, knitr options, helper functions, the pre-estimated model
# objects (loaded from output/), and the analysis-ready data frames.
#
# Heavy estimation is NOT re-run when the book is built: the numbered scripts in
# R/ produce output/*.rds, and we load those. Data preparation is cheap, so it
# is recomputed here (and walked through visibly in the "Data" chapter).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(mlogit)
  library(knitr)
  library(ggplot2)
})

knitr::opts_chunk$set(fig.align = "center")

output_dir <- "output"

# ── Helpers ───────────────────────────────────────────────────────────────────
# fmt(): fixed-decimal inline numbers.
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# gk(): spell out Greek letters for table/caption text. LaTeX text mode (the PDF
# default font) silently drops Unicode Greek inside tabular cells; body prose
# keeps real symbols via LaTeX math (e.g. $\mu$), which render fine.
gk <- function(x) {
  x <- gsub("χ²", "Chi-sq", x, fixed = TRUE)  # chi-squared
  x <- gsub("β", "beta",  x, fixed = TRUE)
  x <- gsub("α", "alpha", x, fixed = TRUE)
  x <- gsub("σ", "sigma", x, fixed = TRUE)
  x <- gsub("μ", "mu",    x, fixed = TRUE)
  x <- gsub("≈", "~",     x, fixed = TRUE)
  x
}

# ── Pre-estimated model objects ───────────────────────────────────────────────
t3  <- readRDS(file.path(output_dir, "table3_mnl.rds"))
t4h <- readRDS(file.path(output_dir, "table4_hev.rds"))
t4r <- readRDS(file.path(output_dir, "table4_rpl.rds"))
t5  <- readRDS(file.path(output_dir, "table5_wtp.rds"))
iia <- readRDS(file.path(output_dir, "iia_ul.rds"))

# ── Analysis-ready data (mirrors R/01_data_prep.R) ────────────────────────────
alt_labels <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")

price_matrix <- tibble::tribble(
  ~scenario, ~GEN,  ~GT,   ~NAT,  ~CHO,  ~CAB,
  1, 3.38, 7.88, 6.75, 6.75, 6.75,   2, 4.50, 9.00, 6.75, 5.63, 9.00,
  3, 6.75, 9.00, 5.63, 9.00, 6.75,   4, 5.63, 9.00, 9.00, 6.75, 5.63,
  5, 6.75, 6.75, 6.75, 7.88, 7.88,   6, 4.50, 7.88, 7.88, 9.00, 9.00,
  7, 4.50, 6.75, 5.63, 6.75, 7.88,   8, 5.63, 5.63, 5.63, 7.88, 9.00,
  9, 5.63, 7.88, 5.63, 7.88, 9.00,  10, 5.63, 6.75, 7.88, 5.63, 9.00,
 11, 6.75, 7.88, 9.00, 5.63, 7.88,  12, 4.50, 5.63, 5.63, 7.88, 7.88,
 13, 3.38, 9.00, 7.88, 7.88, 9.00,  14, 3.38, 5.63, 5.63, 5.63, 5.63,
 15, 3.38, 6.75, 9.00, 9.00, 9.00,  16, 6.75, 5.63, 7.88, 6.75, 6.75,
 17, 5.63, 5.63, 5.63, 5.63, 5.63
)

.demo_names <- c(
  "subject","gender","age","household","child","education","income","student",
  "occupation","race","gbeef","steak_buy","other_beef","gbeef_freq","steak_freq",
  "color","brand","quality_gd","ex_fat","in_fat","price_imp","safety","juiciness",
  "flavor","tenderness","consistency","doneness","production","quality_gd2",
  "slaughter","food_safe","quality_pref","cooked","chance1","chance2","chance3",
  "chance4","chance5")
.sc_names <- paste0("sc", 1:17)

df_wide <- read_excel("harrison.xls", sheet = "trt 10&11",
                      col_names = c(.demo_names, .sc_names, "blank", "session"),
                      skip = 3, col_types = "text") |>
  mutate(subject = suppressWarnings(as.numeric(subject))) |>
  filter(!is.na(subject), subject >= 1001, subject <= 1199) |>
  filter(subject != 1060) |>
  mutate(treatment = factor(if_else(subject < 1100, "nonhyp", "hyp"),
                            levels = c("nonhyp", "hyp"))) |>
  mutate(across(all_of(c(setdiff(.demo_names, "subject"), .sc_names, "session")),
                as.numeric))

df_long <- df_wide |>
  select(subject, treatment, session, all_of(.sc_names)) |>
  pivot_longer(all_of(.sc_names), names_to = "sc_var", values_to = "choice_code") |>
  mutate(scenario = as.integer(sub("sc", "", sc_var)),
         choice   = factor(choice_code, levels = 1:6, labels = alt_labels)) |>
  select(-sc_var, -choice_code) |>
  left_join(price_matrix, by = "scenario") |>
  arrange(subject, scenario)

df_mlogit <- df_long |>
  filter(scenario <= 16) |>
  crossing(alt = factor(alt_labels, levels = alt_labels)) |>
  mutate(
    chosen    = as.character(choice) == as.character(alt),
    price_alt = case_when(
      alt == "GEN" ~ GEN, alt == "GT"  ~ GT,  alt == "NAT" ~ NAT,
      alt == "CHO" ~ CHO, alt == "CAB" ~ CAB, alt == "NONE" ~ NA_real_),
    chid = paste(subject, scenario, sep = "_")) |>
  select(chid, subject, treatment, session, scenario, alt, chosen, price_alt) |>
  arrange(subject, scenario, alt)
