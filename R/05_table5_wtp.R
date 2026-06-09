# =============================================================================
# 05_table5_wtp.R
# Phase 4: Reproduce Table 5 (WTP + Poe test) and Figure 2 (HEV market shares)
#
# WTP formula: WTP_j = -β_j / α_j  (willingness to pay for steak j vs. NONE)
#   MNL: β_j = ASC, α_j = own-price coefficient
#   RPL: β_j = mean ASC (μ_j), α_j = fixed price coefficient
#   HEV: β_j = ASC, α_j = own-price coefficient (σ_j does not enter WTP)
#
# Poe combinatorial test (one-tailed):
#   H0: WTP_NH ≤ WTP_HY  (hypothetical ≥ nonhypothetical)
#   H1: WTP_NH > WTP_HY  (nonhypothetical > hypothetical)
#   Draw B=1000 bootstrap vectors from MVN(coef, vcov) for each group.
#   p-value = Pr(WTP_NH[b] - WTP_HY[b'] < 0) over all 1000×1000 pairs.
#   Reference: Poe, Giraud & Loomis (2005), American Journal of Agricultural Economics.
#
# Figure 2: Predicted market shares from HEV nonhyp and hyp models
#   For each of the 17 scenarios, compute P(choose alt j | HEV params).
#   Uses the same quadrature as 04_table4_adv.R.
#
# NOTE: Run from RStudio. Unicode path fails in terminal Rscript.
# Requires: output/table3_mnl.rds, output/table4_rpl.rds, output/table4_hev.rds,
#           output/df_long.rds, output/price_matrix.rds
# Saves:    output/table5_wtp.rds, output/table5_wtp.csv, output/figure2.pdf
# =============================================================================

library(dplyr)
library(tidyr)
library(MASS)        # mvrnorm for Poe bootstrap
library(Matrix)      # nearPD — fix near-singular vcov matrices
library(knitr)
library(ggplot2)

# ── 0. Paths and data ─────────────────────────────────────────────────────────

project_dir  <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
output_dir   <- file.path(project_dir, "output")

t3   <- readRDS(file.path(output_dir, "table3_mnl.rds"))
t4r  <- readRDS(file.path(output_dir, "table4_rpl.rds"))
t4h  <- readRDS(file.path(output_dir, "table4_hev.rds"))
df_long      <- readRDS(file.path(output_dir, "df_long.rds"))
price_matrix <- readRDS(file.path(output_dir, "price_matrix.rds"))

alts      <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")
steak_alts <- alts[1:5]   # exclude NONE from WTP

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: WTP Computation
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 1: WTP Estimates                 ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Helper: compute WTP from coef vector
# theta_asc[1:5]   = ASCs for GEN, GT, NAT, CHO, CAB (in that order)
# theta_price[6:10] = price params for GEN, ..., CAB
# WTP_j = -asc_j / alpha_j
compute_wtp <- function(theta) {
  asc   <- theta[1:5]
  alpha <- theta[6:10]
  -asc / alpha
}

# ── MNL WTP ───────────────────────────────────────────────────────────────────

asc_names   <- t3$asc_names    # "GEN:(intercept)", ...
price_names <- t3$price_names  # "p_GEN", ...

cf_mnl_nh  <- coef(t3$mnl_nh)
cf_mnl_hyp <- coef(t3$mnl_hyp)

# Reorder: [asc_GEN, asc_GT, ..., alpha_GEN, ..., alpha_CAB]
theta_mnl_nh  <- c(cf_mnl_nh[asc_names],   cf_mnl_nh[price_names])
theta_mnl_hyp <- c(cf_mnl_hyp[asc_names],  cf_mnl_hyp[price_names])

wtp_mnl_nh  <- compute_wtp(theta_mnl_nh)
wtp_mnl_hyp <- compute_wtp(theta_mnl_hyp)

names(wtp_mnl_nh) <- names(wtp_mnl_hyp) <- steak_alts

cat("MNL WTP (nonhyp):", round(wtp_mnl_nh, 3), "\n")
cat("MNL WTP (hyp):   ", round(wtp_mnl_hyp, 3), "\n\n")

# ── RPL WTP ───────────────────────────────────────────────────────────────────
# Mean WTP = -mean_ASC_j / alpha_j

cf_rpl_nh  <- coef(t4r$rpl_nh)
cf_rpl_hyp <- coef(t4r$rpl_hyp)

# mlogit RPL coef: first 5 = price params (p_GEN, ...), next 5 = mean ASCs,
# then 5 sd.ASCs.  Extract carefully by name.
theta_rpl_nh  <- c(cf_rpl_nh[asc_names],   cf_rpl_nh[price_names])
theta_rpl_hyp <- c(cf_rpl_hyp[asc_names],  cf_rpl_hyp[price_names])

wtp_rpl_nh  <- compute_wtp(theta_rpl_nh)
wtp_rpl_hyp <- compute_wtp(theta_rpl_hyp)
names(wtp_rpl_nh) <- names(wtp_rpl_hyp) <- steak_alts

cat("RPL WTP (nonhyp):", round(wtp_rpl_nh, 3), "\n")
cat("RPL WTP (hyp):   ", round(wtp_rpl_hyp, 3), "\n\n")

# ── HEV WTP ───────────────────────────────────────────────────────────────────
# par_nh[1:5] = ASCs, par_nh[6:10] = alpha, par_nh[11:15] = log_sigma

theta_hev_nh  <- t4h$par_nh[1:10]   # [ASC, alpha]
theta_hev_hyp <- t4h$par_hyp[1:10]

wtp_hev_nh  <- compute_wtp(theta_hev_nh)
wtp_hev_hyp <- compute_wtp(theta_hev_hyp)
names(wtp_hev_nh) <- names(wtp_hev_hyp) <- steak_alts

cat("HEV WTP (nonhyp):", round(wtp_hev_nh, 3), "\n")
cat("HEV WTP (hyp):   ", round(wtp_hev_hyp, 3), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Poe Combinatorial Test
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 2: Poe Combinatorial Test        ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Poe test for H0: WTP_NH = WTP_HY
# p-value = Pr(WTP_NH[b] - WTP_HY[b'] < 0) over B×B bootstrap pairs
# A p-value > 0.05 means we cannot reject equality of WTP.

make_pd <- function(S) {
  # If S is already positive definite return as-is; otherwise use nearPD.
  tryCatch(
    { chol(S); S },
    error = function(e) as.matrix(Matrix::nearPD(S, ensureSymmetry = TRUE)$mat)
  )
}

poe_test <- function(vcov_nh, vcov_hyp,
                     theta_nh, theta_hyp,
                     B = 1000, seed = 123) {
  set.seed(seed)
  # Ensure positive definiteness before drawing (HEV vcov can be near-singular
  # when a scale parameter drifts to a very large value during optimisation).
  vcov_nh  <- make_pd(vcov_nh)
  vcov_hyp <- make_pd(vcov_hyp)
  # Draw B coefficient vectors from MVN
  draws_nh  <- MASS::mvrnorm(B, mu = theta_nh,  Sigma = vcov_nh)
  draws_hyp <- MASS::mvrnorm(B, mu = theta_hyp, Sigma = vcov_hyp)

  # Compute WTP for each draw: matrix B × 5
  wtp_nh_draws  <- t(apply(draws_nh,  1, compute_wtp))
  wtp_hyp_draws <- t(apply(draws_hyp, 1, compute_wtp))

  # For each steak, compute Poe p-value:
  # p_j = proportion of (b, b') pairs where WTP_NH[b,j] - WTP_HY[b',j] < 0
  p_vals <- sapply(seq_len(5), function(j) {
    # Efficient: for each b, count fraction of b' where NH[b,j] < HY[b',j]
    mean(sapply(wtp_nh_draws[, j], function(w_nh) mean(w_nh < wtp_hyp_draws[, j])))
  })
  names(p_vals) <- steak_alts
  p_vals
}

# ── MNL Poe test ──────────────────────────────────────────────────────────────
cat("Computing Poe test — MNL (B=1000, 1M comparisons per steak) ...\n")

vcov_mnl_nh  <- vcov(t3$mnl_nh)[c(asc_names, price_names),
                                  c(asc_names, price_names)]
vcov_mnl_hyp <- vcov(t3$mnl_hyp)[c(asc_names, price_names),
                                   c(asc_names, price_names)]

poe_mnl <- poe_test(vcov_mnl_nh, vcov_mnl_hyp,
                    theta_mnl_nh, theta_mnl_hyp)
cat("  MNL Poe p-values:", round(poe_mnl, 4), "\n\n")

# ── RPL Poe test ──────────────────────────────────────────────────────────────
cat("Computing Poe test — RPL (B=1000) ...\n")

# Reorder vcov to match theta ordering: [asc, price]
vcov_rpl_nh  <- vcov(t4r$rpl_nh)
vcov_rpl_hyp <- vcov(t4r$rpl_hyp)

# Subset to ASC and price params only (exclude sd.ASC for WTP bootstrap)
ord_names <- c(asc_names, price_names)
vcov_rpl_nh_sub  <- vcov_rpl_nh[ord_names, ord_names]
vcov_rpl_hyp_sub <- vcov_rpl_hyp[ord_names, ord_names]

poe_rpl <- poe_test(vcov_rpl_nh_sub, vcov_rpl_hyp_sub,
                    theta_rpl_nh, theta_rpl_hyp)
cat("  RPL Poe p-values:", round(poe_rpl, 4), "\n\n")

# ── HEV Poe test ──────────────────────────────────────────────────────────────
cat("Computing Poe test — HEV (B=1000) ...\n")

# HEV vcov is 15×15; use only first 10 rows/cols [ASC, price]
vcov_hev_nh_sub  <- t4h$vcov_nh[1:10,  1:10]
vcov_hev_hyp_sub <- t4h$vcov_hyp[1:10, 1:10]

poe_hev <- poe_test(vcov_hev_nh_sub, vcov_hev_hyp_sub,
                    theta_hev_nh, theta_hev_hyp)
cat("  HEV Poe p-values:", round(poe_hev, 4), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Format Table 5
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 3: Table 5                       ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

steak_labels <- c("Generic", "Guaranteed Tender", "Natural",
                  "USDA Choice", "Certified Angus Beef")

table5 <- data.frame(
  Steak        = steak_labels,
  # MNL
  MNL_NH       = round(wtp_mnl_nh,  2),
  MNL_HY       = round(wtp_mnl_hyp, 2),
  MNL_Poe      = round(poe_mnl,     4),
  # HEV
  HEV_NH       = round(wtp_hev_nh,  2),
  HEV_HY       = round(wtp_hev_hyp, 2),
  HEV_Poe      = round(poe_hev,     4),
  # RPL
  RPL_NH       = round(wtp_rpl_nh,  2),
  RPL_HY       = round(wtp_rpl_hyp, 2),
  RPL_Poe      = round(poe_rpl,     4),
  row.names    = NULL
)

cat("=== TABLE 5: Willingness-to-Pay Estimates ===\n\n")
print(kable(
  table5,
  format    = "simple",
  col.names = c("Steak",
                "NH($)", "HY($)", "Poe p",
                "NH($)", "HY($)", "Poe p",
                "NH($)", "HY($)", "Poe p"),
  caption   = "Table 5. WTP (vs. NONE) by model and treatment"
))
cat("\n  NH = Nonhypothetical  HY = Hypothetical\n")
cat("  Poe p: one-tailed, H0: WTP_NH ≤ WTP_HY, B=1000 draws\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Figure 2 — Predicted Market Shares (HEV)
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 4: Figure 2 — HEV Market Shares  ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Use same quadrature grid as Table 4
n_q    <- 200
w_lo   <- -8.0; w_hi <- 45.0
w_grid <- seq(w_lo, w_hi, length.out = n_q)
dw     <- diff(w_grid)[1]
wt_trap <- c(0.5, rep(1.0, n_q - 2), 0.5) * dw
ev_base <- exp(-w_grid)

# ── HEV choice probability for a single choice situation ──────────────────────
# V_vec: length-J utility vector
# sigma: length-J scale vector
# Returns: length-J probability vector

hev_probs_single <- function(V_vec, sigma, w_grid, ev_base, wt_trap) {
  J   <- length(V_vec)
  n_q <- length(w_grid)
  probs <- numeric(J)

  for (j in seq_len(J)) {
    oth   <- seq_len(J)[-j]
    sig_j <- sigma[j]
    r_k   <- sigma[oth] / sig_j
    d_k   <- V_vec[j] - V_vec[oth]

    exp_const <- exp(-sigma[oth] * d_k)                        # [(J-1)]
    exp_r_w   <- exp(outer(-r_k, w_grid))                      # [(J-1) × n_q]
    inner_sum <- as.vector(exp_const %*% exp_r_w)               # [n_q]
    ev_part   <- ev_base + inner_sum                            # [n_q]
    log_intgrd <- -w_grid - ev_part

    max_logI  <- max(log_intgrd)
    integral  <- sum(exp(log_intgrd - max_logI) * wt_trap)
    probs[j]  <- max_logI + log(max(integral, 1e-300))
  }

  # Normalise (log-space)
  probs <- exp(probs - max(probs))
  probs / sum(probs)
}

# ── Predicted market share for given HEV parameters across all 17 scenarios ──

predict_market_shares_hev <- function(par15, price_matrix, alts) {
  asc   <- c(par15[1:5], 0)
  alpha <- c(par15[6:10], 0)
  sigma <- c(exp(par15[11:15]), 1.0)
  J     <- length(alts)

  shares <- matrix(NA_real_, nrow = 17, ncol = J,
                   dimnames = list(NULL, alts))

  for (s in 1:17) {
    prices <- c(as.numeric(price_matrix[s, c("GEN","GT","NAT","CHO","CAB")]), 0)
    V_vec  <- asc + alpha * prices
    shares[s, ] <- hev_probs_single(V_vec, sigma, w_grid, ev_base, wt_trap)
  }

  as.data.frame(shares) |>
    mutate(scenario = 1:17) |>
    pivot_longer(-scenario, names_to = "alternative", values_to = "share")
}

cat("Computing HEV predicted market shares (17 scenarios × 6 alts × 2 groups) ...\n")

shares_nh  <- predict_market_shares_hev(t4h$par_nh,  price_matrix, alts) |>
  mutate(treatment = "Nonhypothetical")
shares_hyp <- predict_market_shares_hev(t4h$par_hyp, price_matrix, alts) |>
  mutate(treatment = "Hypothetical")

shares_all <- bind_rows(shares_nh, shares_hyp) |>
  mutate(
    alternative = factor(alternative, levels = alts),
    treatment   = factor(treatment,   levels = c("Nonhypothetical","Hypothetical"))
  )

# ── Observed market shares (for comparison) ────────────────────────────────────
obs_shares <- df_long |>
  group_by(treatment, scenario, choice) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(treatment, scenario) |>
  mutate(share = n / sum(n)) |>
  ungroup() |>
  dplyr::select(scenario, alternative = choice, share, treatment) |>
  mutate(
    treatment   = if_else(treatment == "nonhyp", "Nonhypothetical", "Hypothetical"),
    alternative = as.character(alternative)
  )

# ── Plot Figure 2 ─────────────────────────────────────────────────────────────

alt_colours <- c(
  GEN  = "#1f77b4",
  GT   = "#ff7f0e",
  NAT  = "#2ca02c",
  CHO  = "#d62728",
  CAB  = "#9467bd",
  NONE = "#7f7f7f"
)

fig2 <- ggplot(
  shares_all |> filter(scenario <= 16),   # scenarios 1–16 only (used in estimation)
  aes(x = scenario, y = share, colour = alternative, linetype = treatment)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~treatment, ncol = 1) +
  scale_colour_manual(
    values = alt_colours,
    name   = "Alternative",
    labels = c(GEN="Generic", GT="Guar. Tender",
               NAT="Natural", CHO="USDA Choice",
               CAB="CAB", NONE="None")
  ) +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title    = "Figure 2. Predicted Market Shares — HEV Model",
    subtitle = "Scenarios 1–16 (scenario 17 excluded from estimation)",
    x        = "Scenario",
    y        = "Predicted choice probability",
    caption  = "Replication of Lusk & Schroeder (2004) Figure 2"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "right",
    strip.background  = element_rect(fill = "grey90"),
    panel.grid.minor  = element_blank()
  )

fig2_path <- file.path(output_dir, "figure2_hev_market_shares.pdf")
ggsave(fig2_path, fig2, width = 9, height = 7)
cat("Saved: output/figure2_hev_market_shares.pdf\n")

# Also print summary of shares at key scenarios
cat("\nPredicted market shares — Scenario 1 (NH):\n")
print(
  shares_all |>
    filter(scenario == 1, treatment == "Nonhypothetical") |>
    dplyr::select(alternative, share) |>
    mutate(share = round(share, 4))
)
cat("\nPredicted market shares — Scenario 1 (HY):\n")
print(
  shares_all |>
    filter(scenario == 1, treatment == "Hypothetical") |>
    dplyr::select(alternative, share) |>
    mutate(share = round(share, 4))
)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Save
# ─────────────────────────────────────────────────────────────────────────────

results <- list(
  # WTP estimates
  wtp_mnl_nh  = wtp_mnl_nh,
  wtp_mnl_hyp = wtp_mnl_hyp,
  wtp_rpl_nh  = wtp_rpl_nh,
  wtp_rpl_hyp = wtp_rpl_hyp,
  wtp_hev_nh  = wtp_hev_nh,
  wtp_hev_hyp = wtp_hev_hyp,
  # Poe p-values
  poe_mnl     = poe_mnl,
  poe_rpl     = poe_rpl,
  poe_hev     = poe_hev,
  # Formatted table
  table5      = table5,
  # Market shares
  shares_all  = shares_all
)

saveRDS(results, file.path(output_dir, "table5_wtp.rds"))
write.csv(table5, file.path(output_dir, "table5_wtp.csv"), row.names = FALSE)

cat("\nSaved: output/table5_wtp.rds\n")
cat("Saved: output/table5_wtp.csv\n")
cat("\n=== All outputs complete ===\n")
cat("Tables: table3_mnl.csv, table4_hev.csv, table5_wtp.csv\n")
cat("Figure: figure2_hev_market_shares.pdf\n")
