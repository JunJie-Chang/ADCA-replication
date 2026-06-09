# =============================================================================
# 04_table4_adv.R
# Phase 3: Reproduce Table 4 — RPL and HEV estimates
#
# Models:
#   RPL  — Random Parameters Logit via mlogit
#          ASCs random normal; price params fixed
#          Params: 5 mean_ASC + 5 sd_ASC + 5 price = 15 per model
#
#   HEV  — Heteroscedastic Extreme Value via custom log-likelihood
#          Alternative-specific scale params σ_j (σ_NONE = 1 fixed)
#          Params: 5 ASC + 5 price + 5 log_σ = 15 per model
#
#   MNP  — Multinomial Probit (stub — requires apollo; computationally intensive)
#
# Each model: nonhyp (separate) + hyp (separate) + joint (equal taste + scale μ)
#
# NOTE: Run from RStudio. Unicode path fails in terminal Rscript.
# Requires: output/df_mlogit.rds, output/table3_mnl.rds
# Saves:    output/table4_rpl.rds, output/table4_hev.rds
# =============================================================================

library(dplyr)
library(tidyr)
library(mlogit)
library(knitr)

# ── 0. Paths and data ─────────────────────────────────────────────────────────

project_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
output_dir  <- file.path(project_dir, "output")

df0      <- readRDS(file.path(output_dir, "df_mlogit.rds"))
t3       <- readRDS(file.path(output_dir, "table3_mnl.rds"))

alts <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")
J    <- length(alts)

# ── 1. Prepare variables (same as Table 3) ────────────────────────────────────

df0 <- df0 |>
  mutate(
    alt       = factor(as.character(alt), levels = alts),
    price_alt = replace_na(price_alt, 0),
    p_GEN     = if_else(alt == "GEN", price_alt, 0),
    p_GT      = if_else(alt == "GT",  price_alt, 0),
    p_NAT     = if_else(alt == "NAT", price_alt, 0),
    p_CHO     = if_else(alt == "CHO", price_alt, 0),
    p_CAB     = if_else(alt == "CAB", price_alt, 0)
  )

# Sort by (chid, alt) — required for HEV matrix reshape
df_s <- df0 |>
  mutate(alt_ord = match(as.character(alt), alts)) |>
  arrange(chid, alt_ord) |>
  select(-alt_ord)

# ── 2. dfidx for mlogit (RPL) ─────────────────────────────────────────────────

make_mdata <- function(data) {
  mlogit.data(
    data,
    choice     = "chosen",
    shape      = "long",
    alt.var    = "alt",
    chid.var   = "chid",
    id.var     = "subject",
    alt.levels = alts
  )
}

mdata_nh  <- make_mdata(df0 |> filter(treatment == "nonhyp"))
mdata_hyp <- make_mdata(df0 |> filter(treatment == "hyp"))

f_mnl <- chosen ~ p_GEN + p_GT + p_NAT + p_CHO + p_CAB | 1

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: RPL (Random Parameters Logit)
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 1: RPL — Random Parameters Logit ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# ASC names as they appear in the mlogit object
asc_names   <- t3$asc_names    # e.g. "GEN:(intercept)", ...
price_names <- t3$price_names  # "p_GEN", ...

# rpar spec: ASCs are random normal; price params are fixed
rpar_spec <- setNames(rep("n", 5), asc_names)

cat("Estimating RPL — Nonhypothetical (R=500, panel=TRUE) ...\n")
cat("  (This may take 2–5 minutes)\n")
rpl_nh <- mlogit(
  f_mnl, data = mdata_nh, reflevel = "NONE",
  rpar   = rpar_spec,
  R      = 500,
  halton = NA,
  panel  = TRUE
)
ll_rpl_nh <- as.numeric(logLik(rpl_nh))
cat("  Nonhyp RPL LL:", round(ll_rpl_nh, 3), "\n\n")

cat("Estimating RPL — Hypothetical (R=500, panel=TRUE) ...\n")
cat("  (This may take 2–5 minutes)\n")
rpl_hyp <- mlogit(
  f_mnl, data = mdata_hyp, reflevel = "NONE",
  rpar   = rpar_spec,
  R      = 500,
  halton = NA,
  panel  = TRUE
)
ll_rpl_hyp <- as.numeric(logLik(rpl_hyp))
cat("  Hyp RPL LL:", round(ll_rpl_hyp, 3), "\n\n")

# ── RPL Joint model with scale parameter ──────────────────────────────────────
# Uses simulation-based LL with R=200 Halton draws.
# theta: c(mean_ASC[1:5], price[6:10], log_sd_ASC[11:15], log_mu[16])
# HY group: β_j ~ N(mean_j, exp(log_sd_j)); V = β_j + α_j P
# NH group: V = μ × (β_j + α_j P)

cat("Estimating RPL joint model with scale parameter (R=200) ...\n")
cat("  (This may take 10–30 minutes)\n")

# Halton draws for RPL: R×5 matrix (one column per random ASC)
R_sim <- 200
set.seed(42)
halton_mat <- mlogit:::halton(R_sim, 5)  # R×5 uniform Halton
# Transform to standard normal
halton_norm <- qnorm(halton_mat)          # R×5 N(0,1) draws

ll_rpl_joint_fn <- function(theta, df, R_draws) {
  mu_asc   <- theta[1:5]
  alpha    <- c(theta[6:10], 0)          # fixed price; NONE = 0
  sd_asc   <- exp(theta[11:15])          # positive
  mu_scale <- exp(theta[16])             # scale mu

  alt_idx  <- match(as.character(df$alt), alts)
  asc_base <- c(0, 0, 0, 0, 0, 0)       # placeholder; overwritten per draw
  n_q      <- nrow(R_draws)             # number of draws

  # Group individuals
  subjects <- unique(df$subject)
  n_ind    <- length(subjects)

  log_probs <- numeric(n_ind)

  for (i in seq_len(n_ind)) {
    subj_data <- df |> filter(subject == subjects[i])
    is_nh     <- unique(subj_data$treatment) == "nonhyp"

    # For each draw r: compute panel probability
    panel_probs <- numeric(n_q)
    for (r in seq_len(n_q)) {
      # Draw ASCs for this individual
      beta_r <- mu_asc + sd_asc * R_draws[r, ]   # length 5
      asc_r  <- c(beta_r, 0)                     # NONE ASC = 0

      # Compute V per row
      V_r <- asc_r[alt_idx[match(df$subject, subjects) == i]] +
             alpha[alt_idx[match(df$subject, subjects) == i]] * subj_data$price_alt

      if (is_nh) V_r <- mu_scale * V_r

      # Choice probability per scenario
      cs_prob <- subj_data |>
        mutate(V = V_r) |>
        group_by(chid) |>
        mutate(p = exp(V - max(V)) / sum(exp(V - max(V)))) |>
        ungroup() |>
        filter(chosen) |>
        pull(p)

      panel_probs[r] <- prod(cs_prob)
    }
    log_probs[i] <- log(mean(panel_probs) + 1e-300)
  }

  if (!all(is.finite(log_probs))) return(-1e10)
  sum(log_probs)
}

# NOTE: The loop above is slow due to per-individual R iteration.
# For a faster implementation, use mlogit or apollo with custom setup.
# For now, we skip the joint RPL estimation due to computation time.
# The scale test for RPL uses the same LR framework as MNL.

cat("  [SKIPPED] RPL joint model requires long computation.\n")
cat("  Scale test for RPL: use separate LL sum as benchmark.\n")
cat("  To enable: set estimate_rpl_joint <- TRUE and re-run.\n\n")

estimate_rpl_joint <- FALSE  # set TRUE to enable (very slow)

if (estimate_rpl_joint) {
  # Starting values from separate model averages
  cf_nh_r  <- coef(rpl_nh)
  cf_hyp_r <- coef(rpl_hyp)

  theta0_rpl <- c(
    rowMeans(cbind(cf_nh_r[asc_names],   cf_hyp_r[asc_names])),
    rowMeans(cbind(cf_nh_r[price_names], cf_hyp_r[price_names])),
    log(abs(cf_nh_r[paste0("sd.", asc_names)])),
    0  # log(mu) = 0
  )

  opt_rpl_joint <- optim(
    par     = theta0_rpl,
    fn      = ll_rpl_joint_fn,
    df      = df0,
    R_draws = halton_norm,
    method  = "BFGS",
    control = list(fnscale = -1, maxit = 2000, reltol = 1e-8),
    hessian = TRUE
  )
  cat("  RPL joint LL:", round(opt_rpl_joint$value, 3), "\n")
} else {
  opt_rpl_joint <- NULL
}

# ── RPL summary ───────────────────────────────────────────────────────────────
cat("=== RPL Estimates ===\n")
print(summary(rpl_nh))
cat("\nHypothetical RPL:\n")
print(summary(rpl_hyp))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: HEV (Heteroscedastic Extreme Value)
# ─────────────────────────────────────────────────────────────────────────────
cat("╔══════════════════════════════════════════╗\n")
cat("║  SECTION 2: HEV — Custom Log-Likelihood   ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# ── Quadrature grid ───────────────────────────────────────────────────────────
# Trapezoidal rule on [-8, 45] with 200 points.
# Integrand is negligible outside this range for typical parameter values.

n_q    <- 200
w_lo   <- -8.0
w_hi   <- 45.0
w_grid <- seq(w_lo, w_hi, length.out = n_q)
dw     <- diff(w_grid)[1]
wt_trap <- c(0.5, rep(1.0, n_q - 2), 0.5) * dw
ev_base <- exp(-w_grid)   # exp(-w) part of Gumbel density [n_q]

# ── HEV core log-likelihood ───────────────────────────────────────────────────
# theta: c(asc[1:5], alpha[1:5], log_sigma[11:15])
#   asc[j]      = β_j, ASC for steak j (NONE = 0)
#   alpha[j]    = α_j, own-price coefficient (NONE = 0)
#   log_sigma[j]= log(σ_j) for j = GEN, GT, NAT, CHO, CAB (σ_NONE = 1 fixed)
#
# Data df_s must be sorted by (chid ascending, alt in alts order).
# V_mat[n, j] = utility of alternative j for choice situation n.

hev_ll <- function(theta, df_s, alts, J, n_chid,
                   w_grid, ev_base, wt_trap,
                   mu_scale = 1.0) {
  asc   <- c(theta[1:5], 0.0)
  alpha <- c(theta[6:10], 0.0)
  sigma <- c(exp(theta[11:15]), 1.0)   # σ_NONE = 1 fixed

  alt_idx <- rep(seq_len(J), times = n_chid)   # matches row order after sort

  V_all <- asc[alt_idx] + alpha[alt_idx] * df_s$price_alt

  # Apply group scale if estimating joint model
  if (!is.null(df_s$treatment) && mu_scale != 1.0) {
    V_all[df_s$treatment == "nonhyp"] <- mu_scale * V_all[df_s$treatment == "nonhyp"]
  }

  V_mat  <- matrix(V_all,        nrow = n_chid, ncol = J, byrow = TRUE)
  ch_mat <- matrix(df_s$chosen,  nrow = n_chid, ncol = J, byrow = TRUE)

  n_q_       <- length(w_grid)
  log_probs  <- rep(-1e10, n_chid)

  for (j in seq_len(J)) {
    idx <- which(ch_mat[, j])
    if (length(idx) == 0) next

    ni    <- length(idx)
    V_sub <- V_mat[idx, , drop = FALSE]    # [ni × J]
    sig_j <- sigma[j]
    oth   <- seq_len(J)[-j]
    r_k   <- sigma[oth] / sig_j            # [(J-1)]: ratio σ_k / σ_j

    # d_k[i, k] = V_j[i] - V_k[i]  [ni × (J-1)]
    d_k <- V_sub[, j] - V_sub[, oth, drop = FALSE]

    # exp_const[i, k] = exp(-σ_k × d_k[i,k])  [ni × (J-1)]
    exp_const <- exp(-sweep(d_k, 2, sigma[oth], `*`))

    # exp_r_w[k, q] = exp(-r_k[k] × w_q)  [(J-1) × n_q]
    exp_r_w <- exp(-outer(r_k, w_grid))

    # inner_sum[i, q] = Σ_k exp_const[i,k] × exp_r_w[k,q]  [ni × n_q]
    inner_sum <- exp_const %*% exp_r_w

    # ev_part[i, q] = exp(-w_q) + inner_sum[i, q]
    ev_part <- sweep(inner_sum, 2, ev_base, `+`)

    # log_intgrd[i, q] = -w_q - ev_part[i, q]
    log_intgrd <- -sweep(ev_part, 2, w_grid, `+`)

    # Integrate using log-sum-exp trick per chid
    max_logI  <- apply(log_intgrd, 1, max)
    intgrd_sc <- exp(sweep(log_intgrd, 1, max_logI, `-`))   # [ni × n_q]
    integral  <- as.vector(intgrd_sc %*% wt_trap)            # [ni]

    log_probs[idx] <- max_logI + log(pmax(integral, 1e-300))
  }

  if (!all(is.finite(log_probs))) return(-1e10)
  sum(log_probs)
}

# ── hev_ll_joint: treatment-specific σ (21 parameters) ───────────────────────
# theta[1:5]   → β ASC  (GEN GT NAT CHO CAB; NONE ≡ 0)
# theta[6:10]  → α price (GEN GT NAT CHO CAB; NONE ≡ 0)
# theta[11:16] → log_σ_NH [GEN GT NAT CHO CAB NONE]  — all 6 free
# theta[17:21] → log_σ_H  [GEN GT NAT CHO CAB]       — 5 free; σ_H[NONE] = 1
#
# Identification: σ_{NONE,H} = 1.  σ_{NONE,NH} ≈ 0.49 captures the overall
# difference in error precision between treatments at the NONE alternative.
#
# is_nh: logical vector length n_chid; TRUE = nonhyp choice situation.

hev_ll_joint <- function(theta, df_s, alts, J, n_chid,
                          w_grid, ev_base, wt_trap, is_nh) {
  asc      <- c(theta[1:5],  0.0)
  alpha    <- c(theta[6:10], 0.0)
  sigma_NH <- exp(theta[11:16])
  sigma_H  <- c(exp(theta[17:21]), 1.0)   # σ_H[NONE] = 1 fixed

  alt_idx <- rep(seq_len(J), times = n_chid)
  V_all   <- asc[alt_idx] + alpha[alt_idx] * df_s$price_alt

  V_mat  <- matrix(V_all,       nrow = n_chid, ncol = J, byrow = TRUE)
  ch_mat <- matrix(df_s$chosen, nrow = n_chid, ncol = J, byrow = TRUE)

  log_probs <- rep(-1e10, n_chid)

  for (j in seq_len(J)) {
    idx <- which(ch_mat[, j])
    if (length(idx) == 0) next

    for (grp in list(list(idx = idx[ is_nh[idx]], sigma = sigma_NH),
                     list(idx = idx[!is_nh[idx]], sigma = sigma_H))) {
      if (length(grp$idx) == 0) next

      sigma <- grp$sigma
      V_sub <- V_mat[grp$idx, , drop = FALSE]
      sig_j <- sigma[j]
      oth   <- seq_len(J)[-j]
      r_k   <- sigma[oth] / sig_j

      d_k       <- V_sub[, j] - V_sub[, oth, drop = FALSE]
      exp_const <- exp(-sweep(d_k, 2, sigma[oth], `*`))
      exp_r_w   <- exp(-outer(r_k, w_grid))
      inner_sum <- exp_const %*% exp_r_w
      ev_part   <- sweep(inner_sum, 2, ev_base, `+`)
      log_intgrd <- -sweep(ev_part, 2, w_grid, `+`)

      max_logI  <- apply(log_intgrd, 1, max)
      intgrd_sc <- exp(sweep(log_intgrd, 1, max_logI, `-`))
      integral  <- as.vector(intgrd_sc %*% wt_trap)

      log_probs[grp$idx] <- max_logI + log(pmax(integral, 1e-300))
    }
  }

  if (!all(is.finite(log_probs))) return(-1e10)
  sum(log_probs)
}

# ── Sort subsets and compute n_chid ───────────────────────────────────────────

df_nh  <- df_s |> filter(treatment == "nonhyp")
df_hyp <- df_s |> filter(treatment == "hyp")

n_chid_nh  <- length(unique(df_nh$chid))
n_chid_hyp <- length(unique(df_hyp$chid))
n_chid_all <- length(unique(df_s$chid))

cat("HEV: choice situations —  nonhyp:", n_chid_nh,
    " hyp:", n_chid_hyp, " all:", n_chid_all, "\n\n")

# ── Starting values: from MNL estimates ───────────────────────────────────────
cf_mnl_nh  <- coef(t3$mnl_nh)
cf_mnl_hyp <- coef(t3$mnl_hyp)

hev_theta0_nh <- c(
  cf_mnl_nh[t3$asc_names],
  cf_mnl_nh[t3$price_names],
  rep(0, 5)   # log(σ) = 0 → σ = 1 (reduces to MNL initially)
)
hev_theta0_hyp <- c(
  cf_mnl_hyp[t3$asc_names],
  cf_mnl_hyp[t3$price_names],
  rep(0, 5)
)

# ── Estimate HEV — Nonhypothetical ────────────────────────────────────────────
cat("Estimating HEV — Nonhypothetical ...\n")
cat("  (≈ 2–10 minutes depending on convergence)\n")

opt_hev_nh <- optim(
  par     = hev_theta0_nh,
  fn      = hev_ll,
  df_s    = df_nh,
  alts    = alts,
  J       = J,
  n_chid  = n_chid_nh,
  w_grid  = w_grid,
  ev_base = ev_base,
  wt_trap = wt_trap,
  method  = "BFGS",
  control = list(fnscale = -1, maxit = 5000, reltol = 1e-10),
  hessian = TRUE
)

if (opt_hev_nh$convergence != 0)
  warning("HEV NH convergence code: ", opt_hev_nh$convergence)

ll_hev_nh <- opt_hev_nh$value
cat("  HEV nonhyp LL:", round(ll_hev_nh, 3), "\n")
cat("  σ estimates:  ", round(exp(opt_hev_nh$par[11:15]), 4), "\n\n")

# ── Estimate HEV — Hypothetical ───────────────────────────────────────────────
cat("Estimating HEV — Hypothetical ...\n")

opt_hev_hyp <- optim(
  par     = hev_theta0_hyp,
  fn      = hev_ll,
  df_s    = df_hyp,
  alts    = alts,
  J       = J,
  n_chid  = n_chid_hyp,
  w_grid  = w_grid,
  ev_base = ev_base,
  wt_trap = wt_trap,
  method  = "BFGS",
  control = list(fnscale = -1, maxit = 5000, reltol = 1e-10),
  hessian = TRUE
)

if (opt_hev_hyp$convergence != 0)
  warning("HEV HY convergence code: ", opt_hev_hyp$convergence)

ll_hev_hyp <- opt_hev_hyp$value
cat("  HEV hyp LL:", round(ll_hev_hyp, 3), "\n")
cat("  σ estimates:", round(exp(opt_hev_hyp$par[11:15]), 4), "\n\n")

# ── Estimate HEV — Joint model (21 parameters) ────────────────────────────────
# Taste params shared; scale params differ by alternative AND treatment.
# is_nh: one TRUE/FALSE per choice situation (sampled every J rows).

is_nh_all <- df_s$treatment[seq(1, nrow(df_s), by = J)] == "nonhyp"

# Starting values: taste = average of separate models
# σ_NH from nonhyp model (5 free) + 0 for NONE (log(1) = 0)
# σ_H  from hyp model   (5 free); NONE is fixed at 1 during estimation
theta0_hev_joint <- c(
  rowMeans(cbind(opt_hev_nh$par[1:10], opt_hev_hyp$par[1:10])),  # taste (10)
  c(opt_hev_nh$par[11:15], 0),   # log_σ_NH: GEN GT NAT CHO CAB NONE
  opt_hev_hyp$par[11:15]          # log_σ_H:  GEN GT NAT CHO CAB
)                                  # total: 10 + 6 + 5 = 21 params

cat("Estimating HEV — Joint model (21 params: 10 taste + 6 σ_NH + 5 σ_H) ...\n")
cat("  (≈ 15–45 minutes)\n")

opt_hev_joint <- optim(
  par     = theta0_hev_joint,
  fn      = hev_ll_joint,
  df_s    = df_s,
  alts    = alts,
  J       = J,
  n_chid  = n_chid_all,
  w_grid  = w_grid,
  ev_base = ev_base,
  wt_trap = wt_trap,
  is_nh   = is_nh_all,
  method  = "BFGS",
  control = list(fnscale = -1, maxit = 5000, reltol = 1e-10),
  hessian = TRUE
)

if (opt_hev_joint$convergence != 0)
  warning("HEV joint convergence code: ", opt_hev_joint$convergence)

ll_hev_joint <- opt_hev_joint$value
sigma_NH_hat <- exp(opt_hev_joint$par[11:16])          # 6 values incl. NONE
sigma_H_hat  <- c(exp(opt_hev_joint$par[17:21]), 1.0)  # 6 values; NONE = 1

cat(sprintf("  HEV joint LL: %.3f  (paper: −2031)\n", ll_hev_joint))
cat(sprintf("  σ_NH: %s  None=%.2f\n",
    paste(round(sigma_NH_hat[1:5], 2), collapse = ", "), sigma_NH_hat[6]))
cat(sprintf("  σ_H:  %s  None=1.00 (fixed)\n\n",
    paste(round(sigma_H_hat[1:5], 2), collapse = ", ")))

# ── HEV equality test ─────────────────────────────────────────────────────────
# Separate: 15 params × 2 = 30  |  Joint: 21 params  |  df = 30 − 21 = 9
lr_hev <- -2 * (ll_hev_joint - ll_hev_nh - ll_hev_hyp)
df_hev <- (15 + 15) - 21
p_hev  <- pchisq(lr_hev, df = df_hev, lower.tail = FALSE)
cat(sprintf("HEV equality test: χ²(%d) = %.3f  p = %.4f  (paper: χ² = 294.8)\n\n",
            df_hev, lr_hev, p_hev))

# ── SE from Hessian ───────────────────────────────────────────────────────────
vcov_hev_nh <- tryCatch(solve(-opt_hev_nh$hessian),
                         error = function(e) matrix(NA, 15, 15))
vcov_hev_hyp <- tryCatch(solve(-opt_hev_hyp$hessian),
                          error = function(e) matrix(NA, 15, 15))
vcov_hev_joint <- tryCatch(solve(-opt_hev_joint$hessian),
                            error = function(e) matrix(NA, 21, 21))

se_hev_nh    <- sqrt(pmax(diag(vcov_hev_nh),    0))
se_hev_hyp   <- sqrt(pmax(diag(vcov_hev_hyp),   0))
se_hev_joint <- sqrt(pmax(diag(vcov_hev_joint),  0))

# SE for σ via delta method: SE(σ_j) = σ_j × SE(log_σ_j)
sigma_nh_sep     <- c(exp(opt_hev_nh$par[11:15]),  1.0)  # 6 vals, NONE fixed
sigma_hyp_sep    <- c(exp(opt_hev_hyp$par[11:15]), 1.0)
se_sigma_nh_sep  <- c(sigma_nh_sep[1:5]  * se_hev_nh[11:15],  0)
se_sigma_hyp_sep <- c(sigma_hyp_sep[1:5] * se_hev_hyp[11:15], 0)
se_sigma_NH      <- sigma_NH_hat * se_hev_joint[11:16]
se_sigma_H       <- c(sigma_H_hat[1:5] * se_hev_joint[17:21], 0)  # NONE SE = 0

# ── Format HEV Table ──────────────────────────────────────────────────────────
# Layout mirrors Table 4 in the paper:
#   Rows 1–10 : taste params (ASC + price)  — all three columns populated
#   Rows 11–16: σ_NH params                 — NH col + Joint col
#   Rows 17–22: σ_H  params                 — HY col + Joint col

param_labels_taste <- c(
  paste0("beta ", c("GEN","GT","NAT","CHO","CAB"), " (ASC)"),
  paste0("alpha ", c("GEN","GT","NAT","CHO","CAB"), " (price)")
)
param_labels_scale <- c(
  paste0("sigma_NH ", c("GEN","GT","NAT","CHO","CAB","None")),
  paste0("sigma_H  ", c("GEN","GT","NAT","CHO","CAB","None"))
)

taste_tbl <- data.frame(
  Parameter = param_labels_taste,
  NH_Est = round(opt_hev_nh$par[1:10],    4),
  NH_SE  = round(se_hev_nh[1:10],          4),
  HY_Est = round(opt_hev_hyp$par[1:10],   4),
  HY_SE  = round(se_hev_hyp[1:10],         4),
  JT_Est = round(opt_hev_joint$par[1:10], 4),
  JT_SE  = round(se_hev_joint[1:10],       4),
  row.names = NULL
)

scale_tbl <- data.frame(
  Parameter = param_labels_scale,
  NH_Est = c(round(sigma_nh_sep,     4), rep(NA_real_, 6)),
  NH_SE  = c(round(se_sigma_nh_sep,  4), rep(NA_real_, 6)),
  HY_Est = c(rep(NA_real_, 6), round(sigma_hyp_sep,     4)),
  HY_SE  = c(rep(NA_real_, 6), round(se_sigma_hyp_sep,  4)),
  JT_Est = c(round(sigma_NH_hat, 4), round(sigma_H_hat, 4)),
  JT_SE  = c(round(se_sigma_NH,  4), round(se_sigma_H,  4)),
  row.names = NULL
)

table4_hev <- rbind(taste_tbl, scale_tbl)

cat("=== TABLE 4 (HEV) — Taste Parameters ===\n\n")
print(kable(taste_tbl, format = "simple",
            col.names = c("Parameter","NH Est","(SE)","HY Est","(SE)","JT Est","(SE)")))
cat("\n=== TABLE 4 (HEV) — Scale Parameters ===\n\n")
print(kable(scale_tbl, format = "simple",
            col.names = c("Parameter","NH Est","(SE)","HY Est","(SE)","JT Est","(SE)")))
cat(sprintf("\n  LL: NH=%.3f  HY=%.3f  Joint=%.3f\n",
            ll_hev_nh, ll_hev_hyp, ll_hev_joint))
cat(sprintf("  Equality test: chi2(%d) = %.3f  p = %.4f  (paper: 294.8)\n\n",
            df_hev, lr_hev, p_hev))

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: MNP stub
# ─────────────────────────────────────────────────────────────────────────────
cat("╔═════════════════════════════════════╗\n")
cat("║  SECTION 3: MNP — Apollo stub only  ║\n")
cat("╚═════════════════════════════════════╝\n\n")
cat("MNP is computationally intensive (GHK simulator, 5-dim MVN integral).\n")
cat("Recommended implementation: use the `apollo` package.\n")
cat("See apollo documentation: http://www.apollochoicemodelling.com/\n\n")
cat("Key apollo settings for MNP (diagonal covariance, 2 SDs fixed to 1):\n")
cat("  apollo_probabilities → apollo_mnp(mnp_settings, functionality)\n")
cat("  normalCovSd fixed at 1 for GEN and NONE; free for GT, NAT, CHO, CAB\n")
cat("[SKIPPED]\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Save results
# ─────────────────────────────────────────────────────────────────────────────

results_rpl <- list(
  rpl_nh      = rpl_nh,
  rpl_hyp     = rpl_hyp,
  rpl_joint   = opt_rpl_joint,
  ll_nh       = ll_rpl_nh,
  ll_hyp      = ll_rpl_hyp,
  asc_names   = asc_names,
  price_names = price_names
)

results_hev <- list(
  opt_nh       = opt_hev_nh,
  opt_hyp      = opt_hev_hyp,
  opt_joint    = opt_hev_joint,
  vcov_nh      = vcov_hev_nh,
  vcov_hyp     = vcov_hev_hyp,
  vcov_joint   = vcov_hev_joint,
  ll_nh        = ll_hev_nh,
  ll_hyp       = ll_hev_hyp,
  ll_joint     = ll_hev_joint,
  sigma_NH_hat = sigma_NH_hat,   # exp(par[11:16]): GEN GT NAT CHO CAB NONE
  sigma_H_hat  = sigma_H_hat,    # c(exp(par[17:21]), 1.0): GEN GT NAT CHO CAB NONE
  lr_stat      = lr_hev,
  df_test      = df_hev,
  p_lr         = p_hev,
  table4_hev   = table4_hev,
  taste_tbl    = taste_tbl,
  scale_tbl    = scale_tbl,
  # Raw par vectors for Table 5 WTP (separate model estimates, 15 params each)
  par_nh       = opt_hev_nh$par[1:15],
  par_hyp      = opt_hev_hyp$par[1:15]
)

saveRDS(results_rpl, file.path(output_dir, "table4_rpl.rds"))
saveRDS(results_hev, file.path(output_dir, "table4_hev.rds"))

write.csv(table4_hev,
          file.path(output_dir, "table4_hev.csv"), row.names = FALSE)

cat("Saved: output/table4_rpl.rds\n")
cat("Saved: output/table4_hev.rds\n")
cat("Saved: output/table4_hev.csv\n")
