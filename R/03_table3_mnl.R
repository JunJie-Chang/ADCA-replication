# =============================================================================
# 03_table3_mnl.R
# Phase 2: Reproduce Table 3 — Multinomial Logit (MNL) Estimates
#
# Models estimated:
#   (A) Nonhypothetical MNL (separate, 1072 choice obs = 67 × 16)
#   (B) Hypothetical MNL    (separate,  592 choice obs = 37 × 16)
#   (C) Pooled MNL          (equal params, no scale — LR test benchmark)
#   (D) Joint MNL           (equal taste params + free scale μ)
#
# Utility: V_ij = β_j + α_j × P_ij  (NONE = base: β_NONE = α_NONE = 0)
#   - 5 alternative-specific constants (ASCs): β_GEN, β_GT, β_NAT, β_CHO, β_CAB
#   - 5 own-price parameters: α_GEN, α_GT, α_NAT, α_CHO, α_CAB
#   → 10 taste parameters per model
#
# Joint model (D):
#   V_ij^NH = μ × (β_j + α_j × P_ij)     (nonhyp group scaled by μ)
#   V_ij^H  =      β_j + α_j × P_ij      (hyp group, scale = 1)
#   μ = scale parameter (nonhyp precision relative to hyp)
#   → 11 parameters total (10 taste + 1 scale)
#
# Scale test (LR): -2 × (LL_joint - LL_nonhyp - LL_hyp)
#   df = (10 + 10) - 11 = 9
#
# NOTE: Run from RStudio. The project folder contains Unicode characters
# that Rscript in terminal cannot encode properly.
#
# Requires: output/df_mlogit.rds (from 01_data_prep.R)
# Saves:    output/table3_mnl.rds, output/table3_mnl.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(mlogit)
library(knitr)

# ── 0. Paths and data ─────────────────────────────────────────────────────────

project_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
output_dir  <- file.path(project_dir, "output")

df0 <- readRDS(file.path(output_dir, "df_mlogit.rds"))

cat("Loaded df_mlogit:\n")
cat("  nonhyp choice situations:", sum(df0$treatment == "nonhyp" & df0$alt == "GEN"),
    "(paper: 1072)\n")
cat("  hyp    choice situations:", sum(df0$treatment == "hyp"    & df0$alt == "GEN"),
    "(paper:  592)\n\n")

# ── 1. Prepare alternative-specific price variables ───────────────────────────
# NONE has no price → set to 0 (α_NONE = 0 by construction).
# Create 5 alternative-specific price dummies p_j = price × I(alt == j),
# which appear as generic variables in mlogit but are effectively alt-specific.

alts <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")

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

# ── 2. mlogit.data objects ────────────────────────────────────────────────────
# Using mlogit.data() instead of dfidx() to avoid a node stack overflow bug
# in dfidx 0.2.0 + vctrs 0.6.5.
# NONE is last in alt.levels → reference alternative (ASC forced to 0).

make_mdata <- function(data) {
  mlogit.data(
    data,
    choice   = "chosen",
    shape    = "long",
    alt.var  = "alt",
    chid.var = "chid",
    id.var   = "subject",
    alt.levels = alts
  )
}

mdata_nh  <- make_mdata(df0 |> filter(treatment == "nonhyp"))
mdata_hyp <- make_mdata(df0 |> filter(treatment == "hyp"))
mdata_all <- make_mdata(df0)

# ── 3. MNL formula ────────────────────────────────────────────────────────────
# Part 1 (generic coef):  p_GEN + p_GT + p_NAT + p_CHO + p_CAB
#   → 5 own-price parameters (each active for only one alternative)
# Part 2 (alt-specific):  | 1  → 5 ASCs (NONE excluded as reference)

f_mnl <- chosen ~ p_GEN + p_GT + p_NAT + p_CHO + p_CAB | 1

# ── 4. Estimate separate MNL models (A and B) ─────────────────────────────────

cat("Estimating nonhypothetical MNL...\n")
mnl_nh  <- mlogit(f_mnl, data = mdata_nh,  reflevel = "NONE")

cat("Estimating hypothetical MNL...\n")
mnl_hyp <- mlogit(f_mnl, data = mdata_hyp, reflevel = "NONE")

ll_nh  <- as.numeric(logLik(mnl_nh))
ll_hyp <- as.numeric(logLik(mnl_hyp))
ll_unres <- ll_nh + ll_hyp

cat("  Nonhyp LL:", round(ll_nh,  3), "\n")
cat("  Hyp    LL:", round(ll_hyp, 3), "\n")
cat("  Sum LL:   ", round(ll_unres, 3), "(unrestricted)\n\n")

# Identify coefficient names (mlogit names ASCs as "alt:(intercept)")
cf_nh  <- coef(mnl_nh)
cf_hyp <- coef(mnl_hyp)
cat("Coefficient names (nonhyp):\n")
print(names(cf_nh))
cat("\n")

# Detect ASC and price coefficient names robustly
price_names <- c("p_GEN", "p_GT", "p_NAT", "p_CHO", "p_CAB")
asc_names   <- setdiff(names(cf_nh), price_names)

cat("ASC names:  ", paste(asc_names,   collapse = ", "), "\n")
cat("Price names:", paste(price_names, collapse = ", "), "\n\n")

# ── 5. Pooled MNL — equal params, no scale (benchmark for LR test) ────────────

cat("Estimating pooled MNL (restricted: equal params)...\n")
mnl_pool <- mlogit(f_mnl, data = mdata_all, reflevel = "NONE")
ll_pool  <- as.numeric(logLik(mnl_pool))
cat("  Pooled LL:", round(ll_pool, 3), "\n\n")

# ── 6. Joint model with scale parameter ───────────────────────────────────────
# Custom log-likelihood:
#   θ = (β_GEN, β_GT, β_NAT, β_CHO, β_CAB,       [indices 1:5]  ASCs
#         α_GEN, α_GT, α_NAT, α_CHO, α_CAB,       [indices 6:10] price
#         log_μ)                                   [index 11]     scale
#
# V_nonhyp = μ × (β_j + α_j × P)
# V_hyp    =       β_j + α_j × P
# V_NONE   = 0 (always)

ll_joint_fn <- function(theta, df) {
  beta  <- c(theta[1:5], 0)   # ASCs; NONE = 0
  alpha <- c(theta[6:10], 0)  # price params; NONE = 0
  mu    <- exp(theta[11])     # scale (constrained > 0)

  alt_idx <- match(as.character(df$alt), alts)

  V <- beta[alt_idx] + alpha[alt_idx] * df$price_alt

  # Scale nonhypothetical utilities
  V[df$treatment == "nonhyp"] <- mu * V[df$treatment == "nonhyp"]

  # Log-sum-exp for numerical stability; sum log-probs of chosen options
  df_tmp <- df |>
    mutate(V = V) |>
    group_by(chid) |>
    mutate(
      V_s  = V - max(V),
      expV = exp(V_s),
      prob = expV / sum(expV)
    ) |>
    ungroup()

  ll <- sum(log(pmax(df_tmp$prob[df_tmp$chosen], 1e-300)))
  if (!is.finite(ll)) return(-1e10)
  ll
}

# Starting values: pooled model taste params + log(μ) = 0 (μ = 1)
cf_pool <- coef(mnl_pool)
theta0  <- c(
  cf_pool[asc_names],
  cf_pool[price_names],
  0  # log(μ) = 0 → μ = 1
)

cat("Optimizing joint model with scale parameter...\n")
cat("  (This may take ~30–60 seconds)\n")

opt <- optim(
  par     = theta0,
  fn      = ll_joint_fn,
  df      = df0,
  method  = "BFGS",
  control = list(fnscale = -1, maxit = 10000, reltol = 1e-12),
  hessian = TRUE
)

if (opt$convergence != 0) {
  warning("Joint scale model convergence code: ", opt$convergence,
          " — consider re-running with different starting values")
}

ll_joint <- opt$value
mu_hat   <- exp(opt$par[11])

# Variance-covariance from negative inverse Hessian
vcov_joint <- tryCatch(
  solve(-opt$hessian),
  error = function(e) {
    warning("Hessian inversion failed: ", conditionMessage(e))
    matrix(NA_real_, 11, 11)
  }
)
se_joint <- sqrt(pmax(diag(vcov_joint), 0))

# Delta method: SE(μ) = |∂μ/∂log_μ| × SE(log_μ) = μ × SE(log_μ)
se_mu    <- mu_hat * se_joint[11]
t_mu     <- (mu_hat - 1) / se_mu  # t-stat for H0: μ = 1

cat("Joint model results:\n")
cat("  LL:              ", round(ll_joint, 3), "\n")
cat("  μ (scale):        ", round(mu_hat, 4), "\n")
cat("  SE(μ):            ", round(se_mu,  4), "\n")
cat("  t-stat (μ=1):     ", round(t_mu,   3), "\n")
cat("  Paper:  μ = 0.90, SE = 0.06\n\n")

# ── 7. Scale test (Likelihood Ratio) ──────────────────────────────────────────
# H0: μ = 1 AND equal taste params (joint restricted vs. two separate models)
# df = (params in A + params in B) - params in D = (10 + 10) - 11 = 9

lr_stat <- -2 * (ll_joint - ll_unres)
df_test <- (length(cf_nh) + length(cf_hyp)) - length(opt$par)
p_lr    <- pchisq(lr_stat, df = df_test, lower.tail = FALSE)

cat("=== Scale test (Swait-Louviere LR) ===\n")
cat(sprintf("  χ² = %.3f   df = %d   p = %.4f\n", lr_stat, df_test, p_lr))
cat("  Paper: χ² = 106.8, df = 4, p < 0.01\n\n")

# Simple Wald test: H0: μ = 1
p_wald_mu <- 2 * pnorm(-abs(t_mu))
cat(sprintf("  Wald test μ = 1: t = %.3f, p = %.4f\n\n", t_mu, p_wald_mu))

# ── 8. Format Table 3 ────────────────────────────────────────────────────────
# Layout matches Lusk & Schroeder (2004) Table 3.
# For each model: estimate + (standard error) for each parameter.

se_nh   <- sqrt(diag(vcov(mnl_nh)))
se_hyp  <- sqrt(diag(vcov(mnl_hyp)))

# Parameter labels (paper order: ASCs then price params)
param_labels <- c(
  "β GEN  (Generic)",
  "β GT   (Guaranteed Tender)",
  "β NAT  (Natural)",
  "β CHO  (USDA Choice)",
  "β CAB  (Certified Angus Beef)",
  "α GEN  (price)",
  "α GT   (price)",
  "α NAT  (price)",
  "α CHO  (price)",
  "α CAB  (price)"
)

# Extract estimates in order: ASCs then price params
ordered_names <- c(asc_names, price_names)

est_joint_v <- opt$par[1:10]
# Label joint estimates in same order as separate models
# theta order: [1:5] ASCs, [6:10] price
joint_ordered <- c(opt$par[1:5], opt$par[6:10])  # ASCs then price
se_joint_ord  <- c(se_joint[1:5], se_joint[6:10])

table3 <- data.frame(
  Parameter   = param_labels,
  NH_Est      = round(cf_nh[ordered_names],   4),
  NH_SE       = round(se_nh[ordered_names],   4),
  HY_Est      = round(cf_hyp[ordered_names],  4),
  HY_SE       = round(se_hyp[ordered_names],  4),
  JT_Est      = round(joint_ordered,          4),
  JT_SE       = round(se_joint_ord,           4),
  row.names   = NULL
)

cat("=== TABLE 3: MNL Estimates ===\n")
cat("Columns: Nonhypothetical | Hypothetical | Joint (with scale)\n\n")
print(kable(
  table3,
  format    = "simple",
  col.names = c("Parameter",
                "NH Est", "(SE)",
                "HY Est", "(SE)",
                "JT Est", "(SE)"),
  caption   = "Table 3. MNL Estimates (Lusk & Schroeder 2004)"
))

cat("\n--- Summary statistics ---\n")
cat(sprintf("%-22s %10s %10s %10s\n", "Statistic", "Nonhyp", "Hyp", "Joint"))
cat(sprintf("%-22s %10.3f %10.3f %10.3f\n",
            "Log-likelihood", ll_nh, ll_hyp, ll_joint))
cat(sprintf("%-22s %10d %10d %10d\n",
            "Observations",
            as.integer(sum(df0$treatment=="nonhyp" & df0$alt=="GEN")),
            as.integer(sum(df0$treatment=="hyp"    & df0$alt=="GEN")),
            as.integer(sum(df0$alt == "GEN"))))
cat(sprintf("%-22s %10s %10s %10.4f\n", "Scale (μ)", "—", "—", mu_hat))
cat(sprintf("%-22s %10s %10s %10.4f\n", "SE(μ)",     "—", "—", se_mu))
cat(sprintf("%-22s %10s %10s %10.3f\n", "Scale test χ²", "—", "—", lr_stat))
cat(sprintf("%-22s %10s %10s %10d\n",   "df",        "—", "—", df_test))

# ── 9. Save ───────────────────────────────────────────────────────────────────

results <- list(
  # mlogit model objects
  mnl_nh   = mnl_nh,
  mnl_hyp  = mnl_hyp,
  mnl_pool = mnl_pool,
  # Joint model optimisation output
  joint_opt  = opt,
  vcov_joint = vcov_joint,
  # Scalar summaries
  ll_nh    = ll_nh,
  ll_hyp   = ll_hyp,
  ll_pool  = ll_pool,
  ll_joint = ll_joint,
  mu_hat   = mu_hat,
  se_mu    = se_mu,
  lr_stat  = lr_stat,
  df_test  = df_test,
  p_lr     = p_lr,
  # Formatted table
  table3   = table3,
  # Coefficient name mapping (useful for downstream scripts)
  asc_names   = asc_names,
  price_names = price_names
)

saveRDS(results, file.path(output_dir, "table3_mnl.rds"))
write.csv(table3, file.path(output_dir, "table3_mnl.csv"), row.names = FALSE)
cat("\nSaved: output/table3_mnl.rds\n")
cat("Saved: output/table3_mnl.csv\n")
