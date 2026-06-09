# =============================================================================
# 06_iia_universal_logit.R
# IIA Test via Universal (Mother) Logit — Lusk & Schroeder (2004)
#
# Tests the Independence of Irrelevant Alternatives (IIA) assumption by
# estimating the universal logit of McFadden, Train, and Tye and performing
# a likelihood ratio test against the restricted MNL (Table 3).
#
# Universal logit utility (equation 4 in paper):
#   V_ij = β_j + Σ_{k=1}^{5} α_{jk} × P_{ik}
#   j ∈ {GEN, GT, NAT, CHO, CAB};  NONE = base (V_NONE ≡ 0)
#
# Parameters per model: 5 ASCs + 5 × 5 price params = 30
# MNL restriction: α_{jk} = 0 for j ≠ k  (20 constraints)
# IIA LR test df = 30 − 10 = 20
#
# Target LL values (paper p. 476):
#   Nonhyp: −1245.8 | Hyp: −597.7 | Joint: −1908.8
#
# Target IIA LR statistics χ²(20):
#   Nonhyp: 31.3 (p = 0.051) | Hyp: 59.5 (p < 0.01) | Joint: 66.9 (p < 0.01)
#
# Preference equality test (universal logit version):
#   χ²(30) = 130.7  (p < 0.01)
#
# NOTE: Run from RStudio (project folder contains Unicode characters).
#
# Requires: output/df_mlogit.rds, output/price_matrix.rds,
#           output/table3_mnl.rds
# Saves:    output/iia_ul.rds, output/iia_ul.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(mlogit)
library(knitr)

# ── 0. Paths and data ─────────────────────────────────────────────────────────

project_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
output_dir  <- file.path(project_dir, "output")

df0          <- readRDS(file.path(output_dir, "df_mlogit.rds"))
price_matrix <- readRDS(file.path(output_dir, "price_matrix.rds"))
mnl_results  <- readRDS(file.path(output_dir, "table3_mnl.rds"))

# MNL log-likelihood values from Table 3
ll_mnl_nh  <- mnl_results$ll_nh     # paper: −1261.4
ll_mnl_hyp <- mnl_results$ll_hyp    # paper:  −627.4
ll_mnl_jt  <- mnl_results$ll_joint  # paper: −1942.2

cat("MNL LL values loaded:\n")
cat(sprintf("  Nonhyp: %8.3f  (paper: −1261.4)\n", ll_mnl_nh))
cat(sprintf("  Hyp:    %8.3f  (paper:  −627.4)\n", ll_mnl_hyp))
cat(sprintf("  Joint:  %8.3f  (paper: −1942.2)\n", ll_mnl_jt))

# ── 1. Add all 5 scenario-level prices to every row ───────────────────────────
# df_mlogit only has price_alt (own-price). The universal logit needs all 5
# prices in each row so cross-price effects can be computed.

alts <- c("GEN", "GT", "NAT", "CHO", "CAB", "NONE")

df0 <- df0 |>
  mutate(
    alt       = factor(as.character(alt), levels = alts),
    price_alt = replace_na(price_alt, 0)
  ) |>
  left_join(
    price_matrix |>
      select(scenario, P_GEN = GEN, P_GT = GT,
             P_NAT = NAT, P_CHO = CHO, P_CAB = CAB),
    by = "scenario"
  )

# ── 2. Create 25 cross-price interaction dummies ──────────────────────────────
# p_{j}_{k} = P_k × I(alt == j)
# When alt == NONE: all dummies = 0  →  V_NONE = 0  ✓

steak_alts <- c("GEN", "GT", "NAT", "CHO", "CAB")
price_cols <- c("P_GEN", "P_GT", "P_NAT", "P_CHO", "P_CAB")

for (j in steak_alts) {
  for (ki in seq_along(steak_alts)) {
    k     <- steak_alts[ki]
    pc    <- price_cols[ki]
    vname <- paste0("p_", j, "_", k)
    df0[[vname]] <- ifelse(as.character(df0$alt) == j, df0[[pc]], 0)
  }
}

# Variable ordering: j outer (each = 5), k inner — matches theta indexing below
ul_price_names <- paste0("p_", rep(steak_alts, each = 5), "_",
                          rep(steak_alts, 5))

cat(sprintf("\n25 cross-price dummies created (%d total price vars)\n",
            length(ul_price_names)))
cat(sprintf("Sanity check — p_GEN_GEN nonzero rows (nonhyp): %d  (expected: 1072)\n",
    sum(df0$p_GEN_GEN != 0 & df0$treatment == "nonhyp")))

# ── 3. mlogit formula and data objects ────────────────────────────────────────
# 25 generic-coef terms (alt-varying) + | 1 for 5 ASCs (NONE = reference)
# → 30 parameters per model

ul_terms <- paste(ul_price_names, collapse = " + ")
f_ul     <- as.formula(paste("chosen ~", ul_terms, "| 1"))

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
mdata_all <- make_mdata(df0)

# ── 4. Estimate separate universal logit models ───────────────────────────────

cat("\nEstimating nonhypothetical universal logit (30 params) ...\n")
cat("  (expect ~60–120 s)\n")
ul_nh <- tryCatch(
  mlogit(f_ul, data = mdata_nh, reflevel = "NONE", iterlim = 5000),
  error = function(e) { warning("ul_nh failed: ", e$message); NULL }
)
ll_ul_nh <- if (!is.null(ul_nh)) as.numeric(logLik(ul_nh)) else NA_real_
cat(sprintf("  LL: %8.3f  (paper: −1245.8)\n", ll_ul_nh))

cat("\nEstimating hypothetical universal logit (30 params) ...\n")
ul_hyp <- tryCatch(
  mlogit(f_ul, data = mdata_hyp, reflevel = "NONE", iterlim = 5000),
  error = function(e) { warning("ul_hyp failed: ", e$message); NULL }
)
ll_ul_hyp <- if (!is.null(ul_hyp)) as.numeric(logLik(ul_hyp)) else NA_real_
cat(sprintf("  LL: %8.3f  (paper:  −597.7)\n", ll_ul_hyp))

# ── 5. Pooled universal logit (starting values for joint model) ───────────────

cat("\nEstimating pooled universal logit (starting values for joint) ...\n")
ul_pool <- tryCatch(
  mlogit(f_ul, data = mdata_all, reflevel = "NONE", iterlim = 5000),
  error = function(e) { warning("ul_pool failed: ", e$message); NULL }
)
ll_ul_pool <- if (!is.null(ul_pool)) as.numeric(logLik(ul_pool)) else NA_real_
cat(sprintf("  Pooled LL: %8.3f\n", ll_ul_pool))

# ── 6. Joint universal logit with scale parameter (custom optim) ───────────────
# θ layout: [β_GEN … β_CAB (5),  α_{GEN,GEN} … α_{CAB,CAB} (25),  log_μ (1)]
#           = 31 parameters total
#
# θ[6:30] are stored row-major: j outer (each=5), k inner.
# alpha_mat[j, k] = θ[5 + (j-1)*5 + k]  for j,k = 1..5 (steak alts only)
#
# V_ij^{NH} = μ × (β_j + Σ_k α_{jk} P_{ik})
# V_ij^{H}  =       β_j + Σ_k α_{jk} P_{ik}
# V_{NONE}  = 0  (reference)

ll_ul_joint_fn <- function(theta, df) {
  beta      <- c(theta[1:5], 0)                             # 6-element ASC vector
  alpha_mat <- rbind(
    matrix(theta[6:30], nrow = 5, ncol = 5, byrow = TRUE),  # 5×5 steak params
    rep(0, 5)                                                # row 6 = NONE
  )
  mu <- exp(theta[31])

  alt_idx <- match(as.character(df$alt), alts)  # 1–6

  prices <- cbind(df$P_GEN, df$P_GT, df$P_NAT, df$P_CHO, df$P_CAB)
  prices[is.na(prices)] <- 0

  V <- beta[alt_idx] + rowSums(alpha_mat[alt_idx, ] * prices)

  V[df$treatment == "nonhyp"] <- mu * V[df$treatment == "nonhyp"]

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

# Starting values from pooled universal logit (or fallback to zeros)
if (!is.null(ul_pool)) {
  cf_pool    <- coef(ul_pool)
  ul_asc_nm  <- setdiff(names(cf_pool), ul_price_names)

  asc_vals <- sapply(steak_alts, function(j) {
    nm <- ul_asc_nm[grepl(paste0("\\b", j, "\\b"), ul_asc_nm)]
    if (length(nm) > 0) cf_pool[nm[1]] else 0
  })
  price_vals <- cf_pool[ul_price_names]
} else {
  asc_vals   <- rep(0, 5)
  price_vals <- rep(0, 25)
}

theta0_jt <- c(asc_vals, price_vals, 0)  # log_mu = 0  →  mu = 1

cat("\nOptimizing joint universal logit with scale (31 params) ...\n")
cat("  (expect ~3–5 min)\n")

opt_ul <- optim(
  par     = theta0_jt,
  fn      = ll_ul_joint_fn,
  df      = df0,
  method  = "BFGS",
  control = list(fnscale = -1, maxit = 20000, reltol = 1e-12),
  hessian = FALSE
)

if (opt_ul$convergence != 0)
  warning("Joint UL optim convergence code: ", opt_ul$convergence,
          " — consider different starting values or tighter tolerances.")

ll_ul_jt <- opt_ul$value
mu_ul    <- exp(opt_ul$par[31])

cat(sprintf("  Joint UL LL: %8.3f  (paper: −1908.8)\n", ll_ul_jt))
cat(sprintf("  μ (scale):   %8.4f\n", mu_ul))

# ── 7. IIA likelihood ratio tests ─────────────────────────────────────────────
# H0: all 20 cross-price effects = 0  (MNL is the restricted model)
# LR = 2 × (LL_UL − LL_MNL),  df = 30 − 10 = 20

df_iia <- 20

chi2_nh  <- 2 * (ll_ul_nh  - ll_mnl_nh)
chi2_hyp <- 2 * (ll_ul_hyp - ll_mnl_hyp)
chi2_jt  <- 2 * (ll_ul_jt  - ll_mnl_jt)

p_nh  <- pchisq(chi2_nh,  df = df_iia, lower.tail = FALSE)
p_hyp <- pchisq(chi2_hyp, df = df_iia, lower.tail = FALSE)
p_jt  <- pchisq(chi2_jt,  df = df_iia, lower.tail = FALSE)

cat("\n")
cat("=== IIA Test: Universal Logit vs. MNL  [χ²(20)] ===\n")
cat(sprintf("%-14s  %9s  %9s  %6s  %7s  %9s\n",
            "Model", "LL(MNL)", "LL(UL)", "χ²", "p-val", "Paper χ²"))
cat(sprintf("%-14s  %9.3f  %9.3f  %6.2f  %7.4f  %9s\n",
            "Joint",    ll_mnl_jt,  ll_ul_jt,  chi2_jt,  p_jt,  "66.9"))
cat(sprintf("%-14s  %9.3f  %9.3f  %6.2f  %7.4f  %9s\n",
            "Nonhyp",   ll_mnl_nh,  ll_ul_nh,  chi2_nh,  p_nh,  "31.3"))
cat(sprintf("%-14s  %9.3f  %9.3f  %6.2f  %7.4f  %9s\n",
            "Hyp",      ll_mnl_hyp, ll_ul_hyp, chi2_hyp, p_hyp, "59.5"))

# ── 8. Preference equality test (universal logit version) ─────────────────────
# H0: taste params equal across treatments, allowing for scale difference
# χ² = 2 × (LL_nh_ul + LL_hyp_ul − LL_joint_ul),  df = 30

chi2_eq <- 2 * (ll_ul_nh + ll_ul_hyp - ll_ul_jt)
df_eq   <- 30
p_eq    <- pchisq(chi2_eq, df = df_eq, lower.tail = FALSE)

cat("\n=== Preference Equality Test (Universal Logit) ===\n")
cat(sprintf("  χ²(%d) = %.3f   p = %.4f\n", df_eq, chi2_eq, p_eq))
cat("  Paper: χ²(30) = 130.7, p < 0.01\n")

# ── 9. Summary table ──────────────────────────────────────────────────────────

summary_tbl <- data.frame(
  Model      = c("Joint", "Nonhypothetical", "Hypothetical"),
  LL_MNL     = round(c(ll_mnl_jt,  ll_mnl_nh,  ll_mnl_hyp),  1),
  LL_UL      = round(c(ll_ul_jt,   ll_ul_nh,   ll_ul_hyp),   1),
  Chi2_df20  = round(c(chi2_jt,    chi2_nh,    chi2_hyp),    1),
  p_value    = round(c(p_jt,       p_nh,       p_hyp),       3),
  Paper_Chi2 = c(66.9, 31.3, 59.5),
  Paper_p    = c("<0.01", "0.051", "<0.01"),
  stringsAsFactors = FALSE
)

cat("\n=== IIA Test Summary Table ===\n")
print(kable(
  summary_tbl,
  format  = "simple",
  col.names = c("Model", "LL(MNL)", "LL(UL)", "χ²(20)",
                "p-value", "Paper χ²", "Paper p"),
  caption = "IIA Test via Universal Logit (Lusk & Schroeder 2004, p. 476)"
))

# ── 10. Save ──────────────────────────────────────────────────────────────────

ul_results <- list(
  ul_nh      = ul_nh,
  ul_hyp     = ul_hyp,
  ul_pool    = ul_pool,
  opt_ul_jt  = opt_ul,
  ll_ul_nh   = ll_ul_nh,
  ll_ul_hyp  = ll_ul_hyp,
  ll_ul_jt   = ll_ul_jt,
  mu_ul      = mu_ul,
  chi2_nh    = chi2_nh,   chi2_hyp  = chi2_hyp,  chi2_jt   = chi2_jt,
  p_nh       = p_nh,      p_hyp     = p_hyp,      p_jt      = p_jt,
  chi2_eq    = chi2_eq,   df_eq     = df_eq,       p_eq      = p_eq,
  summary_tbl = summary_tbl
)

saveRDS(ul_results, file.path(output_dir, "iia_ul.rds"))
write.csv(summary_tbl, file.path(output_dir, "iia_ul.csv"), row.names = FALSE)
cat("\nSaved: output/iia_ul.rds\n")
cat("Saved: output/iia_ul.csv\n")
