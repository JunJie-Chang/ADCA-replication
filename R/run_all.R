# =============================================================================
# run_all.R  —  Run all five replication scripts from the terminal
#
# Usage (from any directory):
#   Rscript "/path/to/R/run_all.R"
#
# What it does:
#   1. Determines project_dir from this file's own path (no rstudioapi needed)
#   2. Mocks rstudioapi::getActiveDocumentContext so the five main scripts
#      can resolve their own paths correctly when sourced
#   3. Sources 01 → 02 → 03 → 04 → 05 in order
# =============================================================================

# ── 0. R options ─────────────────────────────────────────────────────────────
# mlogit uses deep recursion; increase the expression evaluation limit.
options(expressions = 500000)

# ── 1. Project root ───────────────────────────────────────────────────────────
# Hardcoded because the folder name contains Unicode + spaces, which causes
# commandArgs() path mangling when called via Rscript in the terminal.

project_dir <- paste0(
  "/Users/junjie/Library/CloudStorage/",
  "OneDrive-\u500b\u4eba(2)/AGEC/114-2/",
  "Applied Discrete Choice Analysis\uff08ADCA)"
)
cat("Project dir:", project_dir, "\n\n")

# Ensure output/ exists
output_dir <- file.path(project_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir)

# ── 1. Mock rstudioapi ────────────────────────────────────────────────────────
# Each script calls rstudioapi::getActiveDocumentContext()$path to find itself.
# We hijack the function to return a fake path inside R/ so that
# dirname(dirname(fake_path)) == project_dir.

fake_doc_path <- file.path(project_dir, "R", "fake_active_doc.R")

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  ns <- getNamespace("rstudioapi")
  if (bindingIsLocked("getActiveDocumentContext", ns)) {
    unlockBinding("getActiveDocumentContext", ns)
  }
  assign("getActiveDocumentContext",
         function() list(path = fake_doc_path, id = "", contents = "", selection = NULL),
         envir = ns)
  cat("rstudioapi::getActiveDocumentContext mocked.\n\n")
} else {
  # rstudioapi not installed: create a shim in the global environment
  # The :: operator will still hit the package, so install it first.
  stop("Please install rstudioapi: install.packages('rstudioapi')")
}

# ── 2. Source scripts in order ────────────────────────────────────────────────

scripts <- file.path(project_dir, "R", c(
  "01_data_prep.R",
  "02_table1_table2.R",
  "03_table3_mnl.R",
  "04_table4_adv.R",
  "05_table5_wtp.R"
))

for (s in scripts) {
  cat(strrep("=", 70), "\n")
  cat("RUNNING:", basename(s), "\n")
  cat(strrep("=", 70), "\n\n")
  tryCatch(
    source(s, echo = FALSE),
    error = function(e) {
      cat("\n*** ERROR in", basename(s), "***\n")
      cat(conditionMessage(e), "\n\n")
    }
  )
  cat("\n")
}

cat(strrep("=", 70), "\n")
cat("run_all.R complete.\n")
cat("Check output/ for saved .rds and .csv files.\n")
