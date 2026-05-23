# =============================================================================
# 00_pipeline_config.R
#
# Single source of truth for canonical specification flags, paths, and the
# inherited Paper 1 filter thresholds. Sourced at the top of every pipeline
# script. The pre-specification snapshot (see README.md) locks every flag
# below; deviations must be disclosed in the manuscript.
#
# REVISED 2026-05-17: Paper 2 split into Paper 2a (filter transportability)
# and Paper 2b (Edwards-2023 operationalization).
# Added per-paper output paths. Analysis scripts 00-12
# remain unchanged; their outputs land in shared data_derived/ and are
# consumed by paper-specific output scripts 13a/14a (2a) and 13b/14b (2b).
#
# Canonical specification for Paper 2 (LOCKED at pre-spec commit):
#   - Outcome:        iSNV_detected = 1{n_iSNV >= 1}, per specimen at each tier
#   - Filter ladder:  Looser / Primary / Tighter (inherited from Paper 1)
#   - Frame:          Ghana M0-only (213 specimens / 95 patients);
#                     Florida 5 transmission clusters (342 isolates),
#                     median coverage >= 50x in both
#   - Adjustment:
#       Florida: HIV + lineage + Age + Sex + smear + run_id
#       Ghana:   HIV + lineage + Age + Sex
#   - Standard errors: cluster-robust sandwich (glm + sandwich::vcovCL HC1)
#       Florida: cluster_id (G=5);  Ghana: patientId
#   - Monotonicity test: chi-bar-squared (Bartholomew) on |log PR|,
#       B = 2000 stratified cluster bootstrap, sign rule = absolute value
#   - Transport Test 1: Florida Primary detection rate in [6.5%, 26.5%]
#   - Transport Test 2: within-cluster Jaccard ~ min(DP_A, DP_B), positive slope
#   - Test 5a: variant-level HIV x log2(DP) interaction at Primary
#   - Test 5b: Florida calendar-year placebo + Cochran-Armitage trend
# =============================================================================

# ---- Canonical flags (locked at pre-spec) ----------------------------------
# `exists()` guards let an outer orchestrator override before re-sourcing,
# matching Paper 1's pattern.
if (!exists("M0_ONLY",         envir = .GlobalEnv, inherits = FALSE)) M0_ONLY         <- TRUE
if (!exists("MIN_COV_INCLUDE", envir = .GlobalEnv, inherits = FALSE)) MIN_COV_INCLUDE <- 50L

# ---- Random seed (LOCKED at pre-spec) --------------------------------------
# This seed is used for all stochastic steps including the chi-bar-squared
# cluster bootstrap. Must match the seed declared in README.md.
if (!exists("SEED", envir = .GlobalEnv, inherits = FALSE)) SEED <- 20260509L

# ---- Inherited Paper 1 thresholds ------------------------------------------
# Per-spec README, Paper 1's thr_*.rds files are committed under data_raw/.
# Each is a list with elements: DP_min, AD1_min, MAF_min, MAF_max, tier.
INHERITED_THRESHOLDS_DIR <- file.path("data_raw", "paper1_thresholds")

# Locked threshold values (asserted against on-disk files at load time).
# These are the values pre-specified in methods_rewrite_v1.md and README.md.
LOCKED_THR <- list(
  Looser  = list(DP_min = 40L,  AD1_min = 3L, MAF_min = 0.02, MAF_max = 0.50),
  Primary = list(DP_min = 60L,  AD1_min = 3L, MAF_min = 0.02, MAF_max = 0.50),
  Tighter = list(DP_min = 100L, AD1_min = 6L, MAF_min = 0.02, MAF_max = 0.50)
)
TIERS <- c("Looser", "Primary", "Tighter")

# ---- Bootstrap settings ----------------------------------------------------
N_BOOTSTRAP_MONO <- 2000L  # Chi-bar-squared monotonicity test
N_BOOTSTRAP_GAP  <- 2000L  # If used by sensitivity (gap-widening), kept for parity

# ---- Outcome-ladder thresholds (sensitivity, not pre-spec primary) ----------
OUTCOME_LADDER_K <- c(1L, 2L)  # MAF_detected at k = 1 vs k >= 2 iSNVs

# ---- Restriction-based MeasAdj sensitivity ---------------------------------
HIGH_QC_DEPTH_MED <- 80L
HIGH_QC_CALLABLE  <- 0.95

# ---- Depth-floor sweep sensitivity -----------------------------------------
DEPTH_FLOOR_SWEEP <- seq(50L, 150L, by = 10L)

# ---- Output paths ----------------------------------------------------------
# Shared paths (analysis intermediate outputs consumed by both papers).
# Per-paper paths (manuscript tables and figures) added at end.
PATHS <- list(
  data_raw          = "data_raw",
  inherited_thr     = INHERITED_THRESHOLDS_DIR,
  meta              = file.path("data_derived", "00_metadata"),
  variants          = file.path("data_derived", "00_variants"),
  frames            = file.path("data_derived", "01_frames"),
  ladder            = file.path("data_derived", "02_ladder_applied"),
  primary_pr        = file.path("data_derived", "03_primary_pr"),
  monotonicity      = file.path("data_derived", "04_monotonicity"),
  transport_t1      = file.path("data_derived", "05_transport_t1"),
  transport_t2      = file.path("data_derived", "06_transport_t2"),
  test5a            = file.path("data_derived", "07_test5a"),
  test5b            = file.path("data_derived", "08_test5b"),
  sens_qc_restrict  = file.path("data_derived", "09_sens_qc_restriction"),
  sens_outcome_lad  = file.path("data_derived", "10_sens_outcome_ladder"),
  sens_depth_sweep  = file.path("data_derived", "11_sens_depth_sweep"),
  sens_year         = file.path("data_derived", "09_sens_year"),
  sens_qc           = file.path("data_derived", "10_sens_qc"),
  sens_outcome      = file.path("data_derived", "11_sens_outcome"),
  # IPW paths (Paper 2b, Florida-only selection-bias analysis).
  ipw_selection     = file.path("data_derived", "18_ipw_selection"),
  ipw_florida       = file.path("data_derived", "19_ipw_florida"),
  # Operating characteristics path (Paper 2b).
  op_chars          = file.path("data_derived", "13_operating_characteristics"),
  # Per-paper output paths (manuscript deliverables).
  paper2a_tables    = file.path("outputs", "paper2a", "tables"),
  paper2a_figures   = file.path("outputs", "paper2a", "figures"),
  paper2b_tables    = file.path("outputs", "paper2b", "tables"),
  paper2b_figures   = file.path("outputs", "paper2b", "figures"),
  # Legacy paths (kept for backward compatibility with 13/14/15 if still used).
  tables            = file.path("outputs", "tables"),
  figures           = file.path("outputs", "figures"),
  supplemental      = file.path("outputs", "supplemental"),
  logs              = file.path("outputs", "orchestrator_logs")
)

for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- Validate inherited thresholds against the locked spec -----------------
# Hard-fails if Paper 1 thr_*.rds files differ from the values pre-specified
# in this paper. This is the gate that prevents silent drift between Paper 1
# revisions and Paper 2's pre-registered application.
validate_inherited_thresholds <- function() {
  if (!dir.exists(INHERITED_THRESHOLDS_DIR)) {
    stop(sprintf(
      "Inherited thresholds directory not found: %s\n",
      INHERITED_THRESHOLDS_DIR),
      "Copy thr_looser.rds, thr_primary.rds, thr_tighter.rds from Paper 1's\n",
      "data_derived/01_calibration_snv_only_ppe_retained_lex/ into this dir.")
  }
  thr_files <- c(Looser  = "thr_looser.rds",
                 Primary = "thr_primary.rds",
                 Tighter = "thr_tighter.rds")
  thr <- list()
  for (nm in names(thr_files)) {
    fp <- file.path(INHERITED_THRESHOLDS_DIR, thr_files[[nm]])
    if (!file.exists(fp)) stop("Missing inherited threshold file: ", fp)
    thr[[nm]] <- readRDS(fp)
    locked <- LOCKED_THR[[nm]]
    for (key in c("DP_min", "AD1_min", "MAF_min", "MAF_max")) {
      observed <- thr[[nm]][[key]]
      expected <- locked[[key]]
      if (is.null(observed) || !isTRUE(all.equal(observed, expected,
                                                 tolerance = 1e-9))) {
        stop(sprintf(
          "Inherited %s tier %s mismatch: file = %s, locked = %s.\n",
          nm, key, toString(observed), toString(expected)),
          "Pre-specification was locked on the values in LOCKED_THR; ",
          "if Paper 1 has been revised, the pre-spec is broken and the ",
          "deviation must be disclosed in the manuscript.")
      }
    }
  }
  thr
}

# Eager validation: any pipeline script that sources this config will fail
# fast on threshold drift.
THR <- validate_inherited_thresholds()

# ---- Diagnostic banner -----------------------------------------------------
banner <- function() {
  cat(strrep("=", 78), "\n", sep = "")
  cat("Paper 2 pipeline specification\n")
  cat(sprintf("  M0_ONLY            = %s\n", M0_ONLY))
  cat(sprintf("  MIN_COV_INCLUDE    = %d\n", MIN_COV_INCLUDE))
  cat(sprintf("  SEED               = %d\n", SEED))
  cat(sprintf("  N_BOOTSTRAP_MONO   = %d\n", N_BOOTSTRAP_MONO))
  cat("Inherited thresholds (validated against LOCKED_THR):\n")
  for (nm in TIERS) {
    t <- THR[[nm]]
    cat(sprintf("  %-7s  DP>=%g  AD1>=%g  MAF in [%.2f, %.2f]\n",
                nm, t$DP_min, t$AD1_min, t$MAF_min, t$MAF_max))
  }
  cat(strrep("=", 78), "\n", sep = "")
}
