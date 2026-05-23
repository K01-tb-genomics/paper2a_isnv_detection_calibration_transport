# =============================================================================
# 05_transport_test_1_band.R
#
# Transport Test 1 (per-specimen prevalence band).
#
# Pre-specified criterion (LOCKED at pre-spec commit):
#   Florida's per-specimen detection rate at the Primary filter tier
#   (point estimate) must fall within [6.5%, 26.5%], i.e., 16.5% +/- 10 pp
#   around the Ghana per-patient any-visit prevalence reported in Paper 1
#   (16/97). Falling inside the band declares transport successful.
#   Falling below indicates Primary is too strict for Florida's depth
#   regime. Falling above indicates Primary is too permissive at Florida's
#   depth (admitting noise). Both directions are interpretable failures.
#
# Inputs:
#   data_derived/02_ladder_applied/florida_ladder.rds
#
# Outputs:
#   data_derived/05_transport_t1/transport_t1.rds
#   data_derived/05_transport_t1/transport_t1_summary.txt
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

log_msg("\n--- 05_transport_test_1_band.R ---")

# ---- Pre-specified band (LOCKED) -------------------------------------------
GHANA_ANCHOR <- 0.165   # 16/97 per-patient any-visit prevalence (Paper 1)
DELTA_PP     <- 0.100   # +/- 10 percentage points
BAND_LOWER   <- GHANA_ANCHOR - DELTA_PP   # 0.065
BAND_UPPER   <- GHANA_ANCHOR + DELTA_PP   # 0.265

# ---- Florida Primary detection rate ----------------------------------------
fl_ladder <- readRDS(file.path(PATHS$ladder, "florida_ladder.rds"))

n_total      <- nrow(fl_ladder)
n_events     <- sum(fl_ladder$MAF_detected_primary)
rate         <- n_events / n_total
wci          <- wilson_ci(n_events, n_total)

# ---- Decision rule ---------------------------------------------------------
inside_band <- (rate >= BAND_LOWER) && (rate <= BAND_UPPER)
direction <- if (inside_band) {
  "inside the band"
} else if (rate < BAND_LOWER) {
  "below the band lower bound"
} else {
  "above the band upper bound"
}
declared <- if (inside_band) {
  "TRANSPORT SUCCESSFUL"
} else {
  "TRANSPORT FAILED"
}

# ---- Logging ---------------------------------------------------------------
log_msg("Florida Primary tier detection rate:")
log_msg("  events / n        = %d / %d", n_events, n_total)
log_msg("  point estimate    = %.4f (%.2f%%)", rate, 100 * rate)
log_msg("  Wilson 95%% CI     = [%.4f, %.4f] (%.2f%% - %.2f%%)",
        wci["lower"], wci["upper"], 100 * wci["lower"], 100 * wci["upper"])
log_msg("Pre-specified band:")
log_msg("  Anchor (Paper 1)  = %.4f (16/97 any-visit per-patient)",
        GHANA_ANCHOR)
log_msg("  Tolerance delta   = +/- %.1f pp", 100 * DELTA_PP)
log_msg("  Band              = [%.4f, %.4f]", BAND_LOWER, BAND_UPPER)
log_msg("Decision:")
log_msg("  Florida point estimate is %s.", direction)
log_msg("  Outcome: %s", declared)

interpretation <- if (inside_band) {
  "Florida Primary detection behaves consistently with the calibration regime within the pre-specified tolerance."
} else if (rate < BAND_LOWER) {
  "Florida Primary detection is suppressed below the band, suggesting the Primary tier is too strict at Florida's depth regime - real signal is being filtered out as low-quality calls."
} else {
  "Florida Primary detection exceeds the band, suggesting the Primary tier is too permissive at Florida's depth regime - admitting depth-correlated noise as iSNVs."
}
log_msg("Interpretation:")
log_msg("  %s", interpretation)

# ---- Save ------------------------------------------------------------------
result <- list(
  test                = "Transport Test 1 (per-specimen prevalence band)",
  florida_n           = n_total,
  florida_n_events    = n_events,
  florida_rate        = rate,
  florida_wilson_ci   = wci,
  ghana_anchor        = GHANA_ANCHOR,
  delta_pp            = DELTA_PP,
  band_lower          = BAND_LOWER,
  band_upper          = BAND_UPPER,
  inside_band         = inside_band,
  direction           = direction,
  declared            = declared,
  interpretation      = interpretation
)

save_with_provenance(
  result, file.path(PATHS$transport_t1, "transport_t1.rds"),
  inputs = list(florida_ladder = file.path(PATHS$ladder,
                                           "florida_ladder.rds")),
  note = sprintf("Transport Test 1: Florida Primary rate = %.4f vs band [%.4f, %.4f]; %s",
                 rate, BAND_LOWER, BAND_UPPER, declared)
)

# Human-readable summary
summary_lines <- c(
  sprintf("Transport Test 1 summary  (%s)", format(Sys.time())),
  "",
  sprintf("Florida Primary detection: %d / %d = %.4f (%.2f%%)",
          n_events, n_total, rate, 100 * rate),
  sprintf("  Wilson 95%% CI: [%.4f, %.4f]", wci["lower"], wci["upper"]),
  "",
  sprintf("Pre-specified band: [%.4f, %.4f]", BAND_LOWER, BAND_UPPER),
  sprintf("  Anchor: %.4f (Paper 1 16/97 any-visit per-patient prevalence)",
          GHANA_ANCHOR),
  sprintf("  Tolerance: +/- %.1f percentage points", 100 * DELTA_PP),
  "",
  sprintf("Decision: %s (%s)", declared, direction),
  "",
  "Interpretation:",
  paste0("  ", interpretation)
)
writeLines(summary_lines,
           file.path(PATHS$transport_t1, "transport_t1_summary.txt"))

log_msg("05_transport_test_1_band.R complete.\n")
