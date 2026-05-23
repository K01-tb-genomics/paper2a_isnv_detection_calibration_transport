# =============================================================================
# 02_apply_filter_ladder.R
#
# Apply the inherited 3-tier filter ladder (Looser/Primary/Tighter) to both
# cohorts. For each tier x cohort, compute:
#   - n_iSNV per specimen
#   - iSNV_detected = 1{n_iSNV >= 1} per specimen
#
# These per-specimen outcomes are merged with the analytic frame and saved
# as a single sample-level data frame per cohort, ready for regression.
#
# Inputs:
#   data_derived/01_frames/florida_frame.rds   Sample-level metadata
#   data_derived/01_frames/ghana_m0_frame.rds
#   data_raw/florida_variants.rds              Per-variant rows
#   data_raw/ghana_variants.rds                Per-variant rows (all visits;
#                                              filtered to M0 here)
#   data_raw/paper1_thresholds/thr_*.rds       Inherited thresholds
#
# Outputs:
#   data_derived/02_ladder_applied/florida_ladder.rds
#   data_derived/02_ladder_applied/ghana_m0_ladder.rds
#   data_derived/02_ladder_applied/florida_variants_at_primary.rds
#   data_derived/02_ladder_applied/ghana_m0_variants_at_primary.rds
#       (filtered variant tables at Primary tier; needed for Test 5a and
#        Transport Test 2)
#   data_derived/02_ladder_applied/ladder_summary.txt
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

log_msg("\n--- 02_apply_filter_ladder.R ---")

# ---- Load frames and variant tables ----------------------------------------
fl_frame <- readRDS(file.path(PATHS$frames, "florida_frame.rds"))
gh_frame <- readRDS(file.path(PATHS$frames, "ghana_m0_frame.rds"))

fl_variants_all <- readRDS(file.path(PATHS$data_raw, "florida_gatkVariantTablesFormatted.rds"))
gh_variants_all <- readRDS(file.path(PATHS$data_raw, "ghana_gatkVariantTablesFormatted.rds"))

# Validate variant tables have required columns.
required_v_cols <- c("Sample", "POS", "REF", "ALT",
                     "DP", "AD1", "MAF", "FILTER")
require_cols(fl_variants_all, required_v_cols, "florida_variants.rds")
require_cols(gh_variants_all, required_v_cols, "ghana_variants.rds")

# ---- Restrict variant tables to analytic samples ---------------------------
fl_variants <- fl_variants_all[fl_variants_all$Sample %in% fl_frame$Sample, ,
                               drop = FALSE]
gh_variants <- gh_variants_all[gh_variants_all$Sample %in% gh_frame$Sample, ,
                               drop = FALSE]

log_msg("Florida variants: %d rows (%d distinct samples) restricted to frame.",
        nrow(fl_variants), length(unique(fl_variants$Sample)))
log_msg("Ghana M0 variants: %d rows (%d distinct samples) restricted to frame.",
        nrow(gh_variants), length(unique(gh_variants$Sample)))

# ---- Apply each tier and merge into the frame ------------------------------
# Note: helper count_isnv_per_sample() zero-fills samples with no calls so
# that every sample in the frame gets a row regardless of variant table
# coverage.
apply_ladder_to_cohort <- function(frame, variants, label) {
  log_msg("  [%s] Applying filter ladder...", label)
  out <- frame
  for (tier in TIERS) {
    thr  <- THR[[tier]]
    cnts <- count_isnv_per_sample(variants, thr, sample_ids = frame$Sample)
    cnt_col <- paste0("n_iSNV_",        tolower(tier))
    bin_col <- paste0("MAF_detected_",  tolower(tier))
    out[[cnt_col]] <- cnts$n_iSNV[match(out$Sample, cnts$Sample)]
    out[[bin_col]] <- as.integer(out[[cnt_col]] >= 1L)
    log_msg("    %s: n_iSNV>0 in %d/%d (%.1f%%); mean count=%.2f",
            tier, sum(out[[bin_col]]), nrow(out),
            100 * mean(out[[bin_col]]), mean(out[[cnt_col]]))
  }
  out
}

fl_ladder <- apply_ladder_to_cohort(fl_frame, fl_variants, "Florida")
gh_ladder <- apply_ladder_to_cohort(gh_frame, gh_variants, "Ghana M0")

# ---- Filtered variant tables at Primary (for Test 5a, Transport Test 2) ----
fl_variants_primary <- fl_variants[
  filter_variants_by_tier(fl_variants, THR$Primary), , drop = FALSE]
gh_variants_primary <- gh_variants[
  filter_variants_by_tier(gh_variants, THR$Primary), , drop = FALSE]

log_msg("Variants surviving Primary filter:")
log_msg("  Florida: %d (from %d rows)", nrow(fl_variants_primary),
        nrow(fl_variants))
log_msg("  Ghana M0: %d (from %d rows)", nrow(gh_variants_primary),
        nrow(gh_variants))

# ---- Save -------------------------------------------------------------------
save_with_provenance(
  fl_ladder, file.path(PATHS$ladder, "florida_ladder.rds"),
  inputs = list(florida_frame    = file.path(PATHS$frames, "florida_frame.rds"),
                florida_variants = file.path(PATHS$data_raw,
                                             "florida_variants.rds"),
                thresholds_dir   = INHERITED_THRESHOLDS_DIR),
  note = "Sample-level Florida frame with n_iSNV_<tier> and MAF_detected_<tier> appended for each of three inherited tiers."
)
save_with_provenance(
  gh_ladder, file.path(PATHS$ladder, "ghana_m0_ladder.rds"),
  inputs = list(ghana_frame    = file.path(PATHS$frames, "ghana_m0_frame.rds"),
                ghana_variants = file.path(PATHS$data_raw,
                                           "ghana_variants.rds"),
                thresholds_dir = INHERITED_THRESHOLDS_DIR),
  note = "Sample-level Ghana M0 frame with n_iSNV_<tier> and MAF_detected_<tier> appended."
)
save_with_provenance(
  fl_variants_primary,
  file.path(PATHS$ladder, "florida_variants_at_primary.rds"),
  inputs = list(florida_variants = file.path(PATHS$data_raw,
                                             "florida_variants.rds")),
  note = "Florida variant calls passing the Primary filter; PE/PPE excluded; FILTER=PASS only."
)
save_with_provenance(
  gh_variants_primary,
  file.path(PATHS$ladder, "ghana_m0_variants_at_primary.rds"),
  inputs = list(ghana_variants = file.path(PATHS$data_raw,
                                           "ghana_variants.rds")),
  note = "Ghana M0 variant calls passing the Primary filter; PE/PPE excluded; FILTER=PASS only."
)

# ---- Ladder summary log -----------------------------------------------------
ladder_summary <- function(df, label) {
  lines <- c(sprintf("=== %s (n = %d) ===", label, nrow(df)))
  for (tier in TIERS) {
    bin_col <- paste0("MAF_detected_", tolower(tier))
    cnt_col <- paste0("n_iSNV_",       tolower(tier))
    n_pos   <- sum(df[[bin_col]])
    n_total <- nrow(df)
    wci     <- wilson_ci(n_pos, n_total)
    lines <- c(lines, sprintf(
      "  %-7s  detected %d/%d (%.1f%%, Wilson 95%% CI %.1f-%.1f%%); mean n_iSNV = %.2f",
      tier, n_pos, n_total, 100 * wci["p"], 100 * wci["lower"],
      100 * wci["upper"], mean(df[[cnt_col]])))
  }
  lines
}

summary_lines <- c(
  sprintf("Filter ladder application summary (%s)", format(Sys.time())),
  "",
  "Inherited thresholds:",
  vapply(TIERS, function(t) sprintf(
    "  %-7s  DP>=%g  AD1>=%g  MAF in [%.2f, %.2f]",
    t, THR[[t]]$DP_min, THR[[t]]$AD1_min,
    THR[[t]]$MAF_min, THR[[t]]$MAF_max), character(1)),
  "",
  ladder_summary(fl_ladder, "Florida"),
  "",
  ladder_summary(gh_ladder, "Ghana M0")
)
writeLines(summary_lines, file.path(PATHS$ladder, "ladder_summary.txt"))

log_msg("02_apply_filter_ladder.R complete.\n")
