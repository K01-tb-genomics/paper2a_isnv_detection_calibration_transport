# =============================================================================
# 01_build_analytic_frames.R
#
# Build the sample-level analytic frames for both cohorts.
#
# Ghana:   restrict to M0 visit only (treatment-effect-free); apply
#          MIN_COV_INCLUDE depth floor; require non-missing core covariates.
#          Expected ~213 specimens from 95 of 97 patients (per Paper 1).
# Florida: apply MIN_COV_INCLUDE depth floor; require non-missing core
#          covariates; ensure cluster_id, run_id present and codable.
#          Expected up to 342 isolates across 5 transmission clusters.
#
# Inputs:
#   data_derived/00_metadata/florida_metadata.rds
#   data_derived/00_metadata/ghana_metadata.rds
#
# Outputs:
#   data_derived/01_frames/florida_frame.rds
#   data_derived/01_frames/ghana_m0_frame.rds
#   data_derived/01_frames/frame_summary.txt   (counts and exclusions log)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

log_msg("\n--- 01_build_analytic_frames.R ---")

fl_meta <- readRDS(file.path(PATHS$meta, "florida_metadata.rds"))
gh_meta <- readRDS(file.path(PATHS$meta, "ghana_metadata.rds"))

# ---- Florida frame ----------------------------------------------------------
fl_required_nonmissing <- c("HIV", "lineage", "Age", "Sex", "smear",
                            "run_id", "cluster_id", "depth_med", "callable_frac_20x",
                            "Year")

fl_n_start <- nrow(fl_meta)

# Stepwise filtering with explicit accounting at each step.
fl_steps <- list()

step_apply <- function(df, mask, label) {
  n_before <- nrow(df)
  df_after <- df[mask, , drop = FALSE]
  fl_steps[[label]] <<- list(n_before = n_before,
                             n_dropped = n_before - nrow(df_after),
                             n_after  = nrow(df_after))
  df_after
}

# (i) Coverage floor.
fl_meta_1 <- step_apply(fl_meta,
                        !is.na(fl_meta$depth_med) &
                          fl_meta$depth_med >= MIN_COV_INCLUDE,
                        sprintf("depth_med >= %dx", MIN_COV_INCLUDE))

# (ii) Non-missing core covariates.
fl_meta_2 <- step_apply(
  fl_meta_1,
  complete.cases(fl_meta_1[, fl_required_nonmissing, drop = FALSE]),
  "complete-cases on core covariates"
)

# (iii) Distinct samples (de-duplicate by Sample if any duplicates).
fl_meta_3 <- step_apply(
  fl_meta_2,
  !duplicated(fl_meta_2$Sample),
  "distinct Sample"
)

fl_frame <- fl_meta_3
log_msg("Florida frame: %d -> %d samples after frame construction.",
        fl_n_start, nrow(fl_frame))
log_msg("  Florida HIV+: %d (%.1f%%)",
        sum(fl_frame$HIV == 1L), 100 * mean(fl_frame$HIV == 1L))
log_msg("  Florida cluster_id distribution: %s",
        paste(sprintf("%s=%d", names(table(fl_frame$cluster_id)),
                      as.integer(table(fl_frame$cluster_id))),
              collapse = ", "))
log_msg("  Florida run_id distribution: %s",
        paste(sprintf("%s=%d", names(table(fl_frame$run_id)),
                      as.integer(table(fl_frame$run_id))),
              collapse = ", "))

# ---- Ghana M0 frame ---------------------------------------------------------
gh_required_nonmissing <- c("HIV", "lineage", "Age", "Sex", "patientId",
                            "depth_med", "callable_frac_20x")
gh_n_start <- nrow(gh_meta)
gh_steps <- list()

gh_step_apply <- function(df, mask, label) {
  n_before <- nrow(df)
  df_after <- df[mask, , drop = FALSE]
  gh_steps[[label]] <<- list(n_before = n_before,
                             n_dropped = n_before - nrow(df_after),
                             n_after  = nrow(df_after))
  df_after
}

# (i) Restrict to M0.
gh_meta_1 <- gh_step_apply(gh_meta,
                           !is.na(gh_meta$visit) & gh_meta$visit == "M0",
                           "visit == M0")

# (ii) Coverage floor.
gh_meta_2 <- gh_step_apply(
  gh_meta_1,
  !is.na(gh_meta_1$depth_med) & gh_meta_1$depth_med >= MIN_COV_INCLUDE,
  sprintf("depth_med >= %dx", MIN_COV_INCLUDE)
)

# (iii) Non-missing core covariates.
gh_meta_3 <- gh_step_apply(
  gh_meta_2,
  complete.cases(gh_meta_2[, gh_required_nonmissing, drop = FALSE]),
  "complete-cases on core covariates"
)

# (iv) Distinct samples.
gh_meta_4 <- gh_step_apply(
  gh_meta_3,
  !duplicated(gh_meta_3$Sample),
  "distinct Sample"
)

gh_frame <- gh_meta_4
log_msg("Ghana M0 frame: %d -> %d samples (%d distinct patients).",
        gh_n_start, nrow(gh_frame), length(unique(gh_frame$patientId)))
log_msg("  Ghana HIV+: %d (%.1f%%)",
        sum(gh_frame$HIV == 1L), 100 * mean(gh_frame$HIV == 1L))
log_msg("  Ghana mean specimens/patient: %.2f",
        nrow(gh_frame) / length(unique(gh_frame$patientId)))

# ---- Sanity checks against pre-spec expectations ---------------------------
# These are warnings, not hard fails: the analytic frame may legitimately
# differ slightly from manuscript text counts (e.g., if a covariate is
# missing for a few samples).
expect_msg <- function(actual, expected, label, tol = 0.10) {
  if (abs(actual - expected) / expected > tol) {
    warning(sprintf(
      "%s: observed %d, pre-spec text expected ~%d (%.0f%% deviation).",
      label, actual, expected, 100 * (actual - expected) / expected),
      call. = FALSE)
  }
}
expect_msg(nrow(fl_frame),                    342, "Florida frame size")
expect_msg(nrow(gh_frame),                    213, "Ghana M0 frame size")
expect_msg(length(unique(gh_frame$patientId)), 95, "Ghana M0 patient count")

# ---- Save -------------------------------------------------------------------
save_with_provenance(
  fl_frame, file.path(PATHS$frames, "florida_frame.rds"),
  inputs = list(florida_metadata = file.path(PATHS$meta,
                                             "florida_metadata.rds")),
  note = sprintf("Florida analytic frame; depth_med >= %dx; complete-case core covariates.",
                 MIN_COV_INCLUDE)
)
save_with_provenance(
  gh_frame, file.path(PATHS$frames, "ghana_m0_frame.rds"),
  inputs = list(ghana_metadata = file.path(PATHS$meta,
                                           "ghana_metadata.rds")),
  note = sprintf("Ghana M0-only analytic frame; depth_med >= %dx; complete-case core covariates.",
                 MIN_COV_INCLUDE)
)

# ---- Frame summary log ------------------------------------------------------
write_step_log <- function(steps, label) {
  c(sprintf("=== %s ===", label),
    vapply(names(steps), function(k) {
      s <- steps[[k]]
      sprintf("  %-50s n=%d -> %d  (dropped %d)",
              k, s$n_before, s$n_after, s$n_dropped)
    }, character(1)))
}

summary_lines <- c(
  sprintf("Frame summary written: %s", format(Sys.time())),
  "",
  write_step_log(fl_steps, "Florida"),
  "",
  write_step_log(gh_steps, "Ghana"),
  "",
  sprintf("Final Florida n = %d", nrow(fl_frame)),
  sprintf("Final Ghana M0 n = %d (across %d patients)",
          nrow(gh_frame), length(unique(gh_frame$patientId)))
)
writeLines(summary_lines, file.path(PATHS$frames, "frame_summary.txt"))

log_msg("01_build_analytic_frames.R complete.\n")
