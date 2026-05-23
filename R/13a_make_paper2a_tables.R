# =============================================================================
# 13a_make_paper2a_tables.R
#
# Builds Paper 2a (filter transportability) manuscript tables.
#
# Tables produced:
#   Table 1  Cohort characteristics with depth regime, cluster structure,
#            and pair structure. Florida-only (the transport-application
#            cohort). Ghana calibration anchor cited via Paper 1 in narrative.
#   Table 2  Transport test summary: T1 prevalence band, T2 Jaccard~depth
#            slope, both-empty complementary slope, hurdle decomposition
#            (extensive + intensive margins), drop-cluster-3 sensitivity.
#
# Outputs:
#   outputs/paper2a/tables/table1_cohort_structure.{csv,txt}
#   outputs/paper2a/tables/table2_transport_summary.{csv,txt}
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

log_msg("\n--- 13a_make_paper2a_tables.R ---")

TABLES_DIR <- PATHS$paper2a_tables
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Formatting helpers ----------------------------------------------------
fmt_n_pct <- function(n, total) {
  if (is.na(n) || is.na(total) || total == 0) return("-")
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}
fmt_median_iqr <- function(x, d = 0) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("-")
  q <- quantile(x, c(0.25, 0.5, 0.75))
  sprintf(sprintf("%%.%df (%%.%df-%%.%df)", d, d, d), q[2], q[1], q[3])
}
fmt_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
write_table <- function(df, name) {
  csv_path <- file.path(TABLES_DIR, paste0(name, ".csv"))
  txt_path <- file.path(TABLES_DIR, paste0(name, ".txt"))
  write.csv(df, csv_path, row.names = FALSE, na = "")
  capture.output(print(df, row.names = FALSE), file = txt_path)
  log_msg("  Wrote %s.csv and %s.txt", name, name)
}

# ============================================================================
# Table 1 - Florida cohort structure (transport-application cohort)
#
# Note: Ghana calibration anchor (n_pairs = 169, 67 patients; per-patient
# prevalence 16/97 = 16.5%; pair-level Jaccard slope 0.029 per +10x) is cited
# from Paper 1 in the manuscript narrative rather than reproduced here.
# ============================================================================
make_table1 <- function() {
  log_msg("\nBuilding Table 1 (Florida cohort structure)...")
  fl <- readRDS(file.path(PATHS$ladder, "florida_ladder.rds"))
  t2 <- readRDS(file.path(PATHS$transport_t2, "transport_t2.rds"))
  pairs <- readRDS(file.path(PATHS$transport_t2, "transport_t2_pairs.rds"))

  # Cluster sizes
  cluster_tbl <- table(fl$cluster_id)
  cluster_sizes <- paste(sort(as.integer(cluster_tbl), decreasing = FALSE),
                         collapse = ", ")
  n_clusters <- length(cluster_tbl)

  # Year range
  if ("year" %in% names(fl)) {
    yr_range <- sprintf("%d-%d", min(fl$year, na.rm = TRUE),
                        max(fl$year, na.rm = TRUE))
  } else {
    yr_range <- "[PLACEHOLDER: year range not in ladder file]"
  }

  # Run_id distribution
  if ("run_id" %in% names(fl)) {
    run_tbl <- table(fl$run_id)
    run_str <- paste(sprintf("%s = %d", names(run_tbl), as.integer(run_tbl)),
                     collapse = "; ")
  } else {
    run_str <- "-"
  }

  # Pair-structure counts (from transport_t2 outputs)
  n_pairs_total      <- nrow(pairs)
  n_pairs_evaluable  <- sum(!is.na(pairs$jaccard))
  n_pairs_both_empty <- if ("both_empty" %in% names(pairs)) {
    sum(pairs$both_empty, na.rm = TRUE)
  } else {
    n_pairs_total - n_pairs_evaluable
  }

  rows <- list(
    data.frame(Characteristic = "Specimens, n",
               Value = sprintf("%d", nrow(fl))),
    data.frame(Characteristic = "Patients, n (1 specimen/patient)",
               Value = sprintf("%d", nrow(fl))),
    data.frame(Characteristic = "Year range",
               Value = yr_range),
    data.frame(Characteristic = "Transmission clusters, n",
               Value = sprintf("%d", n_clusters)),
    data.frame(Characteristic = "Cluster sizes",
               Value = cluster_sizes),
    data.frame(Characteristic = "Median sample-level depth (x), median (Q1-Q3)",
               Value = fmt_median_iqr(fl$depth_med)),
    data.frame(Characteristic = "Callable fraction, median (Q1-Q3)",
               Value = fmt_median_iqr(fl$callable, 2)),
    data.frame(Characteristic = "Sequencing runs",
               Value = run_str),
    data.frame(Characteristic = "Within-cluster pairs, total",
               Value = format(n_pairs_total, big.mark = ",")),
    data.frame(Characteristic = "Within-cluster pairs, Jaccard-evaluable",
               Value = format(n_pairs_evaluable, big.mark = ",")),
    data.frame(Characteristic = "Within-cluster pairs, both-empty at Primary",
               Value = format(n_pairs_both_empty, big.mark = ","))
  )

  out <- do.call(rbind, rows)
  attr(out, "footnote") <-
    "Florida = transport-application cohort. Ghana calibration anchor (169 same-visit replicate pairs / 67 patients; per-patient prevalence 16/97 = 16.5%; pair-level Jaccard~depth slope 0.029 per +10x at Primary tier) is cited from [Paper 1] in the manuscript narrative."
  write_table(out, "table1_cohort_structure")
  out
}

# ============================================================================
# Table 2 - Transport test summary
# ============================================================================
make_table2 <- function() {
  log_msg("\nBuilding Table 2 (transport summary)...")

  t1 <- readRDS(file.path(PATHS$transport_t1, "transport_t1.rds"))
  t2 <- readRDS(file.path(PATHS$transport_t2, "transport_t2.rds"))

  # Hurdle decomposition is optional (post-hoc)
  hurdle_path <- file.path(PATHS$transport_t2, "transport_t2_hurdle.rds")
  hurdle <- if (file.exists(hurdle_path)) readRDS(hurdle_path) else NULL

  rows <- list()

  # T1 - prevalence band
  rows[[length(rows) + 1]] <- data.frame(
    Test = "T1: Florida Primary detection rate vs pre-specified band",
    Estimate = sprintf("%.4f (Wilson 95%% CI %.4f-%.4f)",
                       t1$florida_rate,
                       t1$florida_wilson_ci["lower"],
                       t1$florida_wilson_ci["upper"]),
    Criterion = sprintf("Inside band [%.4f, %.4f] (Ghana anchor 0.165 +/- 10pp)",
                        t1$band_lower, t1$band_upper),
    Decision = t1$declared,
    stringsAsFactors = FALSE
  )

  # T2 - Jaccard~depth slope (primary)
  rows[[length(rows) + 1]] <- data.frame(
    Test = "T2: Pair-level Jaccard ~ depth (per +10x)",
    Estimate = sprintf("%+.4f (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d evaluable pairs)",
                       t2$jaccard_slope$estimate,
                       t2$jaccard_slope$ci_lower,
                       t2$jaccard_slope$ci_upper,
                       fmt_p(t2$jaccard_slope$p_value),
                       t2$jaccard_slope$n),
    Criterion = "Positive slope, same direction as Paper 1 Ghana anchor (0.029 per +10x)",
    Decision = t2$declared,
    stringsAsFactors = FALSE
  )

  # T2 - both-empty complementary slope (if available)
  if (!is.null(t2$both_empty_slope)) {
    rows[[length(rows) + 1]] <- data.frame(
      Test = "T2 complementary: P(both-empty) ~ depth (per +10x)",
      Estimate = sprintf("%+.4f (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d pairs)",
                         t2$both_empty_slope$estimate,
                         t2$both_empty_slope$ci_lower,
                         t2$both_empty_slope$ci_upper,
                         fmt_p(t2$both_empty_slope$p_value),
                         t2$both_empty_slope$n),
      Criterion = "Negative slope (joint-detection probability rises with depth)",
      Decision = if (t2$both_empty_slope$estimate < 0 &&
                     t2$both_empty_slope$p_value < 0.05)
                 "CONFIRMS depth-detection coupling" else "INCONCLUSIVE",
      stringsAsFactors = FALSE
    )
  }

  # T2 - drop-cluster-3 sensitivity (if available)
  if (!is.null(t2$drop_cluster3_slope)) {
    rows[[length(rows) + 1]] <- data.frame(
      Test = "T2 sensitivity: Jaccard ~ depth, excluding cluster 3",
      Estimate = sprintf("%+.4f (HC1 95%% CI %+.4f to %+.4f, p=%s)",
                         t2$drop_cluster3_slope$estimate,
                         t2$drop_cluster3_slope$ci_lower,
                         t2$drop_cluster3_slope$ci_upper,
                         fmt_p(t2$drop_cluster3_slope$p_value)),
      Criterion = "Stable sign and similar magnitude to primary T2 result",
      Decision = if (t2$drop_cluster3_slope$estimate > 0 &&
                     t2$drop_cluster3_slope$p_value < 0.05)
                 "STABLE" else "ATTENUATED",
      stringsAsFactors = FALSE
    )
  }

  # Hurdle decomposition (if available)
  if (!is.null(hurdle)) {
    ext <- hurdle$extensive_lpm
    int <- hurdle$intensive_ols
    logit_or <- hurdle$extensive_logit

    rows[[length(rows) + 1]] <- data.frame(
      Test = "Hurdle, extensive margin: P(Jaccard > 0 | union > 0) ~ depth (per +10x, LPM)",
      Estimate = sprintf("%+.4f (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d) | logit OR %.3f (CI %.3f-%.3f)",
                         ext$estimate, ext$ci_lower, ext$ci_upper,
                         fmt_p(ext$p_value), ext$n,
                         logit_or$or, logit_or$or_lower, logit_or$or_upper),
      Criterion = "Positive (any-sharing margin lifts with depth)",
      Decision = if (ext$estimate > 0 && ext$p_value < 0.05)
                 "COUPLING AT EXTENSIVE MARGIN" else "NULL EXTENSIVE",
      stringsAsFactors = FALSE
    )

    rows[[length(rows) + 1]] <- data.frame(
      Test = "Hurdle, intensive margin: E[Jaccard | Jaccard > 0] ~ depth (per +10x, OLS)",
      Estimate = sprintf("%+.4f (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d intersection-positive)",
                         int$estimate, int$ci_lower, int$ci_upper,
                         fmt_p(int$p_value), int$n),
      Criterion = "Directionally consistent with primary T2",
      Decision = if (int$p_value < 0.05) "SIGNIFICANT INTENSIVE COUPLING"
                 else "INTENSIVE NOT DETECTABLE",
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  attr(out, "footnote") <-
    "Decisions report whether each test's pre-specified or expected criterion is met. T1 and T2 are pre-specified (Paper 2 commit). The both-empty complementary diagnostic and the hurdle decomposition are post-hoc and labeled as such in the manuscript."
  write_table(out, "table2_transport_summary")
  out
}

# ============================================================================
# Run all
# ============================================================================
t1 <- make_table1()
t2 <- make_table2()

log_msg("\nAll Paper 2a tables written to %s/", TABLES_DIR)
log_msg("13a_make_paper2a_tables.R complete.\n")
