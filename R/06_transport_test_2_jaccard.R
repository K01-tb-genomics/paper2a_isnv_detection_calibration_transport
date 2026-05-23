# =============================================================================
# 06_transport_test_2_jaccard.R
#
# Transport Test 2 (within-cluster pair-level Jaccard ~ depth).
#
# Pre-specified criterion (LOCKED at pre-spec commit):
#   In Florida, construct all within-cluster sample pairs (across the 5
#   transmission clusters). For each pair, compute (a) the Jaccard index
#   of variant-site sets at the Primary tier and (b) the per-pair minimum
#   sample-level median read depth, min(depth_med_A, depth_med_B).
#   Regress pair-level Jaccard on min-depth (per +10x increment) with
#   cluster-robust HC1 standard errors on cluster_id (G = 5).
#   Pre-specified criterion: positive slope, in the same direction as
#   Paper 1's Ghana Primary-tier slope (0.029 per +10x, p = 0.05; n = 6
#   evaluable pairs).
#
# Complementary diagnostic (Paper 1 parity, not part of formal criterion):
#   P(both samples silent at Primary) ~ min-depth, linear probability
#   model with the same cluster-robust SE.
#
# Inputs:
#   data_derived/02_ladder_applied/florida_ladder.rds
#   data_derived/02_ladder_applied/florida_variants_at_primary.rds
#
# Outputs:
#   data_derived/06_transport_t2/transport_t2.rds
#   data_derived/06_transport_t2/transport_t2_pairs.rds  (per-pair table)
#   data_derived/06_transport_t2/transport_t2_summary.txt
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

if (!requireNamespace("sandwich", quietly = TRUE))
  stop("Package 'sandwich' is required.")

log_msg("\n--- 06_transport_test_2_jaccard.R ---")

# ---- Load ------------------------------------------------------------------
fl_ladder    <- readRDS(file.path(PATHS$ladder, "florida_ladder.rds"))
fl_v_primary <- readRDS(file.path(PATHS$ladder,
                                  "florida_variants_at_primary.rds"))

require_cols(fl_ladder,    c("Sample", "cluster_id", "depth_med"))
# Variants at Primary may legitimately be empty for some samples.

# ---- Build per-sample variant-site sets at Primary -------------------------
fl_v_primary$site <- variant_site_key(fl_v_primary)
sample_sites <- split(fl_v_primary$site, fl_v_primary$Sample)
# Ensure every Florida sample has an entry, even if empty
for (s in fl_ladder$Sample) {
  if (is.null(sample_sites[[s]])) sample_sites[[s]] <- character(0)
}

# Per-sample median depth lookup
sample_depth   <- setNames(fl_ladder$depth_med, fl_ladder$Sample)
sample_cluster <- setNames(as.character(fl_ladder$cluster_id), fl_ladder$Sample)

# ---- Construct within-cluster pairs ----------------------------------------
build_pairs <- function(ladder) {
  pairs <- list()
  for (cid in unique(ladder$cluster_id)) {
    samps <- ladder$Sample[ladder$cluster_id == cid]
    if (length(samps) < 2) next
    combos <- t(combn(samps, 2))
    pairs[[as.character(cid)]] <- data.frame(
      cluster_id = as.character(cid),
      sample_A   = combos[, 1],
      sample_B   = combos[, 2],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pairs)
}

pairs_df <- build_pairs(fl_ladder)
log_msg("Within-cluster pairs constructed: %d pairs across %d clusters",
        nrow(pairs_df), length(unique(pairs_df$cluster_id)))
log_msg("Pairs per cluster:")
for (cid in unique(pairs_df$cluster_id)) {
  log_msg("  %s: %d pairs", cid, sum(pairs_df$cluster_id == cid))
}

# ---- Compute Jaccard, min-depth, both-zero indicator -----------------------
n <- nrow(pairs_df)
pairs_df$n_A          <- NA_integer_
pairs_df$n_B          <- NA_integer_
pairs_df$n_intersect  <- NA_integer_
pairs_df$n_union      <- NA_integer_
pairs_df$jaccard      <- NA_real_
pairs_df$depth_A      <- NA_real_
pairs_df$depth_B      <- NA_real_
pairs_df$min_depth    <- NA_real_
pairs_df$both_empty   <- FALSE

for (i in seq_len(n)) {
  s_a <- pairs_df$sample_A[i]
  s_b <- pairs_df$sample_B[i]
  set_a <- sample_sites[[s_a]]
  set_b <- sample_sites[[s_b]]
  pairs_df$n_A[i]         <- length(set_a)
  pairs_df$n_B[i]         <- length(set_b)
  pairs_df$n_intersect[i] <- length(intersect(set_a, set_b))
  pairs_df$n_union[i]     <- length(union(set_a, set_b))
  pairs_df$jaccard[i]     <- jaccard(set_a, set_b)
  pairs_df$depth_A[i]     <- sample_depth[s_a]
  pairs_df$depth_B[i]     <- sample_depth[s_b]
  pairs_df$min_depth[i]   <- min(sample_depth[s_a], sample_depth[s_b])
  pairs_df$both_empty[i]  <- (length(set_a) == 0L) && (length(set_b) == 0L)
}

n_total      <- nrow(pairs_df)
n_both_empty <- sum(pairs_df$both_empty)
n_evaluable  <- sum(!is.na(pairs_df$jaccard))

log_msg("\nPair characteristics:")
log_msg("  Total pairs:                   %d", n_total)
log_msg("  Both samples empty at Primary: %d (%.1f%%)",
        n_both_empty, 100 * n_both_empty / n_total)
log_msg("  Evaluable for Jaccard slope:   %d", n_evaluable)
log_msg("  Min-depth range: %.1f - %.1f, median %.1f",
        min(pairs_df$min_depth), max(pairs_df$min_depth),
        median(pairs_df$min_depth))

# ---- Jaccard slope regression (HEADLINE TRANSPORT CRITERION) ---------------
reg_data            <- pairs_df[!is.na(pairs_df$jaccard), , drop = FALSE]
reg_data$min_depth_10x <- reg_data$min_depth / 10

run_slope <- function(formula, data, label) {
  if (nrow(data) < 3L) {
    log_msg("  [%s] Insufficient pairs (n = %d) for slope regression.",
            label, nrow(data))
    return(list(n = nrow(data), estimate = NA_real_, se = NA_real_,
                ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
                fit = NULL, vcov = NULL))
  }
  fit <- lm(formula, data = data)
  vcv <- tryCatch(
    sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1"),
    error = function(e) NULL
  )
  if (is.null(vcv)) {
    log_msg("  [%s] vcovCL failed; using model-based SE.", label)
    vcv <- vcov(fit)
  }
  beta <- coef(fit)["min_depth_10x"]
  se   <- sqrt(vcv["min_depth_10x", "min_depth_10x"])
  z    <- beta / se
  list(
    n        = nrow(data),
    estimate = beta,
    se       = se,
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    p_value  = 2 * (1 - pnorm(abs(z))),
    fit      = fit,
    vcov     = vcv
  )
}

log_msg("\nJaccard ~ min_depth (per +10x), Florida within-cluster pairs:")
jacc_slope <- run_slope(jaccard ~ min_depth_10x, reg_data, "Jaccard")
if (!is.na(jacc_slope$estimate)) {
  log_msg("  N evaluable pairs: %d", jacc_slope$n)
  log_msg("  Slope: %.4f per +10x (HC1 cluster-robust SE %.4f)",
          jacc_slope$estimate, jacc_slope$se)
  log_msg("  95%% CI: [%.4f, %.4f]",
          jacc_slope$ci_lower, jacc_slope$ci_upper)
  log_msg("  p-value: %.4f", jacc_slope$p_value)
}

sign_positive <- !is.na(jacc_slope$estimate) && jacc_slope$estimate > 0
declared <- if (is.na(jacc_slope$estimate)) {
  "INDETERMINATE (insufficient pairs)"
} else if (sign_positive) {
  "TRANSPORT CRITERION MET (positive slope, same direction as Ghana)"
} else {
  "TRANSPORT CRITERION FAILED (slope not positive)"
}
log_msg("  Decision: %s", declared)

# ---- Both-zero slope (complementary diagnostic, not formal criterion) ------
bz_data <- pairs_df
bz_data$min_depth_10x <- bz_data$min_depth / 10
bz_data$both_empty_int <- as.integer(bz_data$both_empty)

log_msg("\nBoth-zero slope (complementary diagnostic, Paper 1 parity):")
bz_slope <- run_slope(both_empty_int ~ min_depth_10x, bz_data, "both-zero")
if (!is.na(bz_slope$estimate)) {
  log_msg("  N pairs: %d", bz_slope$n)
  log_msg("  Slope: %.4f per +10x (HC1 cluster-robust SE %.4f)",
          bz_slope$estimate, bz_slope$se)
  log_msg("  95%% CI: [%.4f, %.4f]",
          bz_slope$ci_lower, bz_slope$ci_upper)
  log_msg("  p-value: %.4f", bz_slope$p_value)
  log_msg("  Note: linear probability model on both-empty indicator.")
}

# ---- Both-zero LOGISTIC refit (scale-robustness check, reviewer A3) --------
# Added 2026-05-23 in response to reviewer comment: report both LPM and
# logistic for the binary both-empty outcome. LPM remains the headline
# (same probability scale as Jaccard slope); logistic OR reported alongside.

run_logit_slope <- function(formula, data, label) {
  if (nrow(data) < 3L) {
    log_msg("  [%s] Insufficient pairs (n = %d) for logistic regression.",
            label, nrow(data))
    return(list(n = nrow(data), beta = NA_real_, se = NA_real_,
                or = NA_real_, or_lower = NA_real_, or_upper = NA_real_,
                p_value = NA_real_, fit = NULL, vcov = NULL))
  }
  fit <- glm(formula, family = binomial(link = "logit"), data = data)
  vcv <- tryCatch(
    sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1"),
    error = function(e) NULL
  )
  if (is.null(vcv)) {
    log_msg("  [%s] vcovCL failed; using model-based SE.", label)
    vcv <- vcov(fit)
  }
  beta <- coef(fit)["min_depth_10x"]
  se   <- sqrt(vcv["min_depth_10x", "min_depth_10x"])
  z    <- beta / se
  list(
    n        = nrow(data),
    beta     = beta,
    se       = se,
    or       = exp(beta),
    or_lower = exp(beta - 1.96 * se),
    or_upper = exp(beta + 1.96 * se),
    p_value  = 2 * (1 - pnorm(abs(z))),
    fit      = fit,
    vcov     = vcv
  )
}

log_msg("\nBoth-zero logistic (scale-robustness check):")
bz_logit <- run_logit_slope(both_empty_int ~ min_depth_10x, bz_data,
                            "both-zero-logit")
if (!is.na(bz_logit$or)) {
  log_msg("  N pairs: %d", bz_logit$n)
  log_msg("  Coefficient (log-OR per +10x): %.4f (cluster-robust SE %.4f)",
          bz_logit$beta, bz_logit$se)
  log_msg("  OR per +10x: %.3f (95%% CI %.3f to %.3f)",
          bz_logit$or, bz_logit$or_lower, bz_logit$or_upper)
  log_msg("  p-value: %.4f", bz_logit$p_value)
}

# ---- Anchor reference (for substantive comparison only) --------------------
log_msg("\nGhana Primary anchor (Paper 1 reference): slope = 0.029 per +10x, p = 0.05, n = 6 evaluable pairs.")

# ---- Save ------------------------------------------------------------------
result <- list(
  test               = "Transport Test 2 (within-cluster Jaccard ~ depth in Florida)",
  pair_summary       = list(
    total_pairs      = n_total,
    both_empty_pairs = n_both_empty,
    evaluable_pairs  = n_evaluable,
    pairs_by_cluster = table(pairs_df$cluster_id)
  ),
  jaccard_slope      = jacc_slope,
  both_zero_slope    = bz_slope,
  both_zero_logit    = bz_logit,
  ghana_anchor       = list(slope_per_10x = 0.029, p_value = 0.05,
                            n_pairs = 6),
  pre_spec_criterion = "positive slope on Jaccard ~ min_depth (per +10x), same direction as Paper 1 Ghana Primary",
  sign_positive      = sign_positive,
  declared           = declared
)

save_with_provenance(
  result, file.path(PATHS$transport_t2, "transport_t2.rds"),
  inputs = list(
    florida_ladder      = file.path(PATHS$ladder, "florida_ladder.rds"),
    florida_v_primary   = file.path(PATHS$ladder,
                                    "florida_variants_at_primary.rds")
  ),
  note = sprintf("Transport Test 2: Jaccard slope = %s per +10x (n = %d pairs); %s",
                 if (is.na(jacc_slope$estimate)) "NA"
                 else sprintf("%.4f", jacc_slope$estimate),
                 jacc_slope$n, declared)
)

# Per-pair table (for downstream tables/figures)
saveRDS(pairs_df,
        file.path(PATHS$transport_t2, "transport_t2_pairs.rds"))

# Human-readable summary
fmt <- function(x, d = 4)
  if (is.na(x)) "NA" else sprintf(sprintf("%%.%df", d), x)

summary_lines <- c(
  sprintf("Transport Test 2 summary  (%s)", format(Sys.time())),
  "",
  "Florida within-cluster pair construction:",
  sprintf("  Total pairs:                   %d", n_total),
  sprintf("  Both samples empty at Primary: %d (%.1f%%)",
          n_both_empty, 100 * n_both_empty / n_total),
  sprintf("  Evaluable for Jaccard slope:   %d", n_evaluable),
  "",
  "HEADLINE: Jaccard ~ min_depth (per +10x), HC1 cluster-robust on cluster_id:",
  sprintf("  Slope:    %s per +10x", fmt(jacc_slope$estimate)),
  sprintf("  SE:       %s", fmt(jacc_slope$se)),
  sprintf("  95%% CI:   [%s, %s]",
          fmt(jacc_slope$ci_lower), fmt(jacc_slope$ci_upper)),
  sprintf("  p-value:  %s", fmt(jacc_slope$p_value)),
  sprintf("  Decision: %s", declared),
  "",
  "Complementary diagnostic - P(both empty) ~ min_depth:",
  "  Linear probability model (headline scale):",
  sprintf("    Slope:    %s per +10x", fmt(bz_slope$estimate)),
  sprintf("    95%% CI:   [%s, %s]",
          fmt(bz_slope$ci_lower), fmt(bz_slope$ci_upper)),
  sprintf("    p-value:  %s", fmt(bz_slope$p_value)),
  "  Logistic regression (scale-robustness check):",
  sprintf("    OR:       %s per +10x", fmt(bz_logit$or, 3)),
  sprintf("    95%% CI:   [%s, %s]",
          fmt(bz_logit$or_lower, 3), fmt(bz_logit$or_upper, 3)),
  sprintf("    p-value:  %s", fmt(bz_logit$p_value)),
  "",
  "Reference: Paper 1 Ghana Primary slope = 0.029 per +10x, p = 0.05, n = 6 evaluable pairs."
)
writeLines(summary_lines,
           file.path(PATHS$transport_t2, "transport_t2_summary.txt"))

# ---- Tidy CSV: LPM vs logistic side-by-side for both-empty diagnostic ------
# Output: data_derived/06_transport_t2/transport_t2_both_empty_comparison.csv
# Consumed by 13a_make_paper2a_tables.R / 17_make_eappendix_tables.R if you
# want to surface this in an eTable.
bz_comparison <- data.frame(
  model      = c("Linear probability model", "Logistic regression"),
  scale      = c("slope per +10x (Pr units)", "OR per +10x"),
  n_pairs    = c(bz_slope$n, bz_logit$n),
  estimate   = c(bz_slope$estimate, bz_logit$or),
  ci_lower   = c(bz_slope$ci_lower, bz_logit$or_lower),
  ci_upper   = c(bz_slope$ci_upper, bz_logit$or_upper),
  p_value    = c(bz_slope$p_value, bz_logit$p_value),
  stringsAsFactors = FALSE
)
write.csv(bz_comparison,
          file.path(PATHS$transport_t2,
                    "transport_t2_both_empty_comparison.csv"),
          row.names = FALSE)
log_msg("Wrote both-empty LPM/logistic comparison CSV.")

# Drop-cluster-3 sensitivity for T2 — confirms cluster 3 doesn't drive the result
fl_no3 <- fl_ladder[fl_ladder$cluster_id != 3, ]
# Then re-source 06 with this restricted data, OR pull the pairs_df from
# transport_t2_pairs.rds and refit on subset:
pairs_no3 <- pairs_df[pairs_df$cluster_id != "3" & !is.na(pairs_df$jaccard), ]
pairs_no3$min_depth_10x <- pairs_no3$min_depth / 10
fit_no3 <- lm(jaccard ~ min_depth_10x, data = pairs_no3)
sandwich::vcovCL(fit_no3, cluster = pairs_no3$cluster_id, type = "HC1") |>
  diag() |> sqrt() -> se_no3
cat("No-cluster-3 slope:", coef(fit_no3)["min_depth_10x"],
    "SE:", se_no3["min_depth_10x"], "\n")

log_msg("06_transport_test_2_jaccard.R complete.\n")