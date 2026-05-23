# =============================================================================
# 06b_transport_t2_hurdle.R
#
# POST-HOC DECOMPOSITION of Transport Test 2 Jaccard slope.
# Not part of the pre-specified analysis plan. Surfaced transparently.
#
# Decomposes headline Jaccard slope into:
#   Extensive margin: P(Jaccard > 0 | union > 0) ~ min_depth (per +10x)
#                     fit as both LPM (parity with both-empty slope)
#                     and logit (methodological completeness)
#   Intensive margin: E[Jaccard | Jaccard > 0] ~ min_depth (per +10x)
#                     OLS, HC1 cluster-robust on cluster_id.
#
# Inputs:
#   data_derived/06_transport_t2/transport_t2_pairs.rds
#   data_derived/06_transport_t2/transport_t2.rds   (for cross-reference)
#
# Outputs:
#   data_derived/06_transport_t2/transport_t2_hurdle.rds
#   data_derived/06_transport_t2/transport_t2_hurdle_summary.txt
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

if (!requireNamespace("sandwich", quietly = TRUE))
  stop("Package 'sandwich' is required.")

log_msg("\n--- 06b_transport_t2_hurdle.R (post-hoc decomposition) ---")

pairs_path <- file.path(PATHS$transport_t2, "transport_t2_pairs.rds")
if (!file.exists(pairs_path)) {
  stop("transport_t2_pairs.rds not found; run 06_transport_test_2_jaccard.R first.")
}
pairs   <- readRDS(pairs_path)
parent  <- readRDS(file.path(PATHS$transport_t2, "transport_t2.rds"))

# ---- Helpers ---------------------------------------------------------------
slope_lm <- function(formula, data, label) {
  if (nrow(data) < 3L) {
    log_msg("  [%s] Insufficient pairs (n = %d).", label, nrow(data))
    return(list(n = nrow(data), estimate = NA, se = NA,
                ci_lower = NA, ci_upper = NA, p_value = NA,
                fit = NULL, vcov = NULL))
  }
  fit <- lm(formula, data = data)
  vcv <- tryCatch(
    sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1"),
    error = function(e) { log_msg("  [%s] vcovCL failed; using model-based SE.", label); vcov(fit) }
  )
  beta <- coef(fit)["min_depth_10x"]
  se   <- sqrt(vcv["min_depth_10x", "min_depth_10x"])
  z    <- beta / se
  list(n = nrow(data), estimate = unname(beta), se = unname(se),
       ci_lower = unname(beta - 1.96 * se),
       ci_upper = unname(beta + 1.96 * se),
       p_value  = 2 * (1 - pnorm(abs(z))),
       fit = fit, vcov = vcv)
}

slope_logit <- function(formula, data, label) {
  if (nrow(data) < 3L) {
    return(list(n = nrow(data), estimate = NA, se = NA,
                ci_lower = NA, ci_upper = NA, p_value = NA,
                fit = NULL, vcov = NULL,
                or = NA, or_lower = NA, or_upper = NA))
  }
  fit <- glm(formula, data = data, family = binomial(link = "logit"))
  vcv <- tryCatch(
    sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1"),
    error = function(e) { log_msg("  [%s] vcovCL failed; using model-based SE.", label); vcov(fit) }
  )
  beta <- coef(fit)["min_depth_10x"]
  se   <- sqrt(vcv["min_depth_10x", "min_depth_10x"])
  z    <- beta / se
  list(n = nrow(data), estimate = unname(beta), se = unname(se),
       ci_lower = unname(beta - 1.96 * se),
       ci_upper = unname(beta + 1.96 * se),
       p_value  = 2 * (1 - pnorm(abs(z))),
       or       = exp(unname(beta)),
       or_lower = exp(unname(beta - 1.96 * se)),
       or_upper = exp(unname(beta + 1.96 * se)),
       fit = fit, vcov = vcv)
}

# ---- Construct subsets -----------------------------------------------------
ext_data <- pairs[!is.na(pairs$jaccard), , drop = FALSE]
ext_data$min_depth_10x      <- ext_data$min_depth / 10
ext_data$intersect_positive <- as.integer(ext_data$n_intersect > 0L)

int_data <- ext_data[ext_data$intersect_positive == 1L, , drop = FALSE]

n_union_positive   <- nrow(ext_data)
n_intersect_pos    <- nrow(int_data)
n_intersect_zero   <- n_union_positive - n_intersect_pos

log_msg("\nHurdle decomposition data:")
log_msg("  Union-positive pairs (extensive denominator): %d", n_union_positive)
log_msg("  Intersection-positive pairs (intensive numerator/denominator): %d (%.1f%% of union-positive)",
        n_intersect_pos, 100 * n_intersect_pos / n_union_positive)
log_msg("  Jaccard == 0 pairs (extensive zeros): %d", n_intersect_zero)

# ---- Fit ------------------------------------------------------------------
log_msg("\nExtensive margin (LPM): P(Jaccard > 0 | union > 0) ~ min_depth (per +10x)")
ext_lpm <- slope_lm(intersect_positive ~ min_depth_10x, ext_data, "any-sharing LPM")
if (!is.na(ext_lpm$estimate)) {
  log_msg("  Slope: %+.4f per +10x (HC1 SE %.4f)", ext_lpm$estimate, ext_lpm$se)
  log_msg("  95%% CI: [%+.4f, %+.4f]; p = %.4f",
          ext_lpm$ci_lower, ext_lpm$ci_upper, ext_lpm$p_value)
}

log_msg("\nExtensive margin (logit): same outcome, logistic link")
ext_logit <- slope_logit(intersect_positive ~ min_depth_10x, ext_data, "any-sharing logit")
if (!is.na(ext_logit$estimate)) {
  log_msg("  beta: %+.4f (HC1 SE %.4f); OR per +10x: %.3f (95%% CI %.3f-%.3f); p = %.4f",
          ext_logit$estimate, ext_logit$se,
          ext_logit$or, ext_logit$or_lower, ext_logit$or_upper,
          ext_logit$p_value)
}

log_msg("\nIntensive margin (OLS): E[Jaccard | Jaccard > 0] ~ min_depth (per +10x)")
int_ols <- slope_lm(jaccard ~ min_depth_10x, int_data, "level-of-sharing OLS")
if (!is.na(int_ols$estimate)) {
  log_msg("  Slope: %+.4f per +10x (HC1 SE %.4f)", int_ols$estimate, int_ols$se)
  log_msg("  95%% CI: [%+.4f, %+.4f]; p = %.4f",
          int_ols$ci_lower, int_ols$ci_upper, int_ols$p_value)
}

# ---- Cross-reference parent slopes ----------------------------------------
log_msg("\nFor reference (from pre-spec script 06):")
log_msg("  Headline Jaccard slope (union > 0): %+.4f per +10x (CI [%+.4f, %+.4f])",
        parent$jaccard_slope$estimate,
        parent$jaccard_slope$ci_lower,
        parent$jaccard_slope$ci_upper)
log_msg("  Both-empty slope (all pairs):      %+.4f per +10x (CI [%+.4f, %+.4f])",
        parent$both_zero_slope$estimate,
        parent$both_zero_slope$ci_lower,
        parent$both_zero_slope$ci_upper)

# ---- Save -----------------------------------------------------------------
result <- list(
  test       = "Transport Test 2 hurdle decomposition (POST-HOC)",
  status     = "post-hoc; not part of pre-specified analysis plan",
  pair_counts = list(
    union_positive       = n_union_positive,
    intersect_positive   = n_intersect_pos,
    intersect_zero       = n_intersect_zero
  ),
  extensive_lpm    = ext_lpm,
  extensive_logit  = ext_logit,
  intensive_ols    = int_ols,
  parent_headline  = parent$jaccard_slope,
  parent_both_zero = parent$both_zero_slope
)

save_with_provenance(
  result, file.path(PATHS$transport_t2, "transport_t2_hurdle.rds"),
  inputs = list(
    pairs_rds  = pairs_path,
    parent_rds = file.path(PATHS$transport_t2, "transport_t2.rds")
  ),
  note = sprintf(
    "Post-hoc hurdle: extensive LPM = %+.4f/+10x; intensive OLS = %+.4f/+10x; n_int_pos = %d / %d",
    ext_lpm$estimate, int_ols$estimate, n_intersect_pos, n_union_positive)
)

fmt <- function(x, d = 4)
  if (is.na(x)) "NA" else sprintf(sprintf("%%.%df", d), x)

summary_lines <- c(
  sprintf("Transport T2 hurdle decomposition (POST-HOC)  (%s)",
          format(Sys.time())),
  "",
  "Status: post-hoc, not pre-specified. Decomposes the pre-spec Jaccard",
  "slope into extensive (any sharing) and intensive (level given sharing)",
  "margins on union-positive pairs.",
  "",
  sprintf("Pair counts:"),
  sprintf("  Union-positive (extensive denominator):  %d", n_union_positive),
  sprintf("  Intersection-positive (intensive set):   %d (%.1f%%)",
          n_intersect_pos, 100 * n_intersect_pos / n_union_positive),
  sprintf("  Intersection-zero (extensive zeros):     %d", n_intersect_zero),
  "",
  "Extensive margin (LPM): P(Jaccard > 0 | union > 0) ~ min_depth (per +10x)",
  sprintf("  Slope:    %s per +10x", fmt(ext_lpm$estimate)),
  sprintf("  95%% CI:   [%s, %s]", fmt(ext_lpm$ci_lower), fmt(ext_lpm$ci_upper)),
  sprintf("  p-value:  %s", fmt(ext_lpm$p_value)),
  "",
  "Extensive margin (logit): same outcome, logistic link",
  sprintf("  beta:     %s",       fmt(ext_logit$estimate)),
  sprintf("  OR / +10x: %s (95%% CI %s - %s)",
          fmt(ext_logit$or, 3),
          fmt(ext_logit$or_lower, 3),
          fmt(ext_logit$or_upper, 3)),
  sprintf("  p-value:  %s", fmt(ext_logit$p_value)),
  "",
  "Intensive margin (OLS): E[Jaccard | Jaccard > 0] ~ min_depth (per +10x)",
  sprintf("  Slope:    %s per +10x", fmt(int_ols$estimate)),
  sprintf("  95%% CI:   [%s, %s]", fmt(int_ols$ci_lower), fmt(int_ols$ci_upper)),
  sprintf("  p-value:  %s", fmt(int_ols$p_value)),
  "",
  "For reference (pre-spec, from 06_transport_test_2_jaccard.R):",
  sprintf("  Headline Jaccard:  %s (CI %s - %s)",
          fmt(parent$jaccard_slope$estimate),
          fmt(parent$jaccard_slope$ci_lower),
          fmt(parent$jaccard_slope$ci_upper)),
  sprintf("  Both-empty:        %s (CI %s - %s)",
          fmt(parent$both_zero_slope$estimate),
          fmt(parent$both_zero_slope$ci_lower),
          fmt(parent$both_zero_slope$ci_upper))
)
writeLines(summary_lines,
           file.path(PATHS$transport_t2, "transport_t2_hurdle_summary.txt"))

log_msg("06b_transport_t2_hurdle.R complete.\n")