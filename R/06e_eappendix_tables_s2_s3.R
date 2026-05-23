# =============================================================================
# 06e_eappendix_tables_s2_s3.R
#
# Produces two supplementary tables for the Paper 2a eAppendix:
#
#   eTable S2: Both-empty proportion by depth decile (Florida within-cluster
#              pairs). Wilson 95% CI per bin.
#
#   eTable S3: Leave-one-cluster-out Jaccard slope refits across all 5
#              Florida transmission clusters. Cluster-robust HC1 SEs.
#
# Self-contained: defines its own binning and Wilson CI logic so it does
# not depend on the figure script.
#
# Inputs:
#   data_derived/06_transport_t2/transport_t2_pairs.rds
#
# Outputs:
#   outputs/paper2a/tables/eTable_S2_both_empty_by_depth_decile.csv
#   outputs/paper2a/tables/eTable_S3_leave_one_cluster_out.csv
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

if (!requireNamespace("sandwich", quietly = TRUE))
  stop("Package 'sandwich' is required.")

log_msg("\n--- 06e_eappendix_tables_s2_s3.R ---")

# ---- Load ------------------------------------------------------------------
pairs_path <- file.path(PATHS$transport_t2, "transport_t2_pairs.rds")
if (!file.exists(pairs_path))
  stop("Required input not found: ", pairs_path,
       "\nRun R/06_transport_test_2_jaccard.R first.")

pairs <- readRDS(pairs_path)

# Ensure required columns
require_cols(pairs, c("cluster_id", "min_depth", "both_empty", "jaccard"))

# ============================================================================
# eTable S2: Both-empty proportion by depth decile, Wilson 95% CI
# ============================================================================
log_msg("\nBuilding eTable S2 (both-empty by depth decile)...")

# Decile breaks on min_depth across all 13,076 pairs
brks_s2 <- unique(quantile(pairs$min_depth, probs = seq(0, 1, by = 0.1),
                           na.rm = TRUE))
pairs$decile_bin <- cut(pairs$min_depth, breaks = brks_s2,
                        include.lowest = TRUE)

wilson_ci <- function(p, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, center - half), min(1, center + half))
}

bins <- split(pairs, pairs$decile_bin)
s2_rows <- lapply(seq_along(bins), function(i) {
  d <- bins[[i]]
  n <- nrow(d)
  p <- mean(d$both_empty)
  ci <- wilson_ci(p, n)
  data.frame(
    decile                  = i,
    depth_lower             = min(d$min_depth),
    depth_upper             = max(d$min_depth),
    median_min_depth        = median(d$min_depth),
    n_pairs                 = n,
    n_both_empty            = sum(d$both_empty),
    p_both_empty            = p,
    wilson_95_ci_lower      = ci[1],
    wilson_95_ci_upper      = ci[2],
    stringsAsFactors        = FALSE
  )
})
s2 <- do.call(rbind, s2_rows)

# Display formatting helpers
fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
log_msg("  Decile | median depth | n | P(both-empty) | 95%% CI")
for (i in seq_len(nrow(s2))) {
  log_msg("  %2d     | %5.1fx       | %4d | %s        | (%s, %s)",
          s2$decile[i], s2$median_min_depth[i], s2$n_pairs[i],
          fmt_pct(s2$p_both_empty[i]),
          fmt_pct(s2$wilson_95_ci_lower[i]),
          fmt_pct(s2$wilson_95_ci_upper[i]))
}

out_tables <- PATHS$paper2a_tables
if (!dir.exists(out_tables)) dir.create(out_tables, recursive = TRUE)
write.csv(s2,
          file.path(out_tables, "eTable_S2_both_empty_by_depth_decile.csv"),
          row.names = FALSE)
log_msg("  Wrote eTable_S2_both_empty_by_depth_decile.csv")

# ============================================================================
# eTable S3: Leave-one-cluster-out Jaccard slope refits
# ============================================================================
log_msg("\nBuilding eTable S3 (leave-one-cluster-out refits)...")

clusters <- sort(unique(as.character(pairs$cluster_id)))
log_msg("  Clusters in Florida frame: %s",
        paste(clusters, collapse = ", "))

fit_jaccard_slope <- function(sub_pairs) {
  eval_sub <- sub_pairs[!is.na(sub_pairs$jaccard), ]
  if (nrow(eval_sub) < 3L) {
    return(list(n_pairs = nrow(sub_pairs),
                n_evaluable = nrow(eval_sub),
                beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  eval_sub$min_depth_10x <- eval_sub$min_depth / 10
  fit <- lm(jaccard ~ min_depth_10x, data = eval_sub)
  vcv <- tryCatch(
    sandwich::vcovCL(fit, cluster = eval_sub$cluster_id, type = "HC1"),
    error = function(e) vcov(fit)
  )
  beta <- coef(fit)["min_depth_10x"]
  se   <- sqrt(vcv["min_depth_10x", "min_depth_10x"])
  z    <- beta / se
  list(n_pairs    = nrow(sub_pairs),
       n_evaluable = nrow(eval_sub),
       beta       = beta,
       se         = se,
       p          = 2 * (1 - pnorm(abs(z))))
}

# Headline (all clusters) for reference
all_fit <- fit_jaccard_slope(pairs)

# Leave-one-out refits
s3_rows <- lapply(clusters, function(cid) {
  sub  <- pairs[as.character(pairs$cluster_id) != cid, ]
  res  <- fit_jaccard_slope(sub)
  data.frame(
    excluded_cluster        = cid,
    n_pairs_retained        = res$n_pairs,
    n_evaluable_retained    = res$n_evaluable,
    slope_per_10x           = res$beta,
    cluster_robust_se       = res$se,
    ci_lower                = res$beta - 1.96 * res$se,
    ci_upper                = res$beta + 1.96 * res$se,
    p_value                 = res$p,
    direction               = if (is.na(res$beta)) "NA"
    else if (res$beta > 0) "Positive"
    else "Negative",
    stringsAsFactors        = FALSE
  )
})
s3 <- do.call(rbind, s3_rows)

# Add headline row at the top for context
headline_row <- data.frame(
  excluded_cluster        = "(none, all 5 clusters)",
  n_pairs_retained        = all_fit$n_pairs,
  n_evaluable_retained    = all_fit$n_evaluable,
  slope_per_10x           = all_fit$beta,
  cluster_robust_se       = all_fit$se,
  ci_lower                = all_fit$beta - 1.96 * all_fit$se,
  ci_upper                = all_fit$beta + 1.96 * all_fit$se,
  p_value                 = all_fit$p,
  direction               = if (all_fit$beta > 0) "Positive (headline)" else "Negative (headline)",
  stringsAsFactors        = FALSE
)
s3 <- rbind(headline_row, s3)

log_msg("  Excluded | n pairs | n evaluable | slope per +10x | SE | direction")
for (i in seq_len(nrow(s3))) {
  log_msg("  %-25s | %6d | %6d | %+.4f | %.4f | %s",
          s3$excluded_cluster[i], s3$n_pairs_retained[i],
          s3$n_evaluable_retained[i], s3$slope_per_10x[i],
          s3$cluster_robust_se[i], s3$direction[i])
}

write.csv(s3,
          file.path(out_tables, "eTable_S3_leave_one_cluster_out.csv"),
          row.names = FALSE)
log_msg("  Wrote eTable_S3_leave_one_cluster_out.csv")

log_msg("\n06e_eappendix_tables_s2_s3.R complete.\n")