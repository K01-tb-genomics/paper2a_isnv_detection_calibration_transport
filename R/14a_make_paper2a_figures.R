# =============================================================================
# 14a_make_paper2a_figures.R
#
# Builds Paper 2a (filter transportability) manuscript figures.
#
# Figures produced:
#   Figure 1  Transport diagnostic - two-panel:
#             A. Pair-level Jaccard ~ depth (primary T2)
#             B. P(both-empty at Primary) ~ depth (complementary T2 diagnostic)
#   Figure 2  Hurdle decomposition of T2 (post-hoc) - two-panel:
#             A. Extensive margin: P(Jaccard > 0 | union > 0) ~ depth
#             B. Intensive margin: E[Jaccard | Jaccard > 0] ~ depth
#
# Outputs:
#   outputs/paper2a/figures/figure1_transport.{png,pdf}
#   outputs/paper2a/figures/figure2_hurdle.{png,pdf}
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required.")

log_msg("\n--- 14a_make_paper2a_figures.R ---")

FIGURES_DIR <- PATHS$paper2a_figures
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

library(ggplot2)
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
if (!HAS_PATCHWORK)
  log_msg("  NOTE: patchwork not available; panels saved separately")

# Common theme
fig_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    legend.position     = "bottom",
    plot.title          = element_text(face = "bold", size = 12),
    plot.subtitle       = element_text(size = 10, color = "gray30"),
    plot.caption        = element_text(size = 8, color = "gray40", hjust = 0)
  )

save_fig <- function(p, name, width = 8, height = 9) {
  for (ext in c("png", "pdf")) {
    path <- file.path(FIGURES_DIR, sprintf("%s.%s", name, ext))
    ggsave(path, p, width = width, height = height, dpi = 300, units = "in")
    log_msg("  Wrote %s", path)
  }
}

# HC1 cluster-robust prediction ribbon helper (same form used in
# 15_make_supplementary_figures.R)
cr_ribbon <- function(fit, vcv, depth_grid) {
  X <- cbind(`(Intercept)` = 1, min_depth_10x = depth_grid / 10)
  fit_mu <- as.numeric(X %*% coef(fit))
  fit_se <- sqrt(rowSums((X %*% vcv) * X))
  data.frame(min_depth = depth_grid,
             fit = fit_mu,
             lwr = fit_mu - 1.96 * fit_se,
             upr = fit_mu + 1.96 * fit_se)
}

# Quantile-binned proportion estimator with Wilson CI (used for both-empty
# and extensive-margin binned panels — binary/indicator outcomes).
binned_proportions <- function(pairs, indicator_col, n_bins = 12) {
  brks <- unique(quantile(pairs$min_depth,
                          probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE))
  pairs$depth_bin <- cut(pairs$min_depth, breaks = brks, include.lowest = TRUE)
  do.call(rbind, lapply(split(pairs, pairs$depth_bin), function(d) {
    n <- nrow(d); p <- mean(d[[indicator_col]]); z <- 1.96
    denom  <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denom
    half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
    data.frame(depth_mid = median(d$min_depth),
               p = p, n = n,
               lwr = max(0, center - half),
               upr = min(1, center + half))
  }))
}

# Quantile-binned mean estimator with normal-approx 95% CI on the bin mean.
# Used for the Jaccard panels (1A and 2B) where the outcome is continuous
# in [0, 1]. Boxplots collapse on the zero-inflated 1A distribution and
# overstate the within-bin spread for the 206-pair 2B distribution;
# per-bin mean +/- CI bars match the 1B/2A pattern.
binned_means <- function(pairs, value_col, n_bins = 5) {
  brks <- unique(quantile(pairs$min_depth,
                          probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE))
  pairs$depth_bin <- cut(pairs$min_depth, breaks = brks, include.lowest = TRUE)
  do.call(rbind, lapply(split(pairs, pairs$depth_bin), function(d) {
    n  <- nrow(d)
    m  <- mean(d[[value_col]])
    se <- if (n > 1) sd(d[[value_col]]) / sqrt(n) else NA_real_
    data.frame(depth_mid = median(d$min_depth),
               mean = m, n = n,
               lwr = if (is.na(se)) NA_real_ else max(0, m - 1.96 * se),
               upr = if (is.na(se)) NA_real_ else min(1, m + 1.96 * se))
  }))
}

# ============================================================================
# Figure 1 - Transport diagnostic (Jaccard + both-empty)
# ============================================================================
make_figure1 <- function() {
  log_msg("\nBuilding Figure 1 (transport diagnostic)...")
  
  pairs_path <- file.path(PATHS$transport_t2, "transport_t2_pairs.rds")
  res_path   <- file.path(PATHS$transport_t2, "transport_t2.rds")
  if (!file.exists(pairs_path) || !file.exists(res_path)) {
    log_msg("  WARNING: required RDS not found; skipping Figure 1")
    return(invisible(NULL))
  }
  
  pairs  <- readRDS(pairs_path)
  result <- readRDS(res_path)
  
  # ---- Panel A: Pair-level Jaccard ~ depth (primary T2) -------------------
  # Scatter replaced with per-bin mean Jaccard + 95% CI bars
  # (parallel to the binned-proportion pattern in 1B and 2A). Boxplots
  # were tried first but collapsed because >95% of pairs have Jaccard = 0
  # in every bin, leaving degenerate boxes. Trend line + cluster-robust
  # ribbon overlaid using the saved LM fit/vcov, clipped to binned range.
  evaluable <- pairs[!is.na(pairs$jaccard), ]
  
  jac <- result$jaccard_slope
  ghana_anchor <- result$ghana_anchor
  
  sub_a <- sprintf(
    "Florida T2: %+.4f per +10x (cluster-robust 95%% CI %+.4f to %+.4f, p=%s, n=%d evaluable pairs)",
    jac$estimate, jac$ci_lower, jac$ci_upper, fmt_p(jac$p_value), jac$n
  )
  
  # Per-bin mean Jaccard with normal-approx 95% CI (5 quintiles)
  jac_bins <- binned_means(evaluable, "jaccard", n_bins = 5)
  x_min_a <- min(jac_bins$depth_mid, na.rm = TRUE)
  x_max_a <- max(jac_bins$depth_mid, na.rm = TRUE)
  # Tight y-axis: Jaccard means are O(10^-3) at this dilution
  y_max_a <- max(c(jac_bins$upr, jac_bins$mean), na.rm = TRUE) * 1.20
  
  # Trend line + ribbon from locked LM fit, clipped to binned range
  if (!is.null(jac$fit) && !is.null(jac$vcov)) {
    jac_grid <- seq(x_min_a, x_max_a, length.out = 200)
    jac_rib  <- cr_ribbon(jac$fit, jac$vcov, jac_grid)
    y_max_a  <- max(y_max_a, max(jac_rib$upr, na.rm = TRUE) * 1.10)
  } else {
    jac_rib <- NULL
  }
  
  p_a <- ggplot() +
    geom_pointrange(data = jac_bins,
                    aes(x = depth_mid, y = mean, ymin = lwr, ymax = upr),
                    color = "gray25", size = 0.35, linewidth = 0.5) +
    {if (!is.null(jac_rib)) geom_ribbon(data = jac_rib,
                                        aes(x = min_depth, ymin = lwr, ymax = upr),
                                        fill = "black", alpha = 0.18)} +
    {if (!is.null(jac_rib)) geom_line(data = jac_rib,
                                      aes(x = min_depth, y = fit),
                                      color = "black", linewidth = 0.8)} +
    coord_cartesian(xlim = c(x_min_a, x_max_a), ylim = c(0, y_max_a)) +
    labs(
      subtitle = sub_a,
      x = "Pair minimum sample-level median depth (x)",
      y = "Mean Jaccard index at Primary tier",
      caption = sprintf(
        "Quintile bins (n approx %d/bin); error bars 95%% CI on bin mean. Ghana anchor: %.3f per +10x (n=%d).",
        round(nrow(evaluable) / nrow(jac_bins)),
        ghana_anchor$slope_per_10x, ghana_anchor$n_pairs
      )
    ) +
    fig_theme
  
  # ---- Panel B: P(both-empty at Primary) ~ depth --------------------------
  # Both-empty: pairs where neither sample has any iSNV at Primary tier.
  # Construct indicator if not already in pairs.
  if (!"both_empty" %in% names(pairs)) {
    pairs$both_empty <- as.integer(is.na(pairs$jaccard) &
                                     pairs$n_intersect == 0L &
                                     pairs$n_union == 0L)
  }
  pairs$both_empty_int <- as.integer(pairs$both_empty)
  
  be <- result$both_zero_slope
  sub_b <- if (!is.null(be)) {
    sprintf(
      "Both-empty slope: %+.4f per +10x (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d pairs)",
      be$estimate, be$ci_lower, be$ci_upper, fmt_p(be$p_value), be$n
    )
  } else {
    "Both-empty slope: [not in transport_t2.rds]"
  }
  
  be_bins <- binned_proportions(pairs, "both_empty_int", n_bins = 12)
  
  # Trend line + cluster-robust ribbon, CLIPPED to the binned-data range.
  # Uses the saved LM fit (with min_depth_10x predictor) and HC1 vcov from
  # transport_t2.rds so the line matches the locked headline slope exactly
  # (rather than re-fitting via geom_smooth on a different x-scale).
  x_min_b <- min(be_bins$depth_mid, na.rm = TRUE)
  x_max_b <- max(be_bins$depth_mid, na.rm = TRUE)
  if (!is.null(be$fit) && !is.null(be$vcov)) {
    be_grid <- seq(x_min_b, x_max_b, length.out = 200)
    be_rib  <- cr_ribbon(be$fit, be$vcov, be_grid)
  } else {
    # Fallback: fit on the fly if upstream save did not include $fit/$vcov
    pairs$min_depth_10x <- pairs$min_depth / 10
    fit_be  <- lm(both_empty_int ~ min_depth_10x, data = pairs)
    vcv_be  <- sandwich::vcovCL(fit_be, cluster = pairs$cluster_id,
                                type = "HC1")
    be_grid <- seq(x_min_b, x_max_b, length.out = 200)
    be_rib  <- cr_ribbon(fit_be, vcv_be, be_grid)
  }
  
  p_b <- ggplot() +
    geom_pointrange(data = be_bins,
                    aes(x = depth_mid, y = p, ymin = lwr, ymax = upr),
                    color = "gray25", size = 0.35, linewidth = 0.5) +
    geom_ribbon(data = be_rib,
                aes(x = min_depth, ymin = lwr, ymax = upr),
                fill = "black", alpha = 0.18) +
    geom_line(data = be_rib,
              aes(x = min_depth, y = fit),
              color = "black", linewidth = 0.8) +
    scale_y_continuous(
      labels = function(x) sprintf("%.0f%%", 100 * x)
    ) +
    coord_cartesian(xlim = c(x_min_b, x_max_b), ylim = c(0, 1)) +
    labs(
      subtitle = sub_b,
      x = "Pair minimum sample-level median depth (x)",
      y = "P(both samples empty at Primary)"
    ) +
    fig_theme
  
  # ---- Compose ------------------------------------------------------------
  if (HAS_PATCHWORK) {
    library(patchwork)
    combined <- (p_a / p_b) +
      plot_annotation(
        title = "Figure 1. Transport diagnostic: pair-level Jaccard concordance vs depth, with complementary both-empty panel",
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      ) +
      plot_layout(heights = c(1, 1)) &
      theme(plot.tag = element_text(face = "bold"))
    combined <- combined + plot_annotation(tag_levels = "A")
    save_fig(combined, "figure1_transport", width = 10, height = 10)
  } else {
    save_fig(p_a + labs(title = "Figure 1A. Pair-level Jaccard ~ depth"),
             "figure1A_jaccard", width = 10, height = 5)
    save_fig(p_b + labs(title = "Figure 1B. P(both-empty) ~ depth"),
             "figure1B_both_empty", width = 10, height = 5)
  }
}

# ============================================================================
# Figure 2 - Hurdle decomposition (post-hoc; promoted from supplementary)
# ============================================================================
make_figure2 <- function() {
  log_msg("\nBuilding Figure 2 (hurdle decomposition)...")
  
  pairs_path  <- file.path(PATHS$transport_t2, "transport_t2_pairs.rds")
  hurdle_path <- file.path(PATHS$transport_t2, "transport_t2_hurdle.rds")
  if (!file.exists(pairs_path) || !file.exists(hurdle_path)) {
    log_msg("  WARNING: required RDS not found; skipping Figure 2")
    return(invisible(NULL))
  }
  
  pairs  <- readRDS(pairs_path)
  hurdle <- readRDS(hurdle_path)
  
  ext_lpm   <- hurdle$extensive_lpm
  ext_logit <- hurdle$extensive_logit
  int_ols   <- hurdle$intensive_ols
  
  ext_pairs <- pairs[!is.na(pairs$jaccard), ]
  ext_pairs$intersect_positive <- as.integer(ext_pairs$n_intersect > 0L)
  
  # ---- Panel A: Extensive margin (binned proportions + HC1 ribbon) -------
  ext_bins <- binned_proportions(ext_pairs, "intersect_positive", n_bins = 12)
  # CLIP extensive-margin trend line + ribbon to the binned-data range,
  # not to max(ext_pairs$min_depth). The high-depth tail is sparse and
  # extrapolating into it overstates the range over which the fit is
  # supported by binned data.
  depth_grid_ext <- seq(min(ext_bins$depth_mid, na.rm = TRUE),
                        max(ext_bins$depth_mid, na.rm = TRUE),
                        length.out = 200)
  ext_rib  <- cr_ribbon(ext_lpm$fit, ext_lpm$vcov, depth_grid_ext)
  
  sub_a <- sprintf(
    paste0("Extensive: slope %+.4f per +10x (HC1 95%% CI %+.4f to %+.4f, p%s, n=%d) | ",
           "logit OR %.3f (CI %.3f-%.3f)"),
    ext_lpm$estimate, ext_lpm$ci_lower, ext_lpm$ci_upper,
    if (ext_lpm$p_value < 0.001) "<0.001" else sprintf("=%.3f", ext_lpm$p_value),
    ext_lpm$n,
    ext_logit$or, ext_logit$or_lower, ext_logit$or_upper
  )
  y_max_top <- max(c(ext_bins$upr, ext_rib$upr), na.rm = TRUE) * 1.10
  
  p_a <- ggplot() +
    geom_pointrange(data = ext_bins,
                    aes(x = depth_mid, y = p, ymin = lwr, ymax = upr),
                    color = "gray25", size = 0.35, linewidth = 0.5) +
    geom_ribbon(data = ext_rib,
                aes(x = min_depth, ymin = lwr, ymax = upr),
                fill = "black", alpha = 0.18) +
    geom_line(data = ext_rib,
              aes(x = min_depth, y = fit),
              color = "black", linewidth = 0.8) +
    scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
    coord_cartesian(xlim = range(depth_grid_ext),
                    ylim = c(0, y_max_top)) +
    labs(subtitle = sub_a,
         x = "Pair minimum sample-level median depth (x)",
         y = "P(Jaccard > 0 | union > 0)") +
    fig_theme
  
  # ---- Panel B: Intensive margin (scatter + HC1 ribbon) ------------------
  int_pairs <- ext_pairs[ext_pairs$intersect_positive == 1L, ]
  int_pairs$cluster_id <- factor(as.character(int_pairs$cluster_id))
  
  # CLIP intensive-margin trend line + ribbon to the intersection-positive
  # data range. The intensive fit was estimated on n = 206 pairs whose
  # depth range is much narrower than the full evaluable set; using the
  # full ext_pairs grid would extrapolate the fit far beyond its support
  # and push the ribbon outside the [0,1] axis.
  depth_grid_int <- seq(min(int_pairs$min_depth), max(int_pairs$min_depth),
                        length.out = 200)
  int_rib <- cr_ribbon(int_ols$fit, int_ols$vcov, depth_grid_int)
  
  # Bin int_pairs by depth quartile for per-bin mean +/- CI
  # display (parallel to 1A). Boxplots were tried first but produced
  # boxes covering nearly the full [0,1] y-range in 3 of 4 bins because
  # the conditional Jaccard distribution among intersection-positive
  # pairs is wide; the boxplot conveyed within-bin variance honestly but
  # obscured the trend. Per-bin mean +/- 95% CI matches the 1B/2A grammar.
  int_bins <- binned_means(int_pairs, "jaccard", n_bins = 4)
  
  sub_b <- sprintf(
    "Intensive: slope %+.4f per +10x (HC1 95%% CI %+.4f to %+.4f, p=%s, n=%d intersection-positive)",
    int_ols$estimate, int_ols$ci_lower, int_ols$ci_upper,
    fmt_p(int_ols$p_value), int_ols$n
  )
  
  # Y-axis bounded to data + ribbon range with a small headroom margin
  y_max_b <- max(c(int_bins$upr, int_bins$mean, int_rib$upr),
                 na.rm = TRUE) * 1.10
  y_max_b <- min(y_max_b, 1)
  
  p_b <- ggplot() +
    geom_pointrange(data = int_bins,
                    aes(x = depth_mid, y = mean, ymin = lwr, ymax = upr),
                    color = "gray25", size = 0.35, linewidth = 0.5) +
    geom_ribbon(data = int_rib,
                aes(x = min_depth, ymin = lwr, ymax = upr),
                fill = "black", alpha = 0.18) +
    geom_line(data = int_rib,
              aes(x = min_depth, y = fit),
              color = "black", linewidth = 0.8) +
    # coord_cartesian clips the visible area without dropping out-of-range
    # ribbon/line data, so the cluster-robust ribbon survives near the
    # boundary instead of being silently filtered by scale_y_continuous.
    coord_cartesian(xlim = range(depth_grid_int), ylim = c(0, y_max_b)) +
    labs(subtitle = sub_b,
         x = "Pair minimum sample-level median depth (x)",
         y = "Mean Jaccard (intersection-positive pairs only)",
         caption = sprintf(
           "Quartile bins (approx %d intersection-positive pairs/bin); error bars 95%% CI on bin mean. Cluster 5 contributed no intersection-positive pairs.",
           round(nrow(int_pairs) / nrow(int_bins))
         )) +
    fig_theme
  
  # ---- Compose ------------------------------------------------------------
  if (HAS_PATCHWORK) {
    library(patchwork)
    combined <- (p_a / p_b) +
      plot_annotation(
        title = "Figure 2. Hurdle decomposition of the T2 Jaccard slope (post-hoc)",
        theme = theme(plot.title = element_text(face = "bold", size = 12))
      ) +
      plot_annotation(tag_levels = "A")
    save_fig(combined, "figure2_hurdle", width = 10, height = 9)
  } else {
    save_fig(p_a + labs(title = "Figure 2A. Extensive margin"),
             "figure2A_extensive", width = 10, height = 5)
    save_fig(p_b + labs(title = "Figure 2B. Intensive margin"),
             "figure2B_intensive", width = 10, height = 5)
  }
}

# Local fmt_p (defined here so figure script doesn't depend on the table
# script being sourced first)
fmt_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ============================================================================
# Run all
# ============================================================================
make_figure1()
make_figure2()

log_msg("\nAll Paper 2a figures written to %s/", FIGURES_DIR)
log_msg("14a_make_paper2a_figures.R complete.\n")