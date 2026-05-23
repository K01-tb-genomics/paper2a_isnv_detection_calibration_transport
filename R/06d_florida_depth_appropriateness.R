# =============================================================================
# 06d_florida_depth_appropriateness.R
#
# Auxiliary script (eAppendix worked example): derive a defensible Florida
# depth floor for iSNV-based work from the transport-test outputs. Two
# parallel thresholds are reported so the reader can see whether the
# empirical and modeled answers agree.
#
# Threshold rules (between-host appropriate; Ghana within-host anchor is
# explicitly NOT used here because expected pair-level any-overlap is much
# lower for between-host within-cluster pairs than for within-host replicate
# pairs):
#
#   1. Empirical Q5: the 5th percentile of min(depth_A, depth_B) among
#      intersection-positive pairs. Reading: floor at or above which 95%
#      of empirically observed any-overlap pairs sit.
#
#   2. Modeled 5%: depth at which the extensive-margin logistic fit
#      predicts P(Jaccard > 0 | union > 0) = 0.05. Reading: smoothed
#      version of the empirical answer; uses the fitted curve rather than
#      the raw quantile.
#
# Outputs:
#   outputs/paper2a/figures/eFigure_florida_depth_appropriateness.{png,pdf}
#   outputs/paper2a/tables/eTable_florida_depth_appropriateness.csv
#   data_derived/06_transport_t2/florida_depth_thresholds.rds
#     (consumed by Paper 2b's R/01_build_analytic_frames.R if frame
#      selection is updated to use the derived threshold)
#
# Framing:
#   This script is illustrative, not corrective. The Paper 2a main text
#   keeps the diagnostic-not-corrective stance; this analysis is a worked
#   example of one defensible procedure for converting the diagnostic into
#   a deployment threshold.
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")  # provides log_msg(), require_cols(), etc.

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required.")
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

log_msg("\n--- 06d_florida_depth_appropriateness.R ---")

# ---- Load -------------------------------------------------------------------
pairs_path  <- file.path(PATHS$transport_t2, "transport_t2_pairs.rds")
hurdle_path <- file.path(PATHS$transport_t2, "transport_t2_hurdle.rds")
ladder_path <- file.path(PATHS$ladder,       "florida_ladder.rds")

for (p in c(pairs_path, hurdle_path, ladder_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

pairs     <- readRDS(pairs_path)
hurdle    <- readRDS(hurdle_path)
fl_ladder <- readRDS(ladder_path)

ext_logit <- hurdle$extensive_logit
if (is.null(ext_logit$fit))
  stop("ext_logit$fit not present in transport_t2_hurdle.rds; re-run R/06b.")

# ---- Empirical Q5 of intersection-positive pairs ---------------------------
int_pos <- pairs[!is.na(pairs$jaccard) & pairs$jaccard > 0, ]
log_msg("Intersection-positive pairs: %d", nrow(int_pos))
empirical_q5  <- as.numeric(quantile(int_pos$min_depth, 0.05))
empirical_q10 <- as.numeric(quantile(int_pos$min_depth, 0.10))
log_msg("Empirical Q5 of min_depth (intersection-positive): %.1fx",
        empirical_q5)
log_msg("Empirical Q10 (reference):                          %.1fx",
        empirical_q10)

# ---- Modeled 5% threshold from extensive-margin logistic fit ---------------
# Solve: logit(0.05) = b0 + b1 * (depth / 10) for depth
beta <- coef(ext_logit$fit)
if (!all(c("(Intercept)", "min_depth_10x") %in% names(beta)))
  stop("Unexpected coefficient names in ext_logit$fit: ",
       paste(names(beta), collapse = ", "))

depth_for_p <- function(p_target) {
  10 * (qlogis(p_target) - beta["(Intercept)"]) / beta["min_depth_10x"]
}
modeled_p05 <- as.numeric(depth_for_p(0.05))
log_msg("Modeled depth at P(any overlap) = 0.05: %.1fx", modeled_p05)

# ---- Sample and pair retention at each threshold ---------------------------
n_fl <- nrow(fl_ladder)
retain_at <- function(thresh) {
  surviving_samples <- fl_ladder$Sample[fl_ladder$depth_med >= thresh]
  surviving_pairs   <- pairs[pairs$sample_A %in% surviving_samples &
                               pairs$sample_B %in% surviving_samples, ]
  data.frame(
    threshold              = thresh,
    n_specimens_retained   = length(surviving_samples),
    pct_specimens_retained = 100 * length(surviving_samples) / n_fl,
    n_pairs_retained       = nrow(surviving_pairs),
    n_intersect_pos_retained = sum(surviving_pairs$jaccard > 0,
                                   na.rm = TRUE)
  )
}

thresholds_df <- rbind(
  cbind(method = "Current Florida floor (50x)",      retain_at(50)),
  cbind(method = "Empirical Q5 (intersect-positive)", retain_at(empirical_q5)),
  cbind(method = "Modeled P(any overlap) = 5%",       retain_at(modeled_p05))
)
log_msg("\nRetention at each candidate threshold:")
print(thresholds_df, row.names = FALSE)

# ---- Save tidy outputs -----------------------------------------------------
out_tables <- PATHS$paper2a_tables
if (!dir.exists(out_tables)) dir.create(out_tables, recursive = TRUE)
write.csv(thresholds_df,
          file.path(out_tables, "eTable_florida_depth_appropriateness.csv"),
          row.names = FALSE)

# Machine-readable thresholds for Paper 2b consumption
thresholds_obj <- list(
  computed_at                    = Sys.time(),
  current_florida_floor          = 50,
  empirical_q5_intersect_pos     = empirical_q5,
  empirical_q10_intersect_pos    = empirical_q10,
  modeled_p05_extensive_margin   = modeled_p05,
  n_intersect_positive_pairs     = nrow(int_pos),
  retention                      = thresholds_df,
  note = paste(
    "Empirical and modeled thresholds derived from the Paper 2a Transport",
    "Test 2 hurdle decomposition. Empirical Q5: depth at or above which 95%",
    "of intersection-positive pairs sit. Modeled 5%: depth at which the",
    "extensive-margin logistic fit predicts P(any overlap) = 0.05. Intended",
    "as one defensible candidate for the Paper 2b Florida analytic frame",
    "depth floor, not a universal rule."
  )
)
saveRDS(thresholds_obj,
        file.path(PATHS$transport_t2, "florida_depth_thresholds.rds"))

# ---- Figure: paired density + Jaccard overlay -----------------------------
library(ggplot2)

# Densities of min_depth for two pair sets
all_pairs_df <- data.frame(
  min_depth = pairs$min_depth[!is.na(pairs$jaccard)],
  pair_type = "All evaluable pairs"
)
int_pos_df <- data.frame(
  min_depth = int_pos$min_depth,
  pair_type = "Intersection-positive pairs"
)
dens_df <- rbind(all_pairs_df, int_pos_df)
dens_df$pair_type <- factor(dens_df$pair_type,
                            levels = c("All evaluable pairs",
                                       "Intersection-positive pairs"))

# Common x-range for both panels (data-supported, no extrapolation)
x_lo <- min(dens_df$min_depth, na.rm = TRUE)
x_hi <- quantile(dens_df$min_depth, 0.995, na.rm = TRUE)

# Threshold lines are inline-labeled directly on each panel (no legend
# entry), so the figure is self-contained and the reader does not need
# the caption to identify dashed vs dotted lines.
emp_label <- sprintf("Empirical Q5: %.0fx", empirical_q5)
mod_label <- sprintf("Modeled 5%%: %.0fx",  modeled_p05)

# Small x-offset for label text so it does not sit directly on the line
label_offset <- (x_hi - x_lo) * 0.012

p_density <- ggplot(dens_df, aes(x = min_depth, fill = pair_type)) +
  geom_density(alpha = 0.55, color = NA) +
  geom_vline(xintercept = empirical_q5, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_vline(xintercept = modeled_p05,  linetype = "dotted",
             color = "black", linewidth = 0.6) +
  annotate("text", x = empirical_q5 + label_offset, y = Inf,
           label = emp_label, hjust = 0, vjust = 1.5,
           size = 3.2, fontface = "plain") +
  annotate("text", x = modeled_p05 + label_offset, y = Inf,
           label = mod_label, hjust = 0, vjust = 1.5,
           size = 3.2, fontface = "plain") +
  scale_fill_manual(values = c("All evaluable pairs"         = "gray70",
                               "Intersection-positive pairs" = "black"),
                    name = NULL) +
  coord_cartesian(xlim = c(x_lo, x_hi)) +
  labs(x = NULL,
       y = "Density",
       subtitle = "Where the any-overlap signal lives in depth space") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p_jaccard <- ggplot(int_pos, aes(x = min_depth, y = jaccard)) +
  geom_jitter(width = 1.5, height = 0, alpha = 0.45, size = 1.3,
              color = "gray25") +
  geom_vline(xintercept = empirical_q5, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_vline(xintercept = modeled_p05,  linetype = "dotted",
             color = "black", linewidth = 0.6) +
  annotate("text", x = empirical_q5 + label_offset, y = Inf,
           label = emp_label, hjust = 0, vjust = 1.5,
           size = 3.2, fontface = "plain") +
  annotate("text", x = modeled_p05 + label_offset, y = Inf,
           label = mod_label, hjust = 0, vjust = 1.5,
           size = 3.2, fontface = "plain") +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0, 1)) +
  labs(x = "Pair minimum sample-level median depth (x)",
       y = "Jaccard index",
       subtitle = sprintf(
         "Jaccard among intersection-positive pairs (n = %d)",
         nrow(int_pos))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# Compose
out_figs <- PATHS$paper2a_figures
if (!dir.exists(out_figs)) dir.create(out_figs, recursive = TRUE)

if (HAS_PATCHWORK) {
  library(patchwork)
  # No bottom caption: the inline labels and the eAppendix prose carry
  # all the figure-defining information, so dropping the caption keeps
  # the figure self-contained and prevents text overflow.
  combined <- (p_density / p_jaccard) +
    plot_annotation(
      title = "Florida depth appropriateness: empirical and modeled floors",
      theme = theme(plot.title = element_text(face = "bold", size = 12)),
      tag_levels = "A"
    )
  ggsave(file.path(out_figs, "eFigure_florida_depth_appropriateness.png"),
         combined, width = 8, height = 8, dpi = 300, units = "in")
  ggsave(file.path(out_figs, "eFigure_florida_depth_appropriateness.pdf"),
         combined, width = 8, height = 8, units = "in")
} else {
  ggsave(file.path(out_figs, "eFigure_florida_depth_appropriateness_A.png"),
         p_density, width = 8, height = 4, dpi = 300, units = "in")
  ggsave(file.path(out_figs, "eFigure_florida_depth_appropriateness_B.png"),
         p_jaccard, width = 8, height = 4, dpi = 300, units = "in")
}

log_msg("\nSaved eFigure_florida_depth_appropriateness.{png,pdf}")
log_msg("Saved eTable_florida_depth_appropriateness.csv")
log_msg("Saved florida_depth_thresholds.rds (for Paper 2b consumption)")
log_msg("06d_florida_depth_appropriateness.R complete.\n")