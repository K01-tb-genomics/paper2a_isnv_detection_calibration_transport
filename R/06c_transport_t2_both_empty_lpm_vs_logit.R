# =============================================================================
# 06c_transport_t2_both_empty_lpm_vs_logit.R
#
# Auxiliary script (Hitchings A3 response): produce a small comparison figure
# overlaying the LPM and logistic fits for the both-empty diagnostic on the
# binned proportion of both-empty pairs versus pair-minimum depth.
#
# Purpose:
#   - Visually confirm that the LPM-based slope (headline) and the logistic OR
#     (scale-robustness check) tell the same substantive story across the
#     observed depth range.
#   - Output goes to the eAppendix, not to Figure 1 (Figure 1B remains LPM-only
#     to keep main-text figures simple).
#
# Inputs:
#   data_derived/06_transport_t2/transport_t2_pairs.rds  (pair-level frame,
#       built by R/06_transport_test_2_jaccard.R)
#   data_derived/06_transport_t2/transport_t2.rds        (model fits)
#
# Outputs:
#   outputs/paper2a/figures/eFigure_S_both_empty_lpm_vs_logit.{png,pdf}
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required.")
if (!requireNamespace("sandwich", quietly = TRUE))
  stop("Package 'sandwich' is required.")

cat("--- 06c_transport_t2_both_empty_lpm_vs_logit.R [v4 ordered] ---\n")

# ---- Load -------------------------------------------------------------------
pairs_df <- readRDS(file.path(PATHS$transport_t2, "transport_t2_pairs.rds"))
t2_res   <- readRDS(file.path(PATHS$transport_t2, "transport_t2.rds"))

bz_data <- pairs_df
bz_data$min_depth_10x  <- bz_data$min_depth / 10
bz_data$both_empty_int <- as.integer(bz_data$both_empty)
cat(sprintf("  bz_data: %d rows, min_depth range [%g, %g], both_empty rate %.3f\n",
            nrow(bz_data),
            min(bz_data$min_depth, na.rm = TRUE),
            max(bz_data$min_depth, na.rm = TRUE),
            mean(bz_data$both_empty_int, na.rm = TRUE)))

# ---- Refit models (so this script is self-contained) -----------------------
fit_lpm <- lm(both_empty_int ~ min_depth_10x, data = bz_data)
fit_log <- glm(both_empty_int ~ min_depth_10x,
               family = binomial(link = "logit"), data = bz_data)
cat("  LPM and logit fits complete\n")

# ---- Bin observed data by depth decile for display -------------------------
# Build per-bin summary explicitly via split() + lapply() to avoid the
# aggregate(formula, FUN = function-returning-named-vector) fragility.
# IMPORTANT: this block must run BEFORE the prediction grid below, which
# references bin_summary$mid.
brks <- quantile(bz_data$min_depth, probs = seq(0, 1, by = 0.1),
                 na.rm = TRUE)
brks <- unique(brks)
if (length(brks) < 2L)
  stop("Not enough distinct min_depth values to construct depth bins.")
bz_data$bin <- cut(bz_data$min_depth, breaks = brks,
                   include.lowest = TRUE)
cat(sprintf("  %d depth breaks, %d bin levels, %d non-NA bin rows\n",
            length(brks),
            length(levels(bz_data$bin)),
            sum(!is.na(bz_data$bin))))

.bins_nonna   <- bz_data[!is.na(bz_data$bin), c("min_depth", "both_empty_int", "bin")]
.split_data   <- split(.bins_nonna, .bins_nonna$bin)
.bin_rows     <- lapply(names(.split_data), function(b) {
  d <- .split_data[[b]]
  n <- nrow(d)
  if (n == 0L) return(NULL)
  p <- mean(d$both_empty_int, na.rm = TRUE)
  data.frame(
    bin = b,
    n   = n,
    p   = p,
    mid = median(d$min_depth, na.rm = TRUE),
    se  = sqrt(p * (1 - p) / n),
    stringsAsFactors = FALSE
  )
})
bin_summary <- do.call(rbind, .bin_rows)
if (is.null(bin_summary) || nrow(bin_summary) == 0L)
  stop("bin_summary is empty after split-aggregate; check upstream pairs_df.")
bin_summary$lo <- pmax(bin_summary$p - 1.96 * bin_summary$se, 0)
bin_summary$hi <- pmin(bin_summary$p + 1.96 * bin_summary$se, 1)
cat(sprintf("  bin_summary built: %d rows\n", nrow(bin_summary)))

# ---- Predicted-probability grid CLIPPED to binned data range ---------------
# Bound the prediction lines to the highest binned midpoint instead of max
# observed depth (the high-depth tail is sparse and neither model is well
# identified there). This block depends on bin_summary above, so order
# matters.
x_min <- min(bin_summary$mid, na.rm = TRUE)
x_max <- max(bin_summary$mid, na.rm = TRUE)
depth_grid <- seq(x_min, x_max, length.out = 200)
grid_df <- data.frame(
  min_depth     = depth_grid,
  min_depth_10x = depth_grid / 10
)
grid_df$p_lpm <- predict(fit_lpm, newdata = grid_df)
grid_df$p_lpm <- pmax(pmin(grid_df$p_lpm, 1), 0)  # clip for display only
grid_df$p_log <- predict(fit_log, newdata = grid_df, type = "response")
cat(sprintf("  Prediction grid built over [%g, %g]\n", x_min, x_max))

# ---- Plot -------------------------------------------------------------------
library(ggplot2)

p_fig <- ggplot() +
  geom_point(data = bin_summary,
             aes(x = mid, y = p),
             color = "black", size = 2) +
  geom_errorbar(data = bin_summary,
                aes(x = mid, ymin = lo, ymax = hi),
                width = 0, color = "black") +
  geom_line(data = grid_df,
            aes(x = min_depth, y = p_lpm, linetype = "Linear probability model"),
            color = "black", linewidth = 0.7) +
  geom_line(data = grid_df,
            aes(x = min_depth, y = p_log, linetype = "Logistic regression"),
            color = "black", linewidth = 0.7) +
  scale_linetype_manual(name = "Model",
                        values = c("Linear probability model" = "solid",
                                   "Logistic regression"      = "dashed")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     name = "P(both samples empty at Primary)") +
  scale_x_continuous(name = "Pair minimum sample-level median depth") +
  coord_cartesian(xlim = c(x_min, x_max)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title    = element_blank())

# ---- Save -------------------------------------------------------------------
# Paper 2a manuscript-output convention (see R/14a_make_paper2a_figures.R):
# figures go to PATHS$paper2a_figures = outputs/paper2a/figures.
FIGURES_DIR <- PATHS$paper2a_figures
if (!dir.exists(FIGURES_DIR)) dir.create(FIGURES_DIR, recursive = TRUE)

ggsave(file.path(FIGURES_DIR, "eFigure_S_both_empty_lpm_vs_logit.png"),
       p_fig, width = 6, height = 4, dpi = 300, units = "in")
ggsave(file.path(FIGURES_DIR, "eFigure_S_both_empty_lpm_vs_logit.pdf"),
       p_fig, width = 6, height = 4, units = "in")

cat(sprintf("Saved eFigure_S_both_empty_lpm_vs_logit.{png,pdf} to %s/\n",
            FIGURES_DIR))
cat("06c_transport_t2_both_empty_lpm_vs_logit.R complete.\n")