# =============================================================================
# run_paper2a.R
#
# Master orchestrator for Paper 2a (filter transportability).
#
# Phase 1 (shared infrastructure):
#   00_prep_metadata, 01_build_analytic_frames, 02_apply_filter_ladder
# Phase 2 (Paper 2a primary analyses):
#   05_transport_test_1_band, 06_transport_test_2_jaccard,
#   06b_transport_t2_hurdle
# Phase 3 (Paper 2a eAppendix analyses):
#   06c_transport_t2_both_empty_lpm_vs_logit,
#   06d_florida_depth_appropriateness,
#   06e_eappendix_tables_s2_s3
# Phase 4 (Paper 2a manuscript deliverables):
#   13a_make_paper2a_tables, 14a_make_paper2a_figures
#
# Outputs land in outputs/paper2a/{tables,figures}/.
#
# Halt-on-error policy:
#   - Phase 1 and Phase 2 scripts halt the pipeline on error (results
#     downstream of these depend on their RDS outputs).
#   - Phase 3 (eAppendix) and Phase 4 (manuscript outputs) log-and-continue
#     so a single eAppendix-script failure does not block the main tables
#     and figures.
#   - All script statuses are logged to
#     outputs/orchestrator_logs/run_paper2a_<timestamp>.log.
#
# Usage (from project root):
#   Rscript R/run_paper2a.R
#   OR
#   source("R/run_paper2a.R") in an interactive session
# =============================================================================

# ---- Working directory check ------------------------------------------------
if (!file.exists("R/00_pipeline_config.R")) {
  stop("Run from project root. Expected R/00_pipeline_config.R to exist.")
}

# ---- Bootstrap config and logging -------------------------------------------
source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

LOG_DIR <- PATHS$logs
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(LOG_DIR,
                      sprintf("run_paper2a_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  sink(type = "message"); sink(); close(log_con)
}, add = TRUE)

cat(strrep("=", 78), "\n", sep = "")
cat("Paper 2a orchestrator (filter transportability)\n")
cat(sprintf("Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("R version: %s\n", paste(R.version$major, R.version$minor, sep = ".")))
cat(sprintf("Seed:      %d\n", SEED))
cat(strrep("=", 78), "\n\n", sep = "")
banner()

# ---- Script execution sequence ----------------------------------------------
# Each entry: list(path, phase, halt_on_error).
PAPER2A_SCRIPTS <- list(
  # Phase 1 (shared infrastructure)
  list(path = "R/00_prep_metadata.R",                          phase = 1, halt = TRUE),
  list(path = "R/01_build_analytic_frames.R",                  phase = 1, halt = TRUE),
  list(path = "R/02_apply_filter_ladder.R",                    phase = 1, halt = TRUE),
  # Phase 2 (Paper 2a primary analyses)
  list(path = "R/05_transport_test_1_band.R",                  phase = 2, halt = TRUE),
  list(path = "R/06_transport_test_2_jaccard.R",               phase = 2, halt = TRUE),
  list(path = "R/06b_transport_t2_hurdle.R",                   phase = 2, halt = TRUE),
  # Phase 3 (Paper 2a eAppendix analyses; log-and-continue)
  list(path = "R/06c_transport_t2_both_empty_lpm_vs_logit.R",  phase = 3, halt = FALSE),
  list(path = "R/06d_florida_depth_appropriateness.R",         phase = 3, halt = FALSE),
  list(path = "R/06e_eappendix_tables_s2_s3.R",                phase = 3, halt = FALSE),
  # Phase 4 (Paper 2a manuscript deliverables; log-and-continue)
  list(path = "R/13a_make_paper2a_tables.R",                   phase = 4, halt = FALSE),
  list(path = "R/14a_make_paper2a_figures.R",                  phase = 4, halt = FALSE)
)

run_script <- function(path) {
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat(sprintf("Running: %s\n", path))
  cat(strrep("-", 78), "\n", sep = "")
  t0 <- Sys.time()
  res <- tryCatch({
    source(path, local = FALSE, echo = FALSE)
    list(status = "OK", err = NULL)
  }, error = function(e) {
    list(status = "ERROR", err = conditionMessage(e))
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (res$status == "OK") {
    cat(sprintf("  Completed in %.1fs\n", dt))
  } else {
    cat(sprintf("  FAILED after %.1fs: %s\n", dt, res$err))
  }
  res$elapsed_s <- dt
  res
}

results <- list()
for (entry in PAPER2A_SCRIPTS) {
  r <- run_script(entry$path)
  r$phase <- entry$phase
  r$halt  <- entry$halt
  results[[entry$path]] <- r
  if (r$status == "ERROR" && isTRUE(entry$halt)) {
    cat("\n*** Paper 2a halt-on-error: stopping at ", entry$path, " ***\n", sep = "")
    stop(sprintf("Halting due to error in %s: %s", entry$path, r$err))
  }
}

# ---- Summary ----------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("Paper 2a orchestrator summary\n")
cat(strrep("=", 78), "\n", sep = "")
phase_names <- c("Shared infrastructure",
                 "Primary analyses",
                 "eAppendix analyses",
                 "Manuscript deliverables")
for (ph in 1:4) {
  cat(sprintf("\nPhase %d (%s):\n", ph, phase_names[ph]))
  for (path in names(results)) {
    if (results[[path]]$phase == ph) {
      cat(sprintf("  [%s]  %s  (%.1fs)\n",
                  results[[path]]$status, path, results[[path]]$elapsed_s))
    }
  }
}

n_ok    <- sum(vapply(results, function(x) x$status == "OK",    logical(1)))
n_err   <- sum(vapply(results, function(x) x$status == "ERROR", logical(1)))
n_total <- length(results)
cat(sprintf("\nTotals: %d/%d scripts OK; %d failed.\n", n_ok, n_total, n_err))
cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("Outputs:\n  Tables:  %s\n  Figures: %s\n",
            PATHS$paper2a_tables, PATHS$paper2a_figures))
cat(sprintf("Log:      %s\n", log_file))
cat(strrep("=", 78), "\n", sep = "")

# ---- Expected manuscript artifacts (for human verification) -----------------
cat("\nExpected Paper 2a manuscript artifacts:\n")
cat("  Main tables\n")
cat("    table1_cohort_structure.{csv,txt}\n")
cat("    table2_transport_summary.{csv,txt}\n")
cat("  Main figures\n")
cat("    figure1_transport.{png,pdf}\n")
cat("    figure2_hurdle.{png,pdf}\n")
cat("  eAppendix tables\n")
cat("    eTable_florida_depth_appropriateness.csv\n")
cat("    eTable_S2_both_empty_by_depth_decile.csv\n")
cat("    eTable_S3_leave_one_cluster_out.csv\n")
cat("  eAppendix figures\n")
cat("    eFigure_S_both_empty_lpm_vs_logit.{png,pdf}\n")
cat("    eFigure_florida_depth_appropriateness.{png,pdf}\n")
cat(strrep("=", 78), "\n", sep = "")

invisible(results)
