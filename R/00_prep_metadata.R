# =============================================================================
# 00_prep_metadata.R
#
# Spec-independent metadata harmonization. Runs once per pipeline invocation.
# Reads raw cohort metadata for Florida and Ghana, harmonizes column types,
# applies binary lineage encoding, and writes prepped metadata to
# data_derived/00_metadata/.
#
# Inputs (read-only):
#   data_raw/florida_epiData.rds  Sample-level Florida metadata
#   data_raw/gh_meta.rds    Sample-level Ghana metadata (all visits)
#
# Outputs:
#   data_derived/00_metadata/florida_metadata.rds  Cleaned, type-stabilized
#   data_derived/00_metadata/ghana_metadata.rds    Cleaned, all visits retained
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

log_msg("\n--- 00_prep_metadata.R ---")

# ---- Load raw metadata ------------------------------------------------------
fl_meta_raw <- readRDS(file.path(PATHS$data_raw, "fl_meta.rds"))
gh_meta_raw <- readRDS(file.path(PATHS$data_raw, "gh_meta.rds"))

log_msg("Florida raw: %d samples, %d columns", nrow(fl_meta_raw), ncol(fl_meta_raw))
log_msg("Ghana raw:   %d samples, %d columns", nrow(gh_meta_raw), ncol(gh_meta_raw))

# ---- Required columns -------------------------------------------------------
fl_required <- c("Sample", "cluster_id", "HIV", "lineage",
                 "Age", "Sex", "smear", "run_id",
                 "depth_med", "callable", "Year")
gh_required <- c("Sample", "patientId", "visit", "HIV", "lineage",
                 "Age", "Sex", "depth_med", "callable_frac_20x")

require_cols(fl_meta_raw, fl_required, "florida_epiData.rds")
require_cols(gh_meta_raw, gh_required, "gh_meta.rds")

# ---- Harmonize types --------------------------------------------------------
# Lineage is binary L4 vs non_L4 in both cohorts (per pre-spec). Coerce to a
# factor with reference level "non_L4" so the regression coefficient
# represents the L4 effect.
harmonize_lineage <- function(x) {
  x_chr <- as.character(x)
  out   <- ifelse(x_chr %in% c("L4", "lineage4", "4"), "L4", "L_non_L4")
  factor(out, levels = c("L_non_L4", "L4"))
}

# HIV: 0/1 integer (canonical). Tolerate "Positive"/"Negative" or TRUE/FALSE.
harmonize_hiv <- function(x) {
  if (is.logical(x))  return(as.integer(x))
  if (is.numeric(x))  return(as.integer(x != 0))
  if (is.character(x) || is.factor(x)) {
    x_chr <- as.character(x)
    return(as.integer(x_chr %in% c("1", "Positive", "POSITIVE", "Pos","HIV+",
                                   "TRUE", "true", "Yes", "Y")))
  }
  stop("Unrecognized HIV encoding")
}

# Sex: factor "F"/"M" with M as reference (typical convention).
harmonize_sex <- function(x) {
  x_chr <- as.character(x)
  out   <- ifelse(x_chr %in% c("F", "Female", "female", "0"), "F",
                  ifelse(x_chr %in% c("M", "Male", "male", "1"), "M", NA))
  factor(out, levels = c("M", "F"))
}

# Smear (Florida only): factor with negative as reference.
harmonize_smear <- function(x) {
  if (is.numeric(x))  return(factor(ifelse(x > 0, "pos", "neg"),
                                    levels = c("neg", "pos", "unknown")))
  if (is.logical(x))  return(factor(ifelse(x, "pos", "neg"),
                                    levels = c("neg", "pos", "unknown")))
  x_chr <- as.character(x)
  out <- ifelse(x_chr %in% c("pos", "Positive", "POS", "Pos", "+", "1",
                             "TRUE", "true"), "pos",
                ifelse(x_chr %in% c("neg", "Negative", "NEG", "Neg", "-", "0",
                                    "FALSE", "false"), "neg",
                       ifelse(x_chr %in% c("Unknown", "UNKNOWN", "unknown", "U", "?",
                                           "NA", "N/A"), "unknown", NA)))
  factor(out, levels = c("neg", "pos", "unknown"))
}

# ---- Florida prep -----------------------------------------------------------
fl_meta <- fl_meta_raw
fl_meta$lineage   <- harmonize_lineage(fl_meta$lineage)
fl_meta$HIV       <- harmonize_hiv(fl_meta$HIV)
fl_meta$Sex       <- harmonize_sex(fl_meta$Sex)
fl_meta$smear     <- harmonize_smear(fl_meta$smear)
fl_meta$cluster_id <- as.factor(fl_meta$cluster_id)
fl_meta$run_id    <- as.factor(fl_meta$run_id)
fl_meta$Age       <- as.numeric(fl_meta$Age)
fl_meta$depth_med <- as.numeric(fl_meta$depth_med)
fl_meta$callable  <- as.numeric(fl_meta$callable)
fl_meta$Year      <- as.integer(fl_meta$Year)

# ---- Ghana prep -------------------------------------------------------------
gh_meta <- gh_meta_raw
gh_meta$lineage   <- harmonize_lineage(gh_meta$lineage)
gh_meta$HIV       <- harmonize_hiv(gh_meta$HIV)
gh_meta$Sex       <- harmonize_sex(gh_meta$Sex)
gh_meta$Age       <- as.numeric(gh_meta$Age)
gh_meta$visit     <- factor(as.character(gh_meta$visit),
                            levels = c("M0", "M1", "M2"))
gh_meta$depth_med <- as.numeric(gh_meta$depth_med)
gh_meta$callable  <- as.numeric(gh_meta$callable)

# ---- Quick diagnostic summary -----------------------------------------------
log_msg("Florida prepped: %d samples, %d unique cluster_id, %d unique run_id, year range %d-%d",
        nrow(fl_meta), length(unique(fl_meta$cluster_id)),
        length(unique(fl_meta$run_id)),
        min(fl_meta$Year, na.rm = TRUE), max(fl_meta$Year, na.rm = TRUE))
log_msg("Ghana prepped:   %d samples (%d M0, %d M1, %d M2) from %d patients",
        nrow(gh_meta),
        sum(gh_meta$visit == "M0", na.rm = TRUE),
        sum(gh_meta$visit == "M1", na.rm = TRUE),
        sum(gh_meta$visit == "M2", na.rm = TRUE),
        length(unique(gh_meta$patientId)))

# ---- Save -------------------------------------------------------------------
save_with_provenance(
  fl_meta, file.path(PATHS$meta, "florida_metadata.rds"),
  inputs = list(florida_metadata_raw = file.path(PATHS$data_raw,
                                                 "florida_metadata.rds")),
  note = "Type-harmonized; lineage binarized to L4 vs non_L4."
)
save_with_provenance(
  gh_meta, file.path(PATHS$meta, "ghana_metadata.rds"),
  inputs = list(ghana_metadata_raw = file.path(PATHS$data_raw,
                                               "ghana_metadata.rds")),
  note = "Type-harmonized; lineage binarized; all visits retained (M0 filter applied in 01_build_analytic_frames.R)."
)

log_msg("00_prep_metadata.R complete.\n")
