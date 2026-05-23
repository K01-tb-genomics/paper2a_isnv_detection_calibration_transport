# =============================================================================
# isnv_helpers.R
#
# Shared utilities used across Paper 2 pipeline scripts. No side effects.
# Sourced after 00_pipeline_config.R.
# =============================================================================

# ---- PE/PPE flag extraction ------------------------------------------------
# Mirrors Paper 1's helper: tolerates logical, numeric, character, factor
# encodings of the PE/PPE indicator. Returns 0/1 integer vector.
get_ppe_flag <- function(v) {
  if ("ppe_flag" %in% names(v)) {
    pp <- v$ppe_flag
    if (is.logical(pp))                             return(as.integer(pp))
    if (is.numeric(pp))                             return(as.integer(pp != 0))
    if (is.character(pp) || is.factor(pp))
      return(as.integer(as.character(pp) %in%
                        c("PPE", "pe/ppe", "PE/PPE", "TRUE", "true", "1")))
  }
  if ("PPE" %in% names(v)) {
    pp <- v$PPE
    if (is.logical(pp)) return(as.integer(pp))
    if (is.numeric(pp)) return(as.integer(pp != 0))
    if (is.character(pp))
      return(as.integer(pp %in% c("1", "TRUE", "true", "PPE", "pe/ppe")))
  }
  rep(0L, nrow(v))  # If no PE/PPE column at all, treat all as non-PE/PPE.
}

# ---- Variant-level filter application --------------------------------------
# Returns logical vector indexing rows of `variants` that PASS the given tier
# threshold. PE/PPE variants are excluded prior to threshold application
# (matches Paper 1 calibration convention).
filter_variants_by_tier <- function(variants, thr) {
  required <- c("DP", "AD1", "MAF", "FILTER")
  missing  <- setdiff(required, names(variants))
  if (length(missing) > 0)
    stop("filter_variants_by_tier: missing columns: ",
         paste(missing, collapse = ", "))

  ppe <- get_ppe_flag(variants)
  pass_filter   <- variants$FILTER == "PASS"
  pass_dp       <- variants$DP    >= thr$DP_min
  pass_ad1      <- variants$AD1   >= thr$AD1_min
  pass_maf_low  <- variants$MAF   >= thr$MAF_min
  pass_maf_high <- variants$MAF   <= thr$MAF_max
  not_ppe       <- ppe == 0L

  pass_filter & pass_dp & pass_ad1 & pass_maf_low & pass_maf_high & not_ppe
}

# ---- Per-sample iSNV count at a given tier ---------------------------------
# Returns a data.frame: Sample, n_iSNV. Samples with no variants (after
# filter) appear with n_iSNV = 0; samples present in `sample_ids` but absent
# from `variants` are zero-filled.
count_isnv_per_sample <- function(variants, thr, sample_ids) {
  pass <- filter_variants_by_tier(variants, thr)
  v_pass <- variants[pass, , drop = FALSE]
  cnts <- as.data.frame(table(Sample = v_pass$Sample),
                        stringsAsFactors = FALSE)
  names(cnts) <- c("Sample", "n_iSNV")
  cnts$n_iSNV <- as.integer(cnts$n_iSNV)
  out <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
  out <- merge(out, cnts, by = "Sample", all.x = TRUE, sort = FALSE)
  out$n_iSNV[is.na(out$n_iSNV)] <- 0L
  out
}

# ---- Wilson 95% CI for a binomial proportion -------------------------------
wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0L) return(c(p = NA_real_, lower = NA_real_, upper = NA_real_))
  alpha <- 1 - conf
  z <- qnorm(1 - alpha / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- p + z^2 / (2 * n)
  half   <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)
  c(p = p, lower = (centre - half) / denom, upper = (centre + half) / denom)
}

# ---- Jaccard similarity for two sets ---------------------------------------
# Returns NA if the union is empty (both sets silent).
jaccard <- function(set_a, set_b) {
  if (length(set_a) == 0 && length(set_b) == 0) return(NA_real_)
  inter <- length(intersect(set_a, set_b))
  un    <- length(union(set_a, set_b))
  if (un == 0) return(NA_real_)
  inter / un
}

# ---- Variant-site key for Jaccard computations -----------------------------
# Constructs a stable per-variant identifier (pos:ref:alt). Used to compute
# Jaccard between two specimens' iSNV call sets.
variant_site_key <- function(v) {
  if (nrow(v) == 0L) return(character(0))
  candidates <- list(
    chrom = c("chrom", "CHROM", "chromosome", "Chromosome", "chr", "CHR"),
    pos   = c("pos",   "POS",   "position",   "Position"),
    ref   = c("ref",   "REF"),
    alt   = c("alt",   "ALT")
  )
  resolved <- vapply(candidates, function(cands) {
    hit <- intersect(cands, names(v))
    if (length(hit) == 0L) NA_character_ else hit[1]
  }, character(1))
  if (any(is.na(resolved))) {
    stop(sprintf(
      "variant_site_key: cannot resolve columns: %s. Available: %s",
      paste(names(resolved)[is.na(resolved)], collapse = ", "),
      paste(names(v), collapse = ", ")
    ))
  }
  paste(v[[resolved["chrom"]]], v[[resolved["pos"]]],
        v[[resolved["ref"]]],   v[[resolved["alt"]]],
        sep = "_")
}

# ---- Required-column validation --------------------------------------------
# Hard-fails with a useful message if a data frame is missing any required
# column. Used at script entry points.
require_cols <- function(df, required, df_label = "data") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing required column(s): %s",
                 df_label, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

# ---- Lightweight logger ----------------------------------------------------
log_msg <- function(...) {
  args <- list(...)
  s <- if (length(args) == 1) as.character(args[[1]]) else do.call(sprintf, args)
  message(s)
  invisible(s)
}

# ---- Save with provenance --------------------------------------------------
# Wraps saveRDS with a small TXT sibling logging input file paths and rowcounts.
# Useful for tracing which inputs produced a given derived artifact.
save_with_provenance <- function(obj, path, inputs = list(), note = NULL) {
  saveRDS(obj, path)
  prov_path <- sub("\\.rds$", "_PROVENANCE.txt", path)
  lines <- c(
    sprintf("Saved: %s", path),
    sprintf("Time:  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Inputs:"
  )
  if (length(inputs) > 0) {
    for (nm in names(inputs)) {
      lines <- c(lines, sprintf("  %-20s = %s", nm, inputs[[nm]]))
    }
  } else {
    lines <- c(lines, "  (none specified)")
  }
  if (!is.null(note)) lines <- c(lines, "Note:", paste0("  ", note))
  writeLines(lines, prov_path)
  invisible(path)
}
