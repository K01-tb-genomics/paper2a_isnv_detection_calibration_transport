# Paper 2 split: code reorganization

This documents the reorganization of the Paper 2 code into two paper-specific
pipelines, locked 2026-05-17.

## What changed

The original integrated Paper 2 (HIV-iSNV + filter transport + Edwards
operationalization) is split into two manuscripts that share infrastructure
but produce separate outputs:

- **Paper 2a**: filter transportability across sequencing depth regimes.
- **Paper 2b**: Edwards-2023 operationalization for covariate-differential
  outcome misclassification in TB genomic epidemiology, with HIV-iSNV as the
  worked example. 

Submission order: Paper 1 first, Paper 2a second (with preprint posted on
submission), Paper 2b third (citing Paper 2a preprint).

## Paper 2b cohort scope (Florida-only, locked)

Paper 2b is restricted to the Florida cohort. The Ghana M0 cohort retains
its role as the calibration source for the filter ladder (Paper 1) and the
source cohort for the transport diagnostic (Paper 2a), but does not enter
Paper 2b's HIV-iSNV inference.

Empirical justification: the operating-characteristics simulation
(`R/16_operating_characteristics.R`, eFigure S2 in the eAppendix) showed
that the Ghana cohort has 2.2 percent power against the observed |log PR|
pattern at alpha = 0.05. With 24 HIV+ specimens from 12 HIV+ patients and
cluster-robust standard errors on patientId, the cohort is severely
underpowered for HIV-stratified exposure-outcome inference. Reporting
Ghana results in Paper 2b would invite the (correct) reviewer critique
that we cannot distinguish "framework prediction fails in Ghana" from "test
could not detect anything at this sample size."

Pipeline consequence: scripts 03, 04, 07, 10, 11, 12, and 16 have been
restricted to Florida-only execution. Scripts 13b, 14b, and 17 produce
Florida-only tables and figures. The Ghana M0 frame is still built by the
shared infrastructure (01_build_analytic_frames.R) because Paper 2a's
transport diagnostic uses it; the frame simply is not consumed by Paper 2b
output scripts.

Pre-spec consequence: the original Paper 2b pre-specification anticipated
parallel two-cohort analysis. The Ghana restriction is a transparent
deviation, disclosed in the Paper 2b Limitations and Pre-specification
Disclosure paragraph, with the operating-characteristics simulation cited
as the empirical basis for the restriction.

## Architecture

Analysis scripts 00-12 and isnv_helpers.R are **unchanged**. They produce
intermediate outputs in `data_derived/` that both papers consume.

Output generation is split into paper-specific scripts:

| Output script | Paper | Replaces |
|---|---|---|
| `R/13a_make_paper2a_tables.R` | 2a | New |
| `R/14a_make_paper2a_figures.R` | 2a | New (combines former Fig 3 + S1) |
| `R/13b_make_paper2b_tables.R` | 2b | Replaces 13_make_tables.R |
| `R/14b_make_paper2b_figures.R` | 2b | Replaces 14_make_figures.R |

Orchestrators:

| Orchestrator | Sources |
|---|---|
| `R/run_paper2a.R` | 00, 00_prep, 01, 02, 05, 06, 06b, 06c, 06d, 06e, 13a, 14a |
| `R/run_paper2b.R` | 00, 00_prep, 01, 02, 03, 04, 07, 08, 09-12, 13b, 14b |

The legacy scripts 13_make_tables.R, 14_make_figures.R, and
15_make_supplementary_figures.R can be deleted once the new pipeline is
validated. They are not invoked by either orchestrator.

## Script-to-paper allocation

| Script | Used by 2a | Used by 2b | Notes |
|---|:---:|:---:|---|
| 00_pipeline_config.R | yes | yes | Shared; extended with paper2a/2b paths |
| isnv_helpers.R | yes | yes | Shared |
| 00_prep_metadata.R | yes | yes | Shared |
| 01_build_analytic_frames.R | yes | yes | Shared |
| 02_apply_filter_ladder.R | yes | yes | Shared |
| 03_primary_pr_per_tier.R | no | yes | 2b primary |
| 04_monotonicity_test_chibar.R | no | yes | 2b primary |
| 05_transport_test_1_band.R | yes | no | 2a primary |
| 06_transport_test_2_jaccard.R | yes | no | 2a primary |
| 06b_transport_t2_hurdle.R | yes | no | 2a primary (post-hoc) |
| 06c_transport_t2_both_empty_lpm_vs_logit.R | yes | no | 2a eAppendix (scale-robustness check) |
| 06d_florida_depth_appropriateness.R | yes | no | 2a eAppendix (depth-floor derivation) |
| 06e_eappendix_tables_s2_s3.R | yes | no | 2a eAppendix (eTables S2, S3) |
| 07_test5a_depth_maf_interaction.R | no | yes | 2b primary (Edwards precondition) |
| 08_test5b_placebo_year.R | no | yes | 2b primary (placebo) |
| 09_sensitivity_year_adjusted.R | no | yes | 2b sensitivity |
| 10_sensitivity_qc_restriction.R | no | yes | 2b sensitivity |
| 11_sensitivity_outcome_ladder.R | no | yes | 2b sensitivity |
| 12_sensitivity_depth_floor_sweep.R | no | yes | 2b sensitivity |
| 13a_make_paper2a_tables.R | yes | no | 2a output |
| 13b_make_paper2b_tables.R | no | yes | 2b output |
| 14a_make_paper2a_figures.R | yes | no | 2a output |
| 14b_make_paper2b_figures.R | no | yes | 2b output |

## Output structure

```
outputs/
├── paper2a/
│   ├── tables/
│   │   ├── table1_cohort_structure.csv             (Florida cohort)
│   │   ├── table2_transport_summary.csv            (T1, T2, hurdle)
│   │   ├── eTable_florida_depth_appropriateness.csv (06d)
│   │   ├── eTable_S2_both_empty_by_depth_decile.csv (06e)
│   │   └── eTable_S3_leave_one_cluster_out.csv      (06e)
│   └── figures/
│       ├── figure1_transport.{png,pdf}             (Jaccard + both-empty)
│       ├── figure2_hurdle.{png,pdf}                (extensive + intensive)
│       ├── eFigure_S_both_empty_lpm_vs_logit.{png,pdf}     (06c)
│       └── eFigure_florida_depth_appropriateness.{png,pdf} (06d)
├── paper2b/
│   ├── tables/
│   │   ├── table1_sample_characteristics.csv   (HIV-stratified)
│   │   ├── table2_primary_pr_monotonicity.csv  (PR ladder + chi-bar test)
│   │   ├── table3_precondition_placebo.csv     (5a + 5b)
│   │   └── table4_sensitivities.csv            (year, QC, outcome, depth)
│   └── figures/
│       ├── figure1_forest.{png,pdf}            (PR per tier x cohort)
│       └── figure2_sensitivities.{png,pdf}     (4-panel sensitivities)
└── orchestrator_logs/
    └── run_paper2{a,b}_<timestamp>.log
```

## How to run

From the project root:

```r
# Paper 2a (filter transportability)
source("R/run_paper2a.R")

# Paper 2b (Edwards operationalization)
source("R/run_paper2b.R")
```

Each orchestrator independently sources the shared infrastructure
(00, 00_prep, 01, 02) before its paper-specific scripts. Running them
sequentially does no harm; the shared infrastructure is idempotent and
intermediate outputs are overwritten.

## Citation graph between papers

- Paper 1 (calibration) is cited by both 2a and 2b for the filter ladder.
- Paper 2a (transport) is cited by 2b for filter applicability evidence.
- Paper 2b does not invoke any 2a-specific computations; it inherits T2's
  passing verdict via narrative citation.

## Migration notes

1. The legacy paths `PATHS$tables`, `PATHS$figures`, and `PATHS$supplemental`
   in `00_pipeline_config.R` are preserved for backward compatibility with
   the old 13/14/15 scripts. The new pipelines do not write to these paths.
2. The legacy Figure 2 (forest, integrated paper) becomes Figure 1 in 2b.
3. The legacy Figure 3 (transport scatter, integrated paper) becomes Figure 1
   in 2a, augmented with the both-empty complementary panel.
4. The legacy Figure 4 (sensitivities, integrated paper) becomes Figure 2 in 2b.
5. The legacy Figure S1 (hurdle, supplementary) is promoted to Figure 2 in 2a.

## Expected fields in intermediate RDS objects

The Paper 2a output scripts assume the following structure in
`data_derived/06_transport_t2/transport_t2.rds`:

```
list(
  jaccard_slope        = list(estimate, ci_lower, ci_upper, p_value, n),
  both_zero_slope     = list(estimate, ci_lower, ci_upper, p_value, n),
  drop_cluster3_slope  = list(estimate, ci_lower, ci_upper, p_value),
  ghana_anchor         = list(slope_per_10x, n_pairs),
  declared             = "CRITERION MET" | "INCONCLUSIVE"
)
```

And in `data_derived/06_transport_t2/transport_t2_hurdle.rds`:

```
list(
  extensive_lpm   = list(fit, vcov, estimate, ci_lower, ci_upper, p_value, n),
  extensive_logit = list(or, or_lower, or_upper),
  intensive_ols   = list(fit, vcov, estimate, ci_lower, ci_upper, p_value, n)
)
```

If `both_zero_slope` or `drop_cluster3_slope` are not present in the
existing 06 output, the table/figure scripts will skip those rows/panels
gracefully but log a warning.
