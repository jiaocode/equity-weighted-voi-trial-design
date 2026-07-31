# ==============================================================================
# 00_analysis_parameters.R
# Prespecified inputs and reproducibility settings for the two-group
# hierarchical wEVSI / wENBS analysis.
#
# Run the workflow with:
#   Rscript --vanilla 07_run_all.R
#
# This file is the only place where the case-study inputs, simulation size,
# equity weights, recruitment-cost scenarios, and participation burdens are set.
# ==============================================================================

if (getRversion() < numeric_version("4.0.0")) {
  stop("R version 4.0.0 or newer is required.", call. = FALSE)
}

# Fix the random-number algorithms. These choices make the same seed reproduce
# the same pseudo-random sequence across supported R versions.
RNGversion("4.0.0")
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

# Use one numerical thread. The shell wrapper sets these before R starts; the
# assignments below also protect users who run an individual R script directly.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)
options(mc.cores = 1L)

ANALYSIS_SOURCE_FILES <- c(
  "00_analysis_parameters.R",
  "01_hierarchical_model_functions.R",
  "02_value_functions.R",
  "03_run_main_analysis.R",
  "04_make_main_figures.R",
  "05_run_appendix_analysis.R",
  "06_make_appendix_figures.R",
  "07_run_all.R",
  "08_run_validation_tests.R"
)

REPRODUCIBILITY <- list(
  analysis_version = "2026-07-31-v1",
  minimum_r_version = "4.0.0",
  rng_version = "4.0.0",
  rng_kind = "Mersenne-Twister",
  normal_kind = "Inversion",
  sample_kind = "Rejection",
  numeric_tolerance = 1e-10
)

compute_analysis_signature <- function(root) {
  paths <- file.path(root, ANALYSIS_SOURCE_FILES)
  missing_paths <- paths[!file.exists(paths)]

  if (length(missing_paths) > 0L) {
    stop(
      "Cannot compute the analysis signature. Missing files: ",
      paste(basename(missing_paths), collapse = ", "),
      call. = FALSE
    )
  }

  file_hashes <- unname(tools::md5sum(paths))
  payload <- c(
    paste0("analysis_version=", REPRODUCIBILITY$analysis_version),
    paste0(basename(paths), "=", file_hashes)
  )

  temporary_file <- tempfile(fileext = ".txt")
  on.exit(unlink(temporary_file), add = TRUE)
  writeLines(payload, temporary_file, useBytes = TRUE)
  unname(tools::md5sum(temporary_file))
}

runtime_information <- function() {
  rng <- RNGkind()

  list(
    r_version = R.version.string,
    platform = R.version$platform,
    rng_kind = rng[1],
    normal_kind = rng[2],
    sample_kind = rng[3],
    omp_threads = Sys.getenv("OMP_NUM_THREADS", unset = NA_character_),
    openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA_character_),
    mkl_threads = Sys.getenv("MKL_NUM_THREADS", unset = NA_character_),
    vec_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA_character_),
    blis_threads = Sys.getenv("BLIS_NUM_THREADS", unset = NA_character_)
  )
}

compute_runtime_signature <- function() {
  information <- runtime_information()
  payload <- c(
    paste0("r_version=", information$r_version),
    paste0("platform=", information$platform),
    paste0("rng_kind=", information$rng_kind),
    paste0("normal_kind=", information$normal_kind),
    paste0("sample_kind=", information$sample_kind),
    paste0("omp_threads=", information$omp_threads),
    paste0("openblas_threads=", information$openblas_threads),
    paste0("mkl_threads=", information$mkl_threads),
    paste0("vec_threads=", information$vec_threads),
    paste0("blis_threads=", information$blis_threads)
  )

  temporary_file <- tempfile(fileext = ".txt")
  on.exit(unlink(temporary_file), add = TRUE)
  writeLines(payload, temporary_file, useBytes = TRUE)
  unname(tools::md5sum(temporary_file))
}

SIMULATION <- list(
  nsim = 20000L,
  n_batches = 50L,
  allocation_step = 0.01,
  omega_grid_points = 301L,
  omega_grid_upper_multiplier = 6,
  posterior_chunk_size = 1000L,
  seed = 20260720L,
  use_common_random_numbers = TRUE
)

if (SIMULATION$nsim %% SIMULATION$n_batches != 0L) {
  stop("SIMULATION$nsim must be divisible by SIMULATION$n_batches.", call. = FALSE)
}

APPENDIX_SIMULATION <- SIMULATION
APPENDIX_SIMULATION$seed <- SIMULATION$seed + 500000L

CASE_STUDY <- list(
  group_names = c("A", "B"),
  z = c(A = 0.41, B = -0.15),
  se = c(A = 0.57, B = 0.18),
  population_share = c(A = 0.17, B = 0.83),
  total_sample_size = 600L,
  outcome_sd = c(A = 2.0, B = 2.0),
  treatment_assignment_probability = c(A = 0.5, B = 0.5),
  mu0 = 0,
  sd_mu0 = 1
)

if (abs(sum(CASE_STUDY$population_share) - 1) > 1e-12) {
  stop("Population shares must sum to one.", call. = FALSE)
}

BORROWING <- data.frame(
  borrowing_id = c("near_complete", "moderate", "near_no"),
  borrowing_axis_label = c("Near-complete", "Moderate", "Near-no"),
  borrowing_facet_label = c(
    "(a) Near-complete borrowing",
    "(b) Moderate borrowing",
    "(c) Near-no borrowing"
  ),
  model_type = c("fixed", "half_normal", "fixed"),
  omega_value = c(0.001, NA_real_, 3.0),
  omega_scale = c(NA_real_, 0.30, NA_real_),
  stringsAsFactors = FALSE
)
BORROWING$borrowing_label <- BORROWING$borrowing_facet_label

RAW_EQUITY_WEIGHT_A <- c(1.0, 1.5, 2.0)
RAW_EQUITY_WEIGHT_B <- 1.0

DECISION_STRUCTURES <- c("common", "subgroup")
DECISION_LABELS <- c(
  common = "Common population-level decision",
  subgroup = "Subgroup-specific decisions"
)

# Concise variants used as rotated row strips, where the label must fit within
# one panel height.
DECISION_ROW_LABELS <- c(
  common = "Common decision",
  subgroup = "Subgroup-specific"
)

ALLOCATION_GRID <- seq(
  from = 0.05,
  to = 0.95,
  by = SIMULATION$allocation_step
)

TRIAL <- list(
  annual_incidence = 15000,
  horizon_cohorts = 10L,
  discount_rate = 0.03,
  base_duration_years = 3,
  include_cohort_at_reporting_time = FALSE,
  fixed_cost_usd = 19000000,
  base_variable_cost_usd = 53000,
  recruitment_cost_increment_A_usd = c(0, 3000, 6000),
  first_A_at_base_cost = 50L,
  A_cost_block_size = 50L,
  additional_operating_cost_usd = 0,
  willingness_to_pay_usd_per_qaly = 150000,
  # Positive QALY-loss magnitudes. They are subtracted from wENBS.
  # The illustrative base case applies the same per-participant burden to both
  # subgroups; differential burdens are not assumed without evidence for them.
  burden_A_qaly = 0.010,
  burden_B_qaly = 0.010
)

if (any(c(TRIAL$burden_A_qaly, TRIAL$burden_B_qaly) < 0)) {
  stop(
    "Participation burdens must be stored as nonnegative QALY-loss magnitudes.",
    call. = FALSE
  )
}

# The appendix evidence-profile analyses use only the zero recruitment-cost
# increment. The recruitment-cost gradient is explored in the main analysis.
APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD <- 0

if (!all(APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD %in%
         TRIAL$recruitment_cost_increment_A_usd)) {
  stop(
    "Appendix recruitment costs must be a subset of the main cost scenarios.",
    call. = FALSE
  )
}

APPENDIX_SCENARIOS <- data.frame(
  scenario_id = c("reversed", "near_threshold", "away_threshold"),
  scenario_label = c(
    "Reversed heterogeneous evidence",
    "Same effects, uncertain, near the threshold",
    "Same effects, precise, away from the threshold"
  ),
  z_A = c(-0.15, 0.15, 0.41),
  z_B = c(0.41, 0.15, 0.41),
  se_A = c(0.18, 0.57, 0.18),
  se_B = c(0.57, 0.57, 0.18),
  stringsAsFactors = FALSE
)

OPTIMIZATION <- list(
  numerical_tolerance = 1e-12
)

PLOT <- list(
  colors = c("#0072B2", "#D55E00", "#009E73"),
  line_types = c(1, 1, 1),
  line_width = 2.4,
  maximum_point_character = 16,
  maximum_point_size = 1.25,
  population_share_line_type = 2,
  zero_line_type = 3,
  png_dpi = 300,
  # Base R shrinks all text to 0.66 of nominal size whenever a multi-panel
  # layout has three or more rows or columns. Every figure script therefore
  # resets par(cex = 1) after setting mfrow and sizes text with the values
  # below, so panel counts do not silently change the typography.
  base_pointsize = 15,
  axis_cex = 1.0,
  panel_title_cex = 1.0,
  outer_label_cex = 1.15,
  strip_label_cex = 1.1,
  legend_cex = 1.1
)

PLOT_LABELS <- list(
  equity_legend = "Equity weight for Group A",
  allocation_x = "Proportion of trial participants allocated to Group A",
  wevsi_y = "Per-person wEVSI (QALYs)",
  wenbs_y = "wENBS (QALYs)",
  recruitment_cost_row = "Recruitment-cost increment",
  decision_structure_row = "Decision structure"
)

OUTPUT <- list(
  main_results_directory = file.path("results", "main"),
  appendix_results_directory = file.path("results", "appendix"),
  main_figures_directory = file.path("figures", "main"),
  appendix_figures_directory = file.path("figures", "appendix"),
  logs_directory = file.path("results", "logs"),
  run_manifest = file.path("results", "run_manifest.csv"),
  reproducibility_check = file.path("results", "reproducibility_check.csv"),
  key_hashes = file.path("results", "key_output_hashes.csv"),
  validation_summary = file.path("results", "validation_summary.csv")
)
