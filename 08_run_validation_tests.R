# ==============================================================================
# 08_run_validation_tests.R
# Deterministic smoke tests plus validation of every main and appendix output.
# This script is run automatically after the analyses and figures.
# ==============================================================================

script_directory <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)

  if (length(file_argument) == 1L) {
    # Rscript encodes spaces in --file= as "~+~"; decode before resolving.
    dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_argument), fixed = TRUE)
    ))
  } else {
    getwd()
  }
}

ROOT <- script_directory()
source(file.path(ROOT, "00_analysis_parameters.R"), local = TRUE)
source(file.path(ROOT, "01_hierarchical_model_functions.R"), local = TRUE)
source(file.path(ROOT, "02_value_functions.R"), local = TRUE)

ANALYSIS_SIGNATURE <- compute_analysis_signature(ROOT)
TOL <- REPRODUCIBILITY$numeric_tolerance
validation_rows <- list()
validation_index <- 1L

record_test <- function(test_name, condition, detail = "") {
  condition <- isTRUE(condition)
  validation_rows[[validation_index]] <<- data.frame(
    test = test_name,
    passed = condition,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
  validation_index <<- validation_index + 1L

  if (!condition) {
    stop("Validation failed: ", test_name, ". ", detail, call. = FALSE)
  }
  invisible(TRUE)
}

check_signature <- function(data, label) {
  record_test(
    paste0(label, ": signature column"),
    "analysis_signature" %in% names(data)
  )
  signatures <- unique(as.character(data$analysis_signature))
  signatures <- signatures[!is.na(signatures) & nzchar(signatures)]
  record_test(
    paste0(label, ": current signature"),
    length(signatures) == 1L && signatures[1] == ANALYSIS_SIGNATURE,
    paste(signatures, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Parameter checks
# ------------------------------------------------------------------------------

record_test(
  "Recruitment costs are 0, 3000, 6000",
  identical(as.numeric(TRIAL$recruitment_cost_increment_A_usd), c(0, 3000, 6000))
)
record_test(
  "Group A participation burden is 0.010 QALYs",
  abs(TRIAL$burden_A_qaly - 0.010) < TOL
)
record_test(
  "Group B participation burden is 0.010 QALYs",
  abs(TRIAL$burden_B_qaly - 0.010) < TOL
)
record_test(
  "Base case applies the same participation burden to both subgroups",
  abs(TRIAL$burden_A_qaly - TRIAL$burden_B_qaly) < TOL
)
record_test(
  "Simulation uses 20000 draws",
  identical(as.integer(SIMULATION$nsim), 20000L)
)
record_test(
  "Common random numbers are required",
  isTRUE(SIMULATION$use_common_random_numbers) &&
    isTRUE(APPENDIX_SIMULATION$use_common_random_numbers)
)
record_test(
  "Allocation grid has 91 unique values",
  length(ALLOCATION_GRID) == 91L &&
    length(unique(ALLOCATION_GRID)) == 91L &&
    abs(min(ALLOCATION_GRID) - 0.05) < TOL &&
    abs(max(ALLOCATION_GRID) - 0.95) < TOL
)

# ------------------------------------------------------------------------------
# Function-level deterministic smoke tests
# ------------------------------------------------------------------------------

weights <- normalize_equity_weights(
  CASE_STUDY$population_share,
  c(A = 2, B = 1)
)
record_test(
  "Equity weights normalize to population-weighted mean one",
  abs(sum(CASE_STUDY$population_share * weights$individual) - 1) < 1e-12
)

allocation <- allocate_two_group_sample(600, 0.17)
record_test(
  "17 percent allocation gives 102 and 498 participants",
  identical(as.integer(allocation), c(102L, 498L))
)
record_test(
  "Every allocation is even and sums to 600",
  all(vapply(
    ALLOCATION_GRID,
    function(p) {
      n <- allocate_two_group_sample(600, p)
      all(n %% 2L == 0L) && sum(n) == 600L
    },
    logical(1)
  ))
)

variance_test <- future_variance_two_arm(
  group_n = c(A = 100, B = 500),
  outcome_sd = c(A = 2, B = 2)
)
record_test(
  "Two-arm future variance formula",
  max(abs(variance_test - c(A = 16 / 100, B = 16 / 500))) < 1e-14
)

record_test(
  "Staircase recruitment cost is zero for first 50 Group A participants",
  staircase_extra_recruitment_cost(50, 3000, 50, 50) == 0
)
record_test(
  "Staircase recruitment cost for 51 Group A participants",
  staircase_extra_recruitment_cost(51, 3000, 50, 50) == 3000
)
record_test(
  "Staircase recruitment cost for 150 Group A participants",
  staircase_extra_recruitment_cost(150, 3000, 50, 50) == 450000
)
record_test(
  "Staircase recruitment cost for 570 Group A participants",
  staircase_extra_recruitment_cost(570, 1, 50, 50) == 2970
)

population_manual_times <- 0:(TRIAL$horizon_cohorts - 1L)
population_manual <- sum(
  TRIAL$annual_incidence /
    (1 + TRIAL$discount_rate)^population_manual_times *
    (population_manual_times > TRIAL$base_duration_years)
)
population_function <- effective_population(
  annual_incidence = TRIAL$annual_incidence,
  horizon_cohorts = TRIAL$horizon_cohorts,
  trial_duration_years = TRIAL$base_duration_years,
  discount_rate = TRIAL$discount_rate,
  include_cohort_at_reporting_time = TRIAL$include_cohort_at_reporting_time
)
record_test(
  "Effective population formula",
  abs(population_manual - population_function) < 1e-10
)

fixed_posterior <- posterior_fixed_omega(
  y = CASE_STUDY$z,
  variance = CASE_STUDY$se^2,
  omega = 0.30,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0
)
fixed_many <- posterior_mean_many_fixed_omega(
  y_matrix = rbind(CASE_STUDY$z, CASE_STUDY$z),
  variance = CASE_STUDY$se^2,
  omega = 0.30,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0
)
record_test(
  "Vectorized fixed-omega posterior matches scalar posterior",
  max(abs(fixed_many[1, ] - fixed_posterior$theta_mean)) < 1e-12
)

random_posterior <- posterior_random_omega(
  y = CASE_STUDY$z,
  variance = CASE_STUDY$se^2,
  omega_scale = 0.30,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0,
  n_omega = 101L,
  upper_multiplier = 6
)
random_many <- posterior_mean_many_random_omega(
  y_matrix = matrix(CASE_STUDY$z, nrow = 1L),
  variance = CASE_STUDY$se^2,
  omega_scale = 0.30,
  omega_grid = random_posterior$omega_grid,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0,
  chunk_size = 1L
)
record_test(
  "Random-omega quadrature weights sum to one",
  abs(sum(random_posterior$omega_weight) - 1) < 1e-12
)
record_test(
  "Vectorized random-omega posterior matches scalar posterior",
  max(abs(random_many[1, ] - random_posterior$theta_mean)) < 1e-10
)

state_one <- prepare_current_state(
  z = CASE_STUDY$z,
  se = CASE_STUDY$se,
  model_type = "half_normal",
  omega_scale = 0.30,
  nsim = 200L,
  seed = 12345L,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0,
  n_omega = 101L,
  posterior_chunk_size = 100L
)
state_two <- prepare_current_state(
  z = CASE_STUDY$z,
  se = CASE_STUDY$se,
  model_type = "half_normal",
  omega_scale = 0.30,
  nsim = 200L,
  seed = 12345L,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0,
  n_omega = 101L,
  posterior_chunk_size = 100L
)
state_three <- prepare_current_state(
  z = CASE_STUDY$z,
  se = CASE_STUDY$se,
  model_type = "half_normal",
  omega_scale = 0.30,
  nsim = 200L,
  seed = 12346L,
  mu0 = CASE_STUDY$mu0,
  sd_mu0 = CASE_STUDY$sd_mu0,
  n_omega = 101L,
  posterior_chunk_size = 100L
)
record_test(
  "Same seed reproduces current-state draws",
  identical(state_one$theta_draws_fingerprint, state_two$theta_draws_fingerprint) &&
    identical(state_one$future_error_fingerprint, state_two$future_error_fingerprint)
)
record_test(
  "Different seed changes current-state draws",
  !identical(state_one$theta_draws_fingerprint, state_three$theta_draws_fingerprint) ||
    !identical(state_one$future_error_fingerprint, state_three$future_error_fingerprint)
)

evidence_one <- simulate_design_evidence(
  state_one,
  total_n = 600,
  proportion_A = 0.35,
  outcome_sd = CASE_STUDY$outcome_sd
)
evidence_two <- simulate_design_evidence(
  state_two,
  total_n = 600,
  proportion_A = 0.35,
  outcome_sd = CASE_STUDY$outcome_sd
)
record_test(
  "Same state and design reproduce future evidence",
  identical(
    evidence_one$design_evidence_fingerprint,
    evidence_two$design_evidence_fingerprint
  )
)

# ------------------------------------------------------------------------------
# Full-output validation
# ------------------------------------------------------------------------------

main_result_dir <- file.path(ROOT, OUTPUT$main_results_directory)
appendix_result_dir <- file.path(ROOT, OUTPUT$appendix_results_directory)

result_paths <- c(
  main_wevsi = file.path(main_result_dir, "main_wevsi_summary.csv"),
  main_wenbs = file.path(main_result_dir, "main_wenbs_summary.csv"),
  main_optimal = file.path(main_result_dir, "optimal_allocation_summary.csv"),
  main_reuse = file.path(main_result_dir, "weight_decision_draw_reuse_validation.csv"),
  main_cost_reuse = file.path(main_result_dir, "cost_draw_reuse_validation.csv"),
  appendix_wevsi = file.path(appendix_result_dir, "appendix_wevsi_summary.csv"),
  appendix_wenbs = file.path(appendix_result_dir, "appendix_wenbs_summary.csv"),
  appendix_optimal = file.path(appendix_result_dir, "appendix_optimal_allocation_summary.csv"),
  appendix_reuse = file.path(appendix_result_dir, "appendix_weight_decision_draw_reuse_validation.csv"),
  appendix_cost_reuse = file.path(appendix_result_dir, "appendix_cost_draw_reuse_validation.csv")
)
record_test("All expected result files exist", all(file.exists(result_paths)))

main_wevsi <- utils::read.csv(result_paths["main_wevsi"], stringsAsFactors = FALSE)
main_wenbs <- utils::read.csv(result_paths["main_wenbs"], stringsAsFactors = FALSE)
main_optimal <- utils::read.csv(result_paths["main_optimal"], stringsAsFactors = FALSE)
main_reuse <- utils::read.csv(result_paths["main_reuse"], stringsAsFactors = FALSE)
main_cost_reuse <- utils::read.csv(result_paths["main_cost_reuse"], stringsAsFactors = FALSE)
appendix_wevsi <- utils::read.csv(result_paths["appendix_wevsi"], stringsAsFactors = FALSE)
appendix_wenbs <- utils::read.csv(result_paths["appendix_wenbs"], stringsAsFactors = FALSE)
appendix_optimal <- utils::read.csv(result_paths["appendix_optimal"], stringsAsFactors = FALSE)
appendix_reuse <- utils::read.csv(result_paths["appendix_reuse"], stringsAsFactors = FALSE)
appendix_cost_reuse <- utils::read.csv(result_paths["appendix_cost_reuse"], stringsAsFactors = FALSE)

check_signature(main_wevsi, "main wEVSI")
check_signature(main_wenbs, "main wENBS")
check_signature(main_optimal, "main optima")
check_signature(appendix_wevsi, "appendix wEVSI")
check_signature(appendix_wenbs, "appendix wENBS")
check_signature(appendix_optimal, "appendix optima")

expected_main_wevsi_rows <-
  nrow(BORROWING) * length(ALLOCATION_GRID) *
  length(RAW_EQUITY_WEIGHT_A) * length(DECISION_STRUCTURES)
expected_main_wenbs_rows <-
  expected_main_wevsi_rows * length(TRIAL$recruitment_cost_increment_A_usd)
expected_appendix_wevsi_rows <-
  nrow(APPENDIX_SCENARIOS) * expected_main_wevsi_rows
expected_appendix_wenbs_rows <-
  expected_appendix_wevsi_rows *
  length(APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD)

record_test("Main wEVSI row count", nrow(main_wevsi) == expected_main_wevsi_rows)
record_test("Main wENBS row count", nrow(main_wenbs) == expected_main_wenbs_rows)
record_test("Appendix wEVSI row count", nrow(appendix_wevsi) == expected_appendix_wevsi_rows)
record_test("Appendix wENBS row count", nrow(appendix_wenbs) == expected_appendix_wenbs_rows)
record_test(
  "Main optimum row count",
  nrow(main_optimal) ==
    nrow(BORROWING) * length(DECISION_STRUCTURES) *
    length(RAW_EQUITY_WEIGHT_A) * length(TRIAL$recruitment_cost_increment_A_usd)
)
record_test(
  "Appendix optimum row count",
  nrow(appendix_optimal) ==
    nrow(APPENDIX_SCENARIOS) * nrow(BORROWING) *
    length(DECISION_STRUCTURES) * length(RAW_EQUITY_WEIGHT_A) *
    length(APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD)
)
record_test(
  "Appendix uses only the zero recruitment-cost increment",
  identical(as.numeric(APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD), 0) &&
    identical(sort(unique(as.numeric(appendix_wenbs$cost_increment_A_usd))), 0)
)

record_test(
  "Main reuse validation passes",
  all(main_reuse$observed_rows == main_reuse$expected_rows) &&
    all(main_reuse$draw_set_count == 1L) &&
    all(main_reuse$seed_count == 1L) &&
    all(main_reuse$base_fingerprint_count == 1L) &&
    all(main_reuse$evidence_fingerprint_count == 1L)
)
record_test(
  "Appendix reuse validation passes",
  all(appendix_reuse$observed_rows == appendix_reuse$expected_rows) &&
    all(appendix_reuse$draw_set_count == 1L) &&
    all(appendix_reuse$seed_count == 1L) &&
    all(appendix_reuse$base_fingerprint_count == 1L) &&
    all(appendix_reuse$evidence_fingerprint_count == 1L)
)
record_test(
  "Main cost reuse validation passes",
  all(main_cost_reuse$observed_costs == length(TRIAL$recruitment_cost_increment_A_usd)) &&
    max(abs(main_cost_reuse$wEVSI_range)) < 1e-12 &&
    max(abs(main_cost_reuse$wEVSI_mcse_range)) < 1e-12 &&
    all(main_cost_reuse$wEVSI_batch_fingerprint_count == 1L) &&
    max(main_cost_reuse$maximum_wENBS_increase) <= 1e-10
)
record_test(
  "Appendix cost reuse validation passes",
  all(appendix_cost_reuse$observed_costs ==
        length(APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD)) &&
    max(abs(appendix_cost_reuse$wEVSI_range)) < 1e-12 &&
    max(abs(appendix_cost_reuse$wEVSI_mcse_range)) < 1e-12 &&
    all(appendix_cost_reuse$wEVSI_batch_fingerprint_count == 1L) &&
    max(appendix_cost_reuse$maximum_wENBS_increase) <= 1e-10
)

validate_wenbs_components <- function(data, label, expected_cost_increments) {
  identity_error <- abs(
    data$wENBS -
      (
        data$effective_population * data$wEVSI +
          data$delta_V_TP -
          data$trial_cost_qaly -
          data$participation_burden_qaly
      )
  )
  expected_burden <-
    data$n_A * data$normalized_weight_A * TRIAL$burden_A_qaly +
    data$n_B * data$normalized_weight_B * TRIAL$burden_B_qaly

  record_test(
    paste0(label, ": wENBS identity"),
    max(identity_error) < 1e-8
  )
  record_test(
    paste0(label, ": participation burden"),
    max(abs(data$participation_burden_qaly - expected_burden)) < 1e-10
  )
  record_test(
    paste0(label, ": cost scenarios"),
    identical(
      sort(unique(as.numeric(data$cost_increment_A_usd))),
      sort(as.numeric(expected_cost_increments))
    )
  )
}
validate_wenbs_components(
  main_wenbs, "Main", TRIAL$recruitment_cost_increment_A_usd
)
validate_wenbs_components(
  appendix_wenbs, "Appendix", APPENDIX_RECRUITMENT_COST_INCREMENT_A_USD
)

validate_optima <- function(curves, optima, grouping_columns, label) {
  mismatch <- FALSE
  for (i in seq_len(nrow(optima))) {
    keep <- rep(TRUE, nrow(curves))
    for (column_name in grouping_columns) {
      keep <- keep & (curves[[column_name]] == optima[[column_name]][i])
    }
    subset_data <- curves[keep, , drop = FALSE]
    subset_data <- subset_data[order(subset_data$p_A), , drop = FALSE]
    index <- which.max(subset_data$wENBS)
    mismatch <- mismatch ||
      nrow(subset_data) != length(ALLOCATION_GRID) ||
      abs(subset_data$p_A[index] - optima$p_A_optimal[i]) > 1e-12 ||
      abs(subset_data$wENBS[index] - optima$maximum_wENBS[i]) > 1e-8
  }
  record_test(paste0(label, ": exact optima reproduce curves"), !mismatch)
}

validate_optima(
  main_wenbs,
  main_optimal,
  c("borrowing_id", "decision_structure", "raw_weight_A", "cost_increment_A_usd"),
  "Main"
)
validate_optima(
  appendix_wenbs,
  appendix_optimal,
  c("scenario_id", "borrowing_id", "decision_structure", "raw_weight_A", "cost_increment_A_usd"),
  "Appendix"
)

main_figure_stems <- c(
  "figure1_wevsi_common",
  "figure2_wevsi_subgroup",
  "figure3_wenbs_common_9panel",
  "figure4_wenbs_subgroup_9panel"
)
# One figure per appendix scenario: common and subgroup-specific decisions are
# shown as the two rows of a single panel grid.
appendix_figure_stems <- paste0(
  "appendix_", APPENDIX_SCENARIOS$scenario_id, "_wenbs"
)
figure_paths <- c(
  file.path(
    ROOT,
    OUTPUT$main_figures_directory,
    paste0(rep(main_figure_stems, each = 2L), c(".png", ".pdf"))
  ),
  file.path(
    ROOT,
    OUTPUT$appendix_figures_directory,
    paste0(rep(appendix_figure_stems, each = 2L), c(".png", ".pdf"))
  )
)
record_test("All expected figures exist", all(file.exists(figure_paths)))
record_test(
  "All expected figures are nonempty",
  all(file.info(figure_paths)$size > 1000)
)

key_paths <- c(
  result_paths[c(
    "main_wevsi", "main_wenbs", "main_optimal",
    "appendix_wevsi", "appendix_wenbs", "appendix_optimal"
  )],
  figure_paths
)
key_hashes <- data.frame(
  file = substring(key_paths, nchar(ROOT) + 2L),
  md5 = unname(tools::md5sum(key_paths)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  key_hashes,
  file.path(ROOT, OUTPUT$key_hashes),
  row.names = FALSE,
  na = ""
)

validation_summary <- do.call(rbind, validation_rows)
utils::write.csv(
  validation_summary,
  file.path(ROOT, OUTPUT$validation_summary),
  row.names = FALSE,
  na = ""
)

message(
  "All ", nrow(validation_summary),
  " validation tests passed."
)
