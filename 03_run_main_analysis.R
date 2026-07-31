# ==============================================================================
# 03_run_main_analysis.R
# Main-case wEVSI and wENBS analysis.
# ==============================================================================

if (!exists("ROOT", inherits = FALSE)) {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)
  ROOT <- if (length(file_argument) == 1L) {
    # Rscript encodes spaces in --file= as "~+~"; decode before resolving.
    dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_argument), fixed = TRUE)
    ))
  } else {
    getwd()
  }
}

source(file.path(ROOT, "00_analysis_parameters.R"), local = TRUE)
source(file.path(ROOT, "01_hierarchical_model_functions.R"), local = TRUE)
source(file.path(ROOT, "02_value_functions.R"), local = TRUE)

if (!isTRUE(SIMULATION$use_common_random_numbers)) {
  stop("This analysis requires common random numbers.", call. = FALSE)
}

ANALYSIS_SIGNATURE <- compute_analysis_signature(ROOT)
RESULT_DIR <- file.path(ROOT, OUTPUT$main_results_directory)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

expected_wevsi_rows <-
  nrow(BORROWING) *
  length(ALLOCATION_GRID) *
  length(RAW_EQUITY_WEIGHT_A) *
  length(DECISION_STRUCTURES)
expected_wenbs_rows <-
  expected_wevsi_rows *
  length(TRIAL$recruitment_cost_increment_A_usd)

wevsi_rows <- vector("list", expected_wevsi_rows)
wenbs_rows <- vector("list", expected_wenbs_rows)
borrowing_diagnostics <- vector("list", nrow(BORROWING))
wevsi_index <- 1L
wenbs_index <- 1L

message(
  "Running main analysis: nsim=", SIMULATION$nsim,
  ", batches=", SIMULATION$n_batches,
  ", allocations=", length(ALLOCATION_GRID),
  ", signature=", ANALYSIS_SIGNATURE
)

for (borrowing_index in seq_len(nrow(BORROWING))) {
  borrowing <- BORROWING[borrowing_index, , drop = FALSE]
  state_seed <- SIMULATION$seed + 10000L * borrowing_index

  message("  Borrowing model: ", borrowing$borrowing_axis_label)

  state <- prepare_current_state(
    z = CASE_STUDY$z,
    se = CASE_STUDY$se,
    model_type = borrowing$model_type,
    omega_value = borrowing$omega_value,
    omega_scale = borrowing$omega_scale,
    nsim = SIMULATION$nsim,
    seed = state_seed,
    mu0 = CASE_STUDY$mu0,
    sd_mu0 = CASE_STUDY$sd_mu0,
    n_omega = SIMULATION$omega_grid_points,
    omega_upper_multiplier = SIMULATION$omega_grid_upper_multiplier,
    posterior_chunk_size = SIMULATION$posterior_chunk_size
  )

  current_variance <- CASE_STUDY$se^2
  expected_shrinkage <- vapply(
    seq_along(current_variance),
    function(g) {
      sum(
        state$omega_weight *
          current_variance[g] /
          (current_variance[g] + state$omega_grid^2)
      )
    },
    numeric(1)
  )

  borrowing_diagnostics[[borrowing_index]] <- data.frame(
    borrowing_id = borrowing$borrowing_id,
    borrowing_label = borrowing$borrowing_facet_label,
    model_type = borrowing$model_type,
    omega_value = borrowing$omega_value,
    omega_scale = borrowing$omega_scale,
    omega_current_posterior_mean = state$omega_posterior_mean,
    posterior_expected_shrinkage_A = expected_shrinkage[1],
    posterior_expected_shrinkage_B = expected_shrinkage[2],
    simulation_seed = state_seed,
    base_draw_fingerprint = paste(
      state$theta_draws_fingerprint,
      state$future_error_fingerprint,
      sep = "::"
    ),
    stringsAsFactors = FALSE
  )

  simulation_draw_set_id <- paste0(
    "main_", borrowing$borrowing_id,
    "_seed_", state_seed,
    "_nsim_", SIMULATION$nsim
  )

  for (proportion_A in ALLOCATION_GRID) {
    evidence <- simulate_design_evidence(
      state = state,
      total_n = CASE_STUDY$total_sample_size,
      proportion_A = proportion_A,
      outcome_sd = CASE_STUDY$outcome_sd
    )

    for (raw_weight_A in RAW_EQUITY_WEIGHT_A) {
      raw_weight <- c(A = raw_weight_A, B = RAW_EQUITY_WEIGHT_B)

      for (decision_structure in DECISION_STRUCTURES) {
        evaluated <- evaluate_wevsi(
          evidence = evidence,
          population_share = CASE_STUDY$population_share,
          raw_equity_weight = raw_weight,
          decision_structure = decision_structure,
          treatment_assignment_probability =
            CASE_STUDY$treatment_assignment_probability,
          n_batches = SIMULATION$n_batches
        )

        summary_row <- evaluated$summary
        summary_row$borrowing_id <- borrowing$borrowing_id
        summary_row$borrowing_label <- borrowing$borrowing_facet_label
        summary_row$borrowing_axis_label <- borrowing$borrowing_axis_label
        summary_row$model_type <- borrowing$model_type
        summary_row$omega_value <- borrowing$omega_value
        summary_row$omega_scale <- borrowing$omega_scale
        summary_row$p_A_requested <- proportion_A
        summary_row$simulation_draw_set_id <- simulation_draw_set_id
        summary_row$simulation_seed <- evidence$simulation_seed
        summary_row$base_draw_fingerprint <- evidence$base_draw_fingerprint
        summary_row$design_evidence_fingerprint <-
          evidence$design_evidence_fingerprint
        summary_row$wEVSI_batch_fingerprint <-
          numeric_fingerprint(evaluated$wEVSI_batch)
        summary_row$analysis_version <- REPRODUCIBILITY$analysis_version
        summary_row$analysis_signature <- ANALYSIS_SIGNATURE
        summary_row$simulation_nsim <- SIMULATION$nsim
        summary_row$simulation_n_batches <- SIMULATION$n_batches

        wevsi_rows[[wevsi_index]] <- summary_row
        wevsi_index <- wevsi_index + 1L

        for (cost_increment in TRIAL$recruitment_cost_increment_A_usd) {
          wenbs_row <- add_wenbs_to_summary(
            wevsi_summary = summary_row,
            wEVSI_batch = evaluated$wEVSI_batch,
            trial = TRIAL,
            cost_increment_A_usd = cost_increment
          )
          wenbs_rows[[wenbs_index]] <- wenbs_row
          wenbs_index <- wenbs_index + 1L
        }
      }
    }
  }
}

wevsi_summary <- do.call(rbind, wevsi_rows)
wenbs_summary <- do.call(rbind, wenbs_rows)
row.names(wevsi_summary) <- NULL
row.names(wenbs_summary) <- NULL

if (nrow(wevsi_summary) != expected_wevsi_rows ||
    wevsi_index != expected_wevsi_rows + 1L) {
  stop("Unexpected number of main wEVSI rows.", call. = FALSE)
}
if (nrow(wenbs_summary) != expected_wenbs_rows ||
    wenbs_index != expected_wenbs_rows + 1L) {
  stop("Unexpected number of main wENBS rows.", call. = FALSE)
}

optimization_columns <- c(
  "borrowing_id",
  "borrowing_label",
  "borrowing_axis_label",
  "decision_structure",
  "raw_weight_A",
  "cost_increment_A_usd"
)

optimal_allocation <- summarize_optima_by_group(
  data = wenbs_summary,
  grouping_columns = optimization_columns,
  numerical_tolerance = OPTIMIZATION$numerical_tolerance
)
optimal_allocation$analysis_version <- REPRODUCIBILITY$analysis_version
optimal_allocation$analysis_signature <- ANALYSIS_SIGNATURE
optimal_allocation$simulation_nsim <- SIMULATION$nsim
optimal_allocation$simulation_n_batches <- SIMULATION$n_batches

# ------------------------------------------------------------------------------
# Internal validation before writing outputs
# ------------------------------------------------------------------------------

normalized_weight_mean <-
  CASE_STUDY$population_share["A"] * wevsi_summary$normalized_weight_A +
  CASE_STUDY$population_share["B"] * wevsi_summary$normalized_weight_B
if (max(abs(normalized_weight_mean - 1)) > 1e-10) {
  stop("Equity-weight normalization failed.", call. = FALSE)
}

identity_error <- abs(
  wenbs_summary$wENBS -
    (
      wenbs_summary$effective_population * wenbs_summary$wEVSI +
        wenbs_summary$delta_V_TP -
        wenbs_summary$trial_cost_qaly -
        wenbs_summary$participation_burden_qaly
    )
)
if (max(identity_error) > 1e-8) {
  stop("The wENBS component identity failed.", call. = FALSE)
}

expected_burden <-
  wenbs_summary$n_A *
    wenbs_summary$normalized_weight_A *
    TRIAL$burden_A_qaly +
  wenbs_summary$n_B *
    wenbs_summary$normalized_weight_B *
    TRIAL$burden_B_qaly
if (max(abs(wenbs_summary$participation_burden_qaly - expected_burden)) > 1e-10) {
  stop("The participant-burden calculation failed.", call. = FALSE)
}

observed_costs <- sort(unique(wenbs_summary$cost_increment_A_usd))
expected_costs <- sort(as.numeric(TRIAL$recruitment_cost_increment_A_usd))
if (length(observed_costs) != length(expected_costs) ||
    any(abs(observed_costs - expected_costs) > 1e-12)) {
  stop("The main recruitment-cost scenarios are inconsistent.", call. = FALSE)
}

# Same evidence must be reused across weights and decision structures.
reuse_key <- interaction(
  wevsi_summary$borrowing_id,
  sprintf("%.2f", wevsi_summary$p_A),
  drop = TRUE
)
reuse_groups <- split(wevsi_summary, reuse_key)
reuse_validation <- vector("list", length(reuse_groups))
reuse_names <- names(reuse_groups)

for (i in seq_along(reuse_groups)) {
  group_data <- reuse_groups[[i]]
  reuse_validation[[i]] <- data.frame(
    group_key = reuse_names[i],
    observed_rows = nrow(group_data),
    expected_rows =
      length(RAW_EQUITY_WEIGHT_A) * length(DECISION_STRUCTURES),
    draw_set_count = length(unique(group_data$simulation_draw_set_id)),
    seed_count = length(unique(group_data$simulation_seed)),
    base_fingerprint_count =
      length(unique(group_data$base_draw_fingerprint)),
    evidence_fingerprint_count =
      length(unique(group_data$design_evidence_fingerprint)),
    stringsAsFactors = FALSE
  )
}
reuse_validation <- do.call(rbind, reuse_validation)

if (any(reuse_validation$observed_rows != reuse_validation$expected_rows) ||
    any(reuse_validation$draw_set_count != 1L) ||
    any(reuse_validation$seed_count != 1L) ||
    any(reuse_validation$base_fingerprint_count != 1L) ||
    any(reuse_validation$evidence_fingerprint_count != 1L)) {
  stop("Simulation reuse failed across equity weights or decisions.", call. = FALSE)
}

# Cost scenarios must reuse the same wEVSI scalar and batch vector.
cost_key <- interaction(
  wenbs_summary$borrowing_id,
  wenbs_summary$decision_structure,
  wenbs_summary$raw_weight_A,
  sprintf("%.2f", wenbs_summary$p_A),
  drop = TRUE
)
cost_groups <- split(wenbs_summary, cost_key)
cost_validation <- vector("list", length(cost_groups))
cost_names <- names(cost_groups)

for (i in seq_along(cost_groups)) {
  group_data <- cost_groups[[i]]
  group_data <- group_data[order(group_data$cost_increment_A_usd), , drop = FALSE]
  cost_validation[[i]] <- data.frame(
    group_key = cost_names[i],
    observed_costs = length(unique(group_data$cost_increment_A_usd)),
    wEVSI_range = max(group_data$wEVSI) - min(group_data$wEVSI),
    wEVSI_mcse_range = max(group_data$wEVSI_mcse) - min(group_data$wEVSI_mcse),
    wEVSI_batch_fingerprint_count =
      length(unique(group_data$wEVSI_batch_fingerprint)),
    maximum_wENBS_increase = max(diff(group_data$wENBS)),
    stringsAsFactors = FALSE
  )
}
cost_validation <- do.call(rbind, cost_validation)

if (any(cost_validation$observed_costs != length(expected_costs)) ||
    any(abs(cost_validation$wEVSI_range) > 1e-12) ||
    any(abs(cost_validation$wEVSI_mcse_range) > 1e-12) ||
    any(cost_validation$wEVSI_batch_fingerprint_count != 1L) ||
    any(cost_validation$maximum_wENBS_increase > 1e-10)) {
  stop("Simulation reuse or monotonicity failed across cost scenarios.", call. = FALSE)
}

# With a higher allocation-dependent cost, the lowest exact maximizing p_A
# cannot move upward.
optimal_key <- interaction(
  optimal_allocation$borrowing_id,
  optimal_allocation$decision_structure,
  optimal_allocation$raw_weight_A,
  drop = TRUE
)
optimal_groups <- split(optimal_allocation, optimal_key)
for (group_data in optimal_groups) {
  group_data <- group_data[order(group_data$cost_increment_A_usd), , drop = FALSE]
  if (any(diff(group_data$p_A_optimal) > 1e-12)) {
    stop("An exact optimal allocation increased with recruitment cost.", call. = FALSE)
  }
}

utils::write.csv(
  wevsi_summary,
  file.path(RESULT_DIR, "main_wevsi_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  wenbs_summary,
  file.path(RESULT_DIR, "main_wenbs_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  optimal_allocation,
  file.path(RESULT_DIR, "optimal_allocation_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  do.call(rbind, borrowing_diagnostics),
  file.path(RESULT_DIR, "borrowing_diagnostics.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  reuse_validation,
  file.path(RESULT_DIR, "weight_decision_draw_reuse_validation.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  cost_validation,
  file.path(RESULT_DIR, "cost_draw_reuse_validation.csv"),
  row.names = FALSE,
  na = ""
)

saveRDS(wevsi_summary, file.path(RESULT_DIR, "main_wevsi_summary.rds"))
saveRDS(wenbs_summary, file.path(RESULT_DIR, "main_wenbs_summary.rds"))

message("Main analysis completed successfully: ", RESULT_DIR)
