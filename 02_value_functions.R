# ==============================================================================
# 02_value_functions.R
# Equity weighting, wEVSI, wENBS, trial costs, participant burden, and exact
# grid maximization.
# ==============================================================================

normalize_equity_weights <- function(population_share, raw_weight) {
  assert_two_group(population_share, "population_share")
  assert_two_group(raw_weight, "raw_weight")

  if (any(population_share <= 0) ||
      abs(sum(population_share) - 1) > 1e-12) {
    stop("population_share must be positive and sum to one.", call. = FALSE)
  }
  if (any(raw_weight <= 0)) {
    stop("raw_weight must be positive.", call. = FALSE)
  }

  denominator <- sum(population_share * raw_weight)
  individual_weight <- raw_weight / denominator
  population_weight <- population_share * individual_weight

  list(
    individual = individual_weight,
    population = population_weight,
    denominator = denominator
  )
}

make_batch_id <- function(nsim, n_batches) {
  nsim <- as.integer(nsim)
  n_batches <- as.integer(n_batches)

  if (nsim < 2L || n_batches < 2L || nsim %% n_batches != 0L) {
    stop("nsim must be divisible by n_batches, with both at least two.", call. = FALSE)
  }

  rep(seq_len(n_batches), each = nsim / n_batches)
}

batch_mean <- function(x, batch_id) {
  x <- as.numeric(x)
  if (length(x) != length(batch_id) || any(!is.finite(x))) {
    stop("x and batch_id must have equal finite lengths.", call. = FALSE)
  }

  sums <- rowsum(x, group = batch_id, reorder = FALSE)
  counts <- as.numeric(table(batch_id))
  as.numeric(sums[, 1] / counts)
}

current_subgroup_value <- function(theta_mean, population_weight, decision_structure) {
  decision_structure <- match.arg(decision_structure, c("common", "subgroup"))

  if (decision_structure == "common") {
    common_treat <- sum(population_weight * theta_mean) > 0
    if (common_treat) theta_mean else c(A = 0, B = 0)
  } else {
    pmax(theta_mean, 0)
  }
}

evaluate_wevsi <- function(
    evidence,
    population_share,
    raw_equity_weight,
    decision_structure,
    treatment_assignment_probability,
    n_batches) {

  decision_structure <- match.arg(decision_structure, c("common", "subgroup"))
  assert_two_group(treatment_assignment_probability, "treatment_assignment_probability")

  if (any(treatment_assignment_probability < 0) ||
      any(treatment_assignment_probability > 1)) {
    stop("Treatment-assignment probabilities must lie between zero and one.", call. = FALSE)
  }

  weights <- normalize_equity_weights(
    population_share = population_share,
    raw_weight = raw_equity_weight
  )

  alpha <- weights$population
  current_mean <- evidence$current_theta_mean
  future_mean <- evidence$future_posterior_mean
  nsim <- nrow(future_mean)
  batch_id <- make_batch_id(nsim, n_batches)

  if (decision_structure == "common") {
    current_treat <- sum(alpha * current_mean) > 0
    current_value <- max(sum(alpha * current_mean), 0)

    future_score <- as.vector(future_mean %*% alpha)
    future_treat <- future_score > 0
    future_value <- pmax(future_score, 0)

    current_contribution <- alpha * if (current_treat) {
      current_mean
    } else {
      c(A = 0, B = 0)
    }

    future_contribution <- sweep(
      future_mean,
      MARGIN = 2,
      STATS = alpha,
      FUN = "*"
    )
    future_contribution <- future_contribution * future_treat

    probability_decision_change <- mean(future_treat != current_treat)
    probability_change_A <- NA_real_
    probability_change_B <- NA_real_
  } else {
    current_treat <- current_mean > 0
    current_value <- sum(alpha * pmax(current_mean, 0))

    future_treat <- future_mean > 0
    future_contribution <- sweep(
      pmax(future_mean, 0),
      MARGIN = 2,
      STATS = alpha,
      FUN = "*"
    )
    future_value <- rowSums(future_contribution)

    current_contribution <- alpha * pmax(current_mean, 0)
    probability_change_A <- mean(future_treat[, 1] != current_treat[1])
    probability_change_B <- mean(future_treat[, 2] != current_treat[2])
    probability_decision_change <- mean(
      (future_treat[, 1] != current_treat[1]) |
        (future_treat[, 2] != current_treat[2])
    )
  }

  wEVSI_batch <- batch_mean(future_value, batch_id) - current_value
  wEVSI <- mean(wEVSI_batch)
  wEVSI_mcse <- stats::sd(wEVSI_batch) / sqrt(length(wEVSI_batch))

  subgroup_wEVSI <- colMeans(future_contribution) - current_contribution

  v_current_g <- current_subgroup_value(
    theta_mean = current_mean,
    population_weight = alpha,
    decision_structure = decision_structure
  )

  # Comparator net benefit is zero. The expected value experienced by a trial
  # participant is the treatment-assignment probability times E(theta_g).
  v_trial_g <- treatment_assignment_probability * current_mean

  delta_v_trial_participants <- sum(
    evidence$group_n *
      weights$individual *
      (v_trial_g - v_current_g)
  )

  summary_row <- data.frame(
    total_n = unname(evidence$total_n),
    p_A = unname(evidence$proportion_A_realized),
    n_A = unname(evidence$group_n["A"]),
    n_B = unname(evidence$group_n["B"]),
    decision_structure = decision_structure,
    raw_weight_A = unname(raw_equity_weight[1]),
    raw_weight_B = unname(raw_equity_weight[2]),
    normalized_weight_A = unname(weights$individual[1]),
    normalized_weight_B = unname(weights$individual[2]),
    current_theta_A = unname(current_mean[1]),
    current_theta_B = unname(current_mean[2]),
    current_value = unname(current_value),
    expected_future_value = unname(mean(future_value)),
    wEVSI = unname(wEVSI),
    wEVSI_mcse = unname(wEVSI_mcse),
    wEVSI_A = unname(subgroup_wEVSI[1]),
    wEVSI_B = unname(subgroup_wEVSI[2]),
    probability_decision_change = unname(probability_decision_change),
    probability_change_A = unname(probability_change_A),
    probability_change_B = unname(probability_change_B),
    delta_V_TP = unname(delta_v_trial_participants),
    v_current_A = unname(v_current_g[1]),
    v_current_B = unname(v_current_g[2]),
    v_trial_A = unname(v_trial_g[1]),
    v_trial_B = unname(v_trial_g[2]),
    omega_posterior_mean = unname(evidence$omega_posterior_mean),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_row,
    wEVSI_batch = wEVSI_batch
  )
}

# ------------------------------------------------------------------------------
# wENBS components
# ------------------------------------------------------------------------------

effective_population <- function(
    annual_incidence,
    horizon_cohorts,
    trial_duration_years,
    discount_rate,
    include_cohort_at_reporting_time = FALSE) {

  horizon_cohorts <- as.integer(horizon_cohorts)
  if (horizon_cohorts < 1L || annual_incidence < 0 || discount_rate < 0) {
    stop("Invalid effective-population inputs.", call. = FALSE)
  }

  times <- 0:(horizon_cohorts - 1L)
  eligible <- if (include_cohort_at_reporting_time) {
    times >= trial_duration_years
  } else {
    times > trial_duration_years
  }

  sum(annual_incidence / (1 + discount_rate)^times * eligible)
}

staircase_extra_recruitment_cost <- function(
    n_A,
    increment_usd,
    first_n_at_base_cost,
    block_size) {

  n_A <- as.integer(round(n_A))
  first_n_at_base_cost <- as.integer(first_n_at_base_cost)
  block_size <- as.integer(block_size)

  if (n_A < 0L || increment_usd < 0 || first_n_at_base_cost < 0L ||
      block_size < 1L) {
    stop("Invalid recruitment-cost inputs.", call. = FALSE)
  }

  if (n_A <= first_n_at_base_cost || increment_usd == 0) {
    return(0)
  }

  remaining <- n_A - first_n_at_base_cost
  full_blocks <- remaining %/% block_size
  partial_block <- remaining %% block_size

  full_cost <- if (full_blocks > 0L) {
    block_size * increment_usd * sum(seq_len(full_blocks))
  } else {
    0
  }

  partial_cost <- partial_block * increment_usd * (full_blocks + 1L)
  full_cost + partial_cost
}

trial_cost_usd <- function(
    n_A,
    n_B,
    fixed_cost_usd,
    base_variable_cost_usd,
    cost_increment_A_usd,
    first_A_at_base_cost,
    A_cost_block_size,
    additional_operating_cost_usd = 0) {

  fixed_cost_usd +
    base_variable_cost_usd * (n_A + n_B) +
    staircase_extra_recruitment_cost(
      n_A = n_A,
      increment_usd = cost_increment_A_usd,
      first_n_at_base_cost = first_A_at_base_cost,
      block_size = A_cost_block_size
    ) +
    additional_operating_cost_usd
}

add_wenbs_to_summary <- function(
    wevsi_summary,
    wEVSI_batch,
    trial,
    cost_increment_A_usd) {

  if (nrow(wevsi_summary) != 1L) {
    stop("wevsi_summary must contain exactly one row.", call. = FALSE)
  }

  population <- effective_population(
    annual_incidence = trial$annual_incidence,
    horizon_cohorts = trial$horizon_cohorts,
    trial_duration_years = trial$base_duration_years,
    discount_rate = trial$discount_rate,
    include_cohort_at_reporting_time = trial$include_cohort_at_reporting_time
  )

  cost_usd <- trial_cost_usd(
    n_A = wevsi_summary$n_A,
    n_B = wevsi_summary$n_B,
    fixed_cost_usd = trial$fixed_cost_usd,
    base_variable_cost_usd = trial$base_variable_cost_usd,
    cost_increment_A_usd = cost_increment_A_usd,
    first_A_at_base_cost = trial$first_A_at_base_cost,
    A_cost_block_size = trial$A_cost_block_size,
    additional_operating_cost_usd = trial$additional_operating_cost_usd
  )

  cost_qaly <- cost_usd / trial$willingness_to_pay_usd_per_qaly

  participation_burden_qaly <-
    wevsi_summary$n_A *
      wevsi_summary$normalized_weight_A *
      trial$burden_A_qaly +
    wevsi_summary$n_B *
      wevsi_summary$normalized_weight_B *
      trial$burden_B_qaly

  deterministic_component <-
    wevsi_summary$delta_V_TP -
    cost_qaly -
    participation_burden_qaly

  wENBS_batch <- population * wEVSI_batch + deterministic_component

  extra <- data.frame(
    cost_increment_A_usd = cost_increment_A_usd,
    effective_population = population,
    evidence_value_qaly = population * wevsi_summary$wEVSI,
    delta_V_TP_qaly = wevsi_summary$delta_V_TP,
    trial_cost_usd = cost_usd,
    trial_cost_qaly = cost_qaly,
    participation_burden_qaly = participation_burden_qaly,
    wENBS = mean(wENBS_batch),
    wENBS_mcse = stats::sd(wENBS_batch) / sqrt(length(wENBS_batch)),
    wEVSI_batch_fingerprint = numeric_fingerprint(wEVSI_batch),
    wENBS_batch_fingerprint = numeric_fingerprint(wENBS_batch),
    stringsAsFactors = FALSE
  )

  cbind(wevsi_summary, extra)
}

summarize_optimal_allocation <- function(data, numerical_tolerance = 1e-12) {
  required_columns <- c("p_A", "wENBS", "wENBS_mcse")
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      "Optimal-allocation data are missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(data) < 2L) {
    stop("At least two allocation designs are required.", call. = FALSE)
  }

  data <- data[order(data$p_A), , drop = FALSE]
  max_index <- which.max(data$wENBS)
  maximum <- data$wENBS[max_index]
  tied <- abs(data$wENBS - maximum) <= numerical_tolerance

  data.frame(
    trial_status = if (maximum > 0) "Trial" else "No trial",
    p_A_optimal = data$p_A[max_index],
    maximum_wENBS = maximum,
    maximum_wENBS_mcse = data$wENBS_mcse[max_index],
    n_exact_ties = sum(tied),
    stringsAsFactors = FALSE
  )
}

summarize_optima_by_group <- function(data, grouping_columns, numerical_tolerance = 1e-12) {
  missing_columns <- setdiff(c(grouping_columns, "p_A", "wENBS", "wENBS_mcse"), names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Data are missing columns needed for optimization: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  key_data <- unique(data[grouping_columns])
  output <- vector("list", nrow(key_data))

  for (i in seq_len(nrow(key_data))) {
    keep <- rep(TRUE, nrow(data))
    for (column_name in grouping_columns) {
      keep <- keep & (data[[column_name]] == key_data[[column_name]][i])
    }

    optimum <- summarize_optimal_allocation(
      data = data[keep, , drop = FALSE],
      numerical_tolerance = numerical_tolerance
    )
    output[[i]] <- cbind(key_data[i, , drop = FALSE], optimum)
  }

  do.call(rbind, output)
}
