# ==============================================================================
# 01_hierarchical_model_functions.R
# Posterior and posterior-predictive calculations for:
#   1) fixed omega (near-complete and near-no borrowing), and
#   2) omega ~ Half-Normal(0, scale) (moderate borrowing).
#
# Current evidence:
#   z_g | theta_g ~ Normal(theta_g, V_g0)
# Hierarchical distribution:
#   theta_g | mu, omega ~ Normal(mu, omega^2)
#   mu ~ Normal(mu0, sd_mu0^2)
# Future evidence:
#   x_g | theta_g, d ~ Normal(theta_g, V_g(d))
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_two_group <- function(x, name) {
  if (length(x) != 2L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector of length two.", call. = FALSE)
  }
  invisible(TRUE)
}

make_omega_grid <- function(
    omega_scale,
    n_grid = 301L,
    upper_multiplier = 6) {

  if (!is.finite(omega_scale) || omega_scale <= 0) {
    stop("omega_scale must be positive.", call. = FALSE)
  }

  n_grid <- as.integer(n_grid)
  if (n_grid < 51L) {
    stop("n_grid must be at least 51.", call. = FALSE)
  }

  seq(
    from = 0,
    to = upper_multiplier * omega_scale,
    length.out = n_grid
  )
}

half_normal_log_density <- function(x, scale) {
  ifelse(
    x < 0,
    -Inf,
    log(2) + stats::dnorm(x, mean = 0, sd = scale, log = TRUE)
  )
}

trapezoid_log_weights <- function(grid) {
  if (length(grid) < 2L || any(diff(grid) <= 0)) {
    stop("grid must be strictly increasing.", call. = FALSE)
  }

  increments <- diff(grid)
  if (max(abs(increments - increments[1])) > 1e-10) {
    stop("The quadrature grid must be equally spaced.", call. = FALSE)
  }

  weights <- rep(1, length(grid))
  weights[c(1, length(weights))] <- 0.5
  log(weights)
}

# ------------------------------------------------------------------------------
# Fixed-omega posterior
# ------------------------------------------------------------------------------

posterior_fixed_omega <- function(
    y,
    variance,
    omega,
    mu0 = 0,
    sd_mu0 = 1) {

  assert_two_group(y, "y")
  assert_two_group(variance, "variance")

  if (any(variance <= 0)) {
    stop("variance must be positive.", call. = FALSE)
  }
  if (!is.finite(omega) || omega <= 0) {
    stop("omega must be positive.", call. = FALSE)
  }
  if (!is.finite(sd_mu0) || sd_mu0 <= 0) {
    stop("sd_mu0 must be positive.", call. = FALSE)
  }

  omega2 <- omega^2
  prior_mu_variance <- sd_mu0^2
  marginal_variance <- variance + omega2

  mu_variance <- 1 / (
    1 / prior_mu_variance +
      sum(1 / marginal_variance)
  )

  mu_mean <- mu_variance * (
    mu0 / prior_mu_variance +
      sum(y / marginal_variance)
  )

  data_weight <- omega2 / (omega2 + variance)
  mu_weight <- variance / (omega2 + variance)

  theta_mean <- data_weight * y + mu_weight * mu_mean
  theta_conditional_variance <- variance * omega2 / (variance + omega2)

  list(
    theta_mean = theta_mean,
    mu_mean = mu_mean,
    mu_variance = mu_variance,
    theta_conditional_variance = theta_conditional_variance,
    y = y,
    variance = variance,
    omega = omega,
    mu0 = mu0,
    sd_mu0 = sd_mu0
  )
}

posterior_mean_many_fixed_omega <- function(
    y_matrix,
    variance,
    omega,
    mu0 = 0,
    sd_mu0 = 1) {

  y_matrix <- as.matrix(y_matrix)
  if (ncol(y_matrix) != 2L || any(!is.finite(y_matrix))) {
    stop("y_matrix must be a finite matrix with two columns.", call. = FALSE)
  }

  assert_two_group(variance, "variance")
  if (any(variance <= 0)) {
    stop("variance must be positive.", call. = FALSE)
  }

  omega2 <- omega^2
  prior_mu_variance <- sd_mu0^2
  marginal_variance <- variance + omega2

  mu_variance <- 1 / (
    1 / prior_mu_variance +
      sum(1 / marginal_variance)
  )

  mu_mean <- mu_variance * (
    mu0 / prior_mu_variance +
      y_matrix[, 1] / marginal_variance[1] +
      y_matrix[, 2] / marginal_variance[2]
  )

  data_weight <- omega2 / (omega2 + variance)
  mu_weight <- variance / (omega2 + variance)

  result <- cbind(
    data_weight[1] * y_matrix[, 1] + mu_weight[1] * mu_mean,
    data_weight[2] * y_matrix[, 2] + mu_weight[2] * mu_mean
  )

  colnames(result) <- c("A", "B")
  result
}

sample_theta_fixed_omega <- function(
    posterior,
    nsim,
    seed) {

  nsim <- as.integer(nsim)
  set.seed(seed)

  mu <- stats::rnorm(
    n = nsim,
    mean = posterior$mu_mean,
    sd = sqrt(posterior$mu_variance)
  )

  theta <- matrix(NA_real_, nrow = nsim, ncol = 2L)

  for (g in 1:2) {
    omega2 <- posterior$omega^2
    variance_g <- posterior$variance[g]

    data_weight <- omega2 / (omega2 + variance_g)
    mu_weight <- variance_g / (omega2 + variance_g)

    conditional_mean <-
      data_weight * posterior$y[g] +
      mu_weight * mu

    conditional_sd <- sqrt(
      posterior$theta_conditional_variance[g]
    )

    theta[, g] <- stats::rnorm(
      n = nsim,
      mean = conditional_mean,
      sd = conditional_sd
    )
  }

  colnames(theta) <- c("A", "B")
  theta
}

# ------------------------------------------------------------------------------
# Random-omega posterior: omega ~ Half-Normal(0, omega_scale)
# ------------------------------------------------------------------------------

posterior_random_omega <- function(
    y,
    variance,
    omega_scale,
    mu0 = 0,
    sd_mu0 = 1,
    omega_grid = NULL,
    n_omega = 301L,
    upper_multiplier = 6) {

  assert_two_group(y, "y")
  assert_two_group(variance, "variance")

  if (any(variance <= 0)) {
    stop("variance must be positive.", call. = FALSE)
  }

  omega_grid <- omega_grid %||% make_omega_grid(
    omega_scale = omega_scale,
    n_grid = n_omega,
    upper_multiplier = upper_multiplier
  )

  omega2 <- omega_grid^2
  prior_mu_variance <- sd_mu0^2

  s1 <- variance[1] + omega2
  s2 <- variance[2] + omega2

  a <- s1 + prior_mu_variance
  d <- s2 + prior_mu_variance
  b <- prior_mu_variance
  determinant <- a * d - b^2

  e1 <- y[1] - mu0
  e2 <- y[2] - mu0

  quadratic <- (
    d * e1^2 -
      2 * b * e1 * e2 +
      a * e2^2
  ) / determinant

  log_likelihood <- -0.5 * (
    2 * log(2 * pi) +
      log(determinant) +
      quadratic
  )

  log_weight <-
    log_likelihood +
    half_normal_log_density(omega_grid, omega_scale) +
    trapezoid_log_weights(omega_grid)

  log_weight <- log_weight - max(log_weight)
  omega_weight <- exp(log_weight)
  omega_weight <- omega_weight / sum(omega_weight)

  mu_variance <- 1 / (
    1 / prior_mu_variance +
      1 / s1 +
      1 / s2
  )

  mu_mean <- mu_variance * (
    mu0 / prior_mu_variance +
      y[1] / s1 +
      y[2] / s2
  )

  data_weight_1 <- omega2 / (omega2 + variance[1])
  data_weight_2 <- omega2 / (omega2 + variance[2])
  mu_weight_1 <- variance[1] / (omega2 + variance[1])
  mu_weight_2 <- variance[2] / (omega2 + variance[2])

  theta_mean_1 <-
    data_weight_1 * y[1] +
    mu_weight_1 * mu_mean

  theta_mean_2 <-
    data_weight_2 * y[2] +
    mu_weight_2 * mu_mean

  list(
    theta_mean = c(
      A = sum(omega_weight * theta_mean_1),
      B = sum(omega_weight * theta_mean_2)
    ),
    omega_mean = sum(omega_weight * omega_grid),
    omega_grid = omega_grid,
    omega_weight = omega_weight,
    mu_mean_by_omega = mu_mean,
    mu_variance_by_omega = mu_variance,
    y = y,
    variance = variance,
    omega_scale = omega_scale,
    mu0 = mu0,
    sd_mu0 = sd_mu0
  )
}

sample_theta_random_omega <- function(
    posterior,
    nsim,
    seed) {

  nsim <- as.integer(nsim)
  set.seed(seed)

  index <- sample.int(
    n = length(posterior$omega_grid),
    size = nsim,
    replace = TRUE,
    prob = posterior$omega_weight
  )

  omega <- posterior$omega_grid[index]
  omega2 <- omega^2

  mu <- stats::rnorm(
    n = nsim,
    mean = posterior$mu_mean_by_omega[index],
    sd = sqrt(posterior$mu_variance_by_omega[index])
  )

  theta <- matrix(NA_real_, nrow = nsim, ncol = 2L)

  for (g in 1:2) {
    variance_g <- posterior$variance[g]

    data_weight <- omega2 / (omega2 + variance_g)
    mu_weight <- variance_g / (omega2 + variance_g)

    conditional_mean <-
      data_weight * posterior$y[g] +
      mu_weight * mu

    conditional_variance <-
      variance_g * omega2 /
      (variance_g + omega2)

    theta[, g] <- stats::rnorm(
      n = nsim,
      mean = conditional_mean,
      sd = sqrt(conditional_variance)
    )
  }

  colnames(theta) <- c("A", "B")
  theta
}

posterior_mean_many_random_omega_chunk <- function(
    y_matrix,
    variance,
    omega_scale,
    omega_grid,
    mu0,
    sd_mu0) {

  y_matrix <- as.matrix(y_matrix)
  omega2 <- omega_grid^2
  prior_mu_variance <- sd_mu0^2

  s1 <- variance[1] + omega2
  s2 <- variance[2] + omega2

  a <- s1 + prior_mu_variance
  d <- s2 + prior_mu_variance
  b <- prior_mu_variance
  determinant <- a * d - b^2

  e1 <- y_matrix[, 1] - mu0
  e2 <- y_matrix[, 2] - mu0

  quadratic <-
    outer(e1^2, d / determinant) -
    2 * outer(e1 * e2, b / determinant) +
    outer(e2^2, a / determinant)

  normalizing_constant <-
    2 * log(2 * pi) +
    log(determinant)

  log_likelihood <- -0.5 * sweep(
    quadratic,
    MARGIN = 2,
    STATS = normalizing_constant,
    FUN = "+"
  )

  log_prior_quadrature <-
    half_normal_log_density(omega_grid, omega_scale) +
    trapezoid_log_weights(omega_grid)

  log_weight <- sweep(
    log_likelihood,
    MARGIN = 2,
    STATS = log_prior_quadrature,
    FUN = "+"
  )

  row_max <- apply(log_weight, 1, max)
  weight <- exp(log_weight - row_max)
  weight <- weight / rowSums(weight)

  mu_variance <- 1 / (
    1 / prior_mu_variance +
      1 / s1 +
      1 / s2
  )

  mu_numerator <-
    mu0 / prior_mu_variance +
    outer(y_matrix[, 1], 1 / s1) +
    outer(y_matrix[, 2], 1 / s2)

  mu_mean <- sweep(
    mu_numerator,
    MARGIN = 2,
    STATS = mu_variance,
    FUN = "*"
  )

  data_weight_1 <- omega2 / (omega2 + variance[1])
  data_weight_2 <- omega2 / (omega2 + variance[2])
  mu_weight_1 <- variance[1] / (omega2 + variance[1])
  mu_weight_2 <- variance[2] / (omega2 + variance[2])

  theta_mean_1 <-
    outer(y_matrix[, 1], data_weight_1) +
    sweep(mu_mean, 2, mu_weight_1, FUN = "*")

  theta_mean_2 <-
    outer(y_matrix[, 2], data_weight_2) +
    sweep(mu_mean, 2, mu_weight_2, FUN = "*")

  result <- cbind(
    rowSums(weight * theta_mean_1),
    rowSums(weight * theta_mean_2)
  )

  colnames(result) <- c("A", "B")
  result
}

posterior_mean_many_random_omega <- function(
    y_matrix,
    variance,
    omega_scale,
    omega_grid,
    mu0 = 0,
    sd_mu0 = 1,
    chunk_size = 1000L) {

  y_matrix <- as.matrix(y_matrix)

  if (ncol(y_matrix) != 2L || any(!is.finite(y_matrix))) {
    stop("y_matrix must be a finite matrix with two columns.", call. = FALSE)
  }

  n <- nrow(y_matrix)
  chunk_size <- as.integer(chunk_size)
  if (chunk_size < 1L) {
    stop("chunk_size must be at least one.", call. = FALSE)
  }
  if (length(omega_grid) < 2L || any(!is.finite(omega_grid)) ||
      any(diff(omega_grid) <= 0)) {
    stop("omega_grid must be finite and strictly increasing.", call. = FALSE)
  }
  result <- matrix(NA_real_, nrow = n, ncol = 2L)

  starts <- seq.int(1L, n, by = chunk_size)

  for (start in starts) {
    end <- min(start + chunk_size - 1L, n)
    index <- start:end

    result[index, ] <- posterior_mean_many_random_omega_chunk(
      y_matrix = y_matrix[index, , drop = FALSE],
      variance = variance,
      omega_scale = omega_scale,
      omega_grid = omega_grid,
      mu0 = mu0,
      sd_mu0 = sd_mu0
    )
  }

  colnames(result) <- c("A", "B")
  result
}


# ------------------------------------------------------------------------------
# Deterministic numeric fingerprints used to verify simulation reuse
# ------------------------------------------------------------------------------

numeric_fingerprint <- function(x) {
  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop("Fingerprint input must be nonempty and finite.", call. = FALSE)
  }

  index_weight <- (seq_along(x) %% 104729L) + 1L
  statistics <- c(
    length = length(x),
    sum = sum(x),
    sum_squares = sum(x^2),
    weighted_sum = sum(x * index_weight),
    minimum = min(x),
    maximum = max(x)
  )

  paste(
    formatC(statistics, digits = 17, format = "g"),
    collapse = "|"
  )
}

# ------------------------------------------------------------------------------
# Trial-design calculations
# ------------------------------------------------------------------------------

allocate_two_group_sample <- function(
    total_n,
    proportion_A) {

  total_n <- as.integer(round(total_n))

  if (total_n < 4L || total_n %% 2L != 0L) {
    stop("total_n must be an even integer of at least four.", call. = FALSE)
  }

  if (!is.finite(proportion_A) ||
      proportion_A <= 0 ||
      proportion_A >= 1) {
    stop("proportion_A must be strictly between zero and one.", call. = FALSE)
  }

  n_A <- 2L * as.integer(round(total_n * proportion_A / 2))
  n_A <- max(2L, min(total_n - 2L, n_A))
  n_B <- total_n - n_A

  c(A = n_A, B = n_B)
}

future_variance_two_arm <- function(
    group_n,
    outcome_sd) {

  assert_two_group(group_n, "group_n")
  assert_two_group(outcome_sd, "outcome_sd")

  if (any(group_n <= 0) || any(outcome_sd <= 0)) {
    stop("group_n and outcome_sd must be positive.", call. = FALSE)
  }

  # Equal treatment allocation within each subgroup:
  # Var(mean_treatment - mean_control) = 4 * sigma_g^2 / n_g.
  4 * outcome_sd^2 / group_n
}

combine_current_and_future <- function(
    current_mean,
    current_variance,
    future_mean_matrix,
    future_variance) {

  combined_variance <- 1 / (
    1 / current_variance +
      1 / future_variance
  )

  combined_mean <-
    sweep(
      future_mean_matrix,
      MARGIN = 2,
      STATS = 1 / future_variance,
      FUN = "*"
    )

  combined_mean <- sweep(
    combined_mean,
    MARGIN = 2,
    STATS = current_mean / current_variance,
    FUN = "+"
  )

  combined_mean <- sweep(
    combined_mean,
    MARGIN = 2,
    STATS = combined_variance,
    FUN = "*"
  )

  list(
    mean = combined_mean,
    variance = combined_variance
  )
}

prepare_current_state <- function(
    z,
    se,
    model_type,
    omega_value = NA_real_,
    omega_scale = NA_real_,
    nsim,
    seed,
    mu0 = 0,
    sd_mu0 = 1,
    n_omega = 301L,
    omega_upper_multiplier = 6,
    posterior_chunk_size = 1000L) {

  assert_two_group(z, "z")
  assert_two_group(se, "se")

  current_variance <- se^2

  if (model_type == "fixed") {
    posterior <- posterior_fixed_omega(
      y = z,
      variance = current_variance,
      omega = omega_value,
      mu0 = mu0,
      sd_mu0 = sd_mu0
    )

    theta_draws <- sample_theta_fixed_omega(
      posterior = posterior,
      nsim = nsim,
      seed = seed
    )

    omega_grid <- omega_value
    omega_weight <- 1
    omega_posterior_mean <- omega_value
  } else if (model_type == "half_normal") {
    omega_grid <- make_omega_grid(
      omega_scale = omega_scale,
      n_grid = n_omega,
      upper_multiplier = omega_upper_multiplier
    )

    posterior <- posterior_random_omega(
      y = z,
      variance = current_variance,
      omega_scale = omega_scale,
      mu0 = mu0,
      sd_mu0 = sd_mu0,
      omega_grid = omega_grid
    )

    theta_draws <- sample_theta_random_omega(
      posterior = posterior,
      nsim = nsim,
      seed = seed
    )

    omega_weight <- posterior$omega_weight
    omega_posterior_mean <- posterior$omega_mean
  } else {
    stop("model_type must be 'fixed' or 'half_normal'.", call. = FALSE)
  }

  set.seed(seed + 1L)
  standard_normal_future <- matrix(
    stats::rnorm(nsim * 2L),
    ncol = 2L
  )
  colnames(standard_normal_future) <- c("A", "B")

  list(
    model_type = model_type,
    z = z,
    se = se,
    current_variance = current_variance,
    current_theta_mean = posterior$theta_mean,
    theta_draws = theta_draws,
    standard_normal_future = standard_normal_future,
    omega_value = omega_value,
    omega_scale = omega_scale,
    omega_grid = omega_grid,
    omega_weight = omega_weight,
    omega_posterior_mean = omega_posterior_mean,
    mu0 = mu0,
    sd_mu0 = sd_mu0,
    nsim = nsim,
    seed = as.integer(seed),
    theta_draws_fingerprint = numeric_fingerprint(theta_draws),
    future_error_fingerprint = numeric_fingerprint(standard_normal_future),
    posterior_chunk_size = posterior_chunk_size
  )
}

simulate_design_evidence <- function(
    state,
    total_n,
    proportion_A,
    outcome_sd) {

  group_n <- allocate_two_group_sample(
    total_n = total_n,
    proportion_A = proportion_A
  )

  future_variance <- future_variance_two_arm(
    group_n = group_n,
    outcome_sd = outcome_sd
  )

  future_estimate <-
    state$theta_draws +
    sweep(
      state$standard_normal_future,
      MARGIN = 2,
      STATS = sqrt(future_variance),
      FUN = "*"
    )

  combined <- combine_current_and_future(
    current_mean = state$z,
    current_variance = state$current_variance,
    future_mean_matrix = future_estimate,
    future_variance = future_variance
  )

  if (state$model_type == "fixed") {
    future_posterior_mean <- posterior_mean_many_fixed_omega(
      y_matrix = combined$mean,
      variance = combined$variance,
      omega = state$omega_value,
      mu0 = state$mu0,
      sd_mu0 = state$sd_mu0
    )
  } else {
    future_posterior_mean <- posterior_mean_many_random_omega(
      y_matrix = combined$mean,
      variance = combined$variance,
      omega_scale = state$omega_scale,
      omega_grid = state$omega_grid,
      mu0 = state$mu0,
      sd_mu0 = state$sd_mu0,
      chunk_size = state$posterior_chunk_size
    )
  }

  list(
    total_n = sum(group_n),
    proportion_A_requested = proportion_A,
    proportion_A_realized = unname(group_n["A"] / sum(group_n)),
    group_n = group_n,
    future_variance = future_variance,
    future_posterior_mean = future_posterior_mean,
    current_theta_mean = state$current_theta_mean,
    omega_posterior_mean = state$omega_posterior_mean,
    simulation_seed = state$seed,
    base_draw_fingerprint = paste(
      state$theta_draws_fingerprint,
      state$future_error_fingerprint,
      sep = "::"
    ),
    design_evidence_fingerprint = numeric_fingerprint(future_posterior_mean)
  )
}
