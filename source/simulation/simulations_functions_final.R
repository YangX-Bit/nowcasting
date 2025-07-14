library(dplyr)
library(tidyr)

simulateData <- function(
    params = list(
      data = list(
        alpha_lamb = c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6)),
        beta_lamb  = 0.5,
        T       = 30,
        date_start = as.Date("2024-01-01"),
        D = 20
      ),
      q_model = list(
        method       = "b_constant",  # "q_constant", "b_constant", "b_rw", "b_ou"
        method_params= list(
          # q_constant: b, phi
          # b_constant: b, phi
          # b_rw: mu_logb, sigma_logb, mu_logitphi, sigma_logitphi
          # b_ou: init_logb, mu_logb, sigma_logb, theta_logb,
          #     theta_logitphi, init_logitphi, mu_logitphi, sigma_logitphi
          b = 0.5, phi = 0.2
        )
      )
    )
){
  # (A) data
  alpha_lamb  <- params$data$alpha_lamb
  beta_lamb   <- params$data$beta_lamb
  T       <- params$data$T
  date_start  <- params$data$date_start
  D           <- params$data$D

  # (B) model for qd
  method          <- params$q_model$method
  method_params   <- params$q_model$method_params

  if (!(length(alpha_lamb) == T)) {
    stop("The length of `alpha_lamb` should be equal to the length of `T`！")
  }

  simsQ_out <- generateQ(method = method, params = method_params, T = T, D = D)

  simulation_result <- runSimulation(
    alpha_lamb = alpha_lamb, beta_lamb  = beta_lamb, T = T, date_start = date_start,
    simsQ_out  = simsQ_out, D = D)

  return(simulation_result)
}

runSimulation <- function(
    alpha_lamb,
    beta_lamb,
    T,
    date_start,
    simsQ_out,
    D
) {
  # Generate the date sequence
  date_seq <- seq.Date(from = date_start, by = "day", length.out = T)

  # Initialize variables
  lambda     <- numeric(T)        # Disease intensity
  case_true    <- integer(T)        # True number of cases
  case_reported<- matrix(0, nrow = T, ncol = D + 1)  # Reported cases
  rownames(case_reported) <- as.character(date_seq)

  # Simulation process
  for (tt in seq_len(T)){

    # 1) λ_t ~ Gamma
    lambda[tt] <- rgamma(1, shape = alpha_lamb[tt], rate = beta_lamb)

    # 2) True number of cases
    case_true[tt] <- rpois(1, lambda = lambda[tt])

    # 3) Get the reporting proportion for the current time point
    if (is.vector(simsQ_out$q)){
      # If q is a vector (constant model)
      prob_temp <- simsQ_out$q
    } else {
      # If q is a matrix (time-varying model, etc.)
      prob_temp <- simsQ_out$q[tt, ]
    }

    # 4) Calculate single-day reporting proportions
    p_temp <- c(prob_temp[1], diff(prob_temp), 1 - prob_temp[D+1])

    # 5) Distribute the true cases to each delay day
    case_reported[tt, ] <- rmultinom(n = 1, size = case_true[tt], prob = p_temp)[1:(D+1)]

  }

  qd_out <- if(is.vector(simsQ_out$q)){
    simsQ_out$q[1:(D+1)]
  } else {
    simsQ_out$q[, 1:(D+1)]
  }

  # Convert true cases to matrix and set row names
  case_true = as.matrix(case_true)
  rownames(case_true) = as.character(date_seq)

  # cumulated reported cases
  case_reported_cumulated <- t(apply(case_reported, 1, cumsum))

  # Return the final result list
  return(list(
    # Reported cases
    case_reported = case_reported[, c(1:(D+1))],
    case_reported_cumulated = case_reported_cumulated[, c(1:(D+1))],
    # True cases
    case_true  = case_true,
    # Disease intensity
    lambda   = round(lambda, 4),
    # b(t) parameter
    b        = round(simsQ_out$b, 4),
    # intercept
    phi = round(simsQ_out$phi, 4),
    # q(d) reporting proportion
    q         = round(qd_out, 4),
    # Date sequence
    date_seq   = date_seq,
    # Used D value
    D     = D
  ))
}

generateQ <- function(method, params, T, D) {
  #----------------------------------------------------------------
  # method:         "q_constant", "b_constant", "b_rw", "b_ou"
  # params:  a list of parameters required by each method
  # T:          number of observations (time points)
  # D:              maximum delay used in certain methods
  #
  # Returns a list:
  #   $q    : either a matrix (T x (max_delay+1)) or a vector
  #   $b_t   : vector (length T) of b(t)
  #   $phi   : vector (length T) of phi(t)
  #----------------------------------------------------------------

  if (method %in% c("q_constant", "b_constant")) {

    if (!all(c("b", "phi") %in% names(params))) {
      stop("method=q_constant,b_constant requires 'b' and 'phi' in params!")
    }
    b <- params$b
    phi <- params$phi
    q <- 1 - (1 - phi) * exp(-b * (0:D))

  } else if (method == "b_rw") {

    if (!all(c("mu_logb", "sigma_logb", "mu_logitphi", "sigma_logitphi") %in% names(params))) {
      stop("method=b_rw requires mu_logb, sigma_logb, mu_logitphi, sigma_logitphi!")
    }

    logb <- arima.sim(model = list(order = c(0, 1, 0)), sd = params$sigma_logb,
            n = T, n.start = 100) |> as.numeric()
    logitphi <- arima.sim(model = list(order = c(0, 1, 0)), sd = params$sigma_logitphi,
            n = T, n.start = 100) |> as.numeric()

    b <- exp(params$mu_logb + logb - mean(logb))
    phi <- plogis(params$mu_logitphi + logitphi - mean(logitphi))

    q <- matrix(NA, nrow = T, ncol = D + 1)
    for (i in seq_len(T)) {
      q[i, ] <- 1 - (1 - phi[i]) * exp(-b[i] * (0:D))
    }

  } else if (method == "b_ou") {

    if (!all(c("init_logb", "mu_logb", "sigma_logb", "theta_logitphi",
            "init_logitphi", "mu_logitphi", "sigma_logitphi", "theta_logb") %in% names(params))) {
      stop("method=b_ou requires init_logb, mu_logb, sigma_logb, theta_logb, init_logitphi, mu_logitphi, sigma_logitphi, theta_logitphi!")
    }

    logb <- logitphi <- numeric(T)
    logb[1] <- params$init_logb
    logitphi[1] <- params$init_logitphi

    for (i in 2:T) {
      drift_b_star <- params$theta_logb * (params$mu_logb - logb[i-1])
      logb[i] <- logb[i-1] + drift_b_star + rnorm(1, 0, params$sigma_logb)

      drift_phi_star <- params$theta_logitphi * (params$mu_logitphi - logitphi[i-1])
      logitphi[i] <- logitphi[i-1] + drift_phi_star + rnorm(1, 0, params$sigma_logitphi)
    }

    b <- exp(logb)
    phi <- plogis(logitphi)

    q <- matrix(NA, nrow = T, ncol = D + 1)
    for (i in seq_len(T)) {
      q[i, ] <- 1 - (1 - phi[i]) * exp(-b[i] * (0:D))
    }

  } else if (method == "sin_b") {
    #-------------------------------------------------------------
    # sin_b: Sine + noise model (final-scale baseline) for b(t) and phi(t)
    #
    # Required parameters in 'params':
    #   1) b_min, b_max        : Lower and upper bounds for b(t)
    #   2) phi_min, phi_max    : Lower and upper bounds for phi(t)
    #
    #   3) b_baseline, phi_baseline : Baseline (offset) in the final scale
    #
    #   4) freq                : Sine wave frequency (shared by both b and phi)
    #   5) amp_b, amp_phi      : Amplitudes for b(t) and phi(t) in final scale
    #   6) sigma_b, sigma_phi  : Std dev of Gaussian noise for b(t) and phi(t)
    #
    # The final b(t) and phi(t) are clamped to [b_min,b_max], [phi_min,phi_max].
    # This approach avoids the [-1,1] -> [x_min,x_max] mapping,
    # so your baseline is directly in the final scale.
    #-------------------------------------------------------------

    required_params <- c(
      "b_min", "b_max",
      "phi_min", "phi_max",
      "b_baseline", "phi_baseline",
      "freq",
      "amp_b", "amp_phi",
      "sigma_b", "sigma_phi"
    )

    if (!all(required_params %in% names(params))) {
      stop(
        "method=sin_b requires the following parameters in params:\n",
        paste(required_params, collapse = ", ")
      )
    }

    # Extract parameters
    b_min       <- params$b_min
    b_max       <- params$b_max
    phi_min     <- params$phi_min
    phi_max     <- params$phi_max

    b_baseline  <- params$b_baseline
    phi_baseline<- params$phi_baseline

    freq        <- params$freq
    amp_b       <- params$amp_b
    amp_phi     <- params$amp_phi

    sigma_b     <- params$sigma_b
    sigma_phi   <- params$sigma_phi

    # Initialize output containers
    b_out   <- numeric(T)
    phi_out <- numeric(T)
    qd_out  <- matrix(NA_real_, nrow = T, ncol = D + 1)

    # Loop over time
    for (t in seq_len(T)) {
      # 1) Compute the raw b(t) and phi(t) in final scale
      #    baseline + sine + noise
      b_raw   <- b_baseline   + amp_b   * sin(2 * pi * freq * t) +
        rnorm(1, mean = 0, sd = sigma_b)
      phi_raw <- phi_baseline + amp_phi * sin(2 * pi * freq * t) +
        rnorm(1, mean = 0, sd = sigma_phi)

      # 2) Clamp to [b_min, b_max], [phi_min, phi_max] to avoid out-of-bound values
      b_out[t]   <- max(b_min, min(b_max, b_raw))
      phi_out[t] <- max(phi_min, min(phi_max, phi_raw))

      # 3) Build q(d): q[d] = 1 - (1 - phi(t)) * exp(-b(t)*d)
      qd_out[t, ] <- 1 - (1 - phi_out[t]) * exp(-b_out[t] * (0:D))
    }

  } else {
    stop("method must be one of: 'q_constant', 'b_constant', 'b_rw', 'b_ou'!")
  }

  return(list(q = q, b = b, phi = phi))
}


## generate exponential decay q. A quick easy way
# - why lambda? use another name
# - you should add an intercept here
generate_exponential_q <- function(D, lambda = 0.3) {
  q <- exp(-lambda * (1:(D+1)))
  q <- q / sum(q)  # Normalize to sum to 1
  q_out <- cumsum(q)
  return(q_out)
}

# Helper functions for transformations
logistic_transform <- function(x, lower, upper) {
  # Transform from real line to (lower, upper) interval
  lower + (upper - lower) * exp(x) / (1 + exp(x))
}

inverse_logistic_transform <- function(y, lower, upper) {
  # Transform from (lower, upper) to real line
  log((y - lower) / (upper - y))
}




### functions to transfer the simulation data to the form of data in the paper
# Parameters:
#  data - Matrix of cases in each day with delays
###
triangleToList <- function(data,
                          now = as.Date("2011-07-04")){

  # sequence of the start date
  admission_dates <- sort(now - 0:(nrow(data) - 1), descending = T)


  df <- as.data.frame(data)
  D <- ncol(data) - 1
  colnames(df) <- paste0("delay", 0:D)

  # long data
  long_df <- df %>%
    mutate(admission_date = admission_dates) %>%
    pivot_longer(cols = starts_with("delay"),
                 names_to = "delay",
                 names_prefix = "delay",
                 values_to = "reported_cases") %>%
    filter(reported_cases > 0) %>%
    mutate(delay = as.numeric(delay),
           report_date = admission_date + delay) %>%
    uncount(reported_cases)

  data_out <- long_df %>% mutate( dHosp= admission_date,
                                     dReport= report_date) %>%
    select(dHosp, dReport) %>%
    as.data.frame()

  return(data_out)
}

### functions to transfer the list data from the paper to triangular data
# Parameters:
#  data - List of data from the form of the paper.
###
listToTriangle <- function(data, now = NULL, D = NULL) {
  # Ensure the data has the required columns
  if (!all(c("dHosp", "dReport") %in% colnames(data))) {
    stop("Input data must have columns 'dHosp' and 'dReport'.")
  }

  # Calculate delay
  data <- data %>%
    mutate(delay = as.numeric(dReport - dHosp))  # Compute delay in days

  # Find the maximum delay (D) and range of admission dates
  if(is.null(D)){ D <- max(data$delay) }  # Maximum delay
  if(is.null(now)){now <- max(data$dHosp)}

  admission_dates <- seq(min(data$dHosp), now, by = "1 day")  # Generate dates from 'now' backward

  # Initialize a matrix for triangular data
  triangular_matrix <- matrix(0, nrow = length(admission_dates), ncol = D + 1)
  rownames(triangular_matrix) <- as.character(admission_dates)
  colnames(triangular_matrix) <- paste0("delay", 0:D)

  # Populate the matrix directly
  for (i in 1:nrow(data)) {
    hosp_idx <- which(admission_dates == data$dHosp[i])  # Row index for admission date
    delay_idx <- data$delay[i] + 1  # Column index for delay (convert 0-based to 1-based indexing)
    triangular_matrix[hosp_idx, delay_idx] <- triangular_matrix[hosp_idx, delay_idx] + 1
  }

  return(triangular_matrix)
}

########## output the cumulatived matrix by col
cumulative_matrix <- function(mat) {
  # Check if input is a matrix
  if (!is.matrix(mat)) {
    stop("Input must be a matrix.")
  }

  # Get the number of rows and columns
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)

  # Loop through each column starting from the second
  for (j in 2:n_cols) {
    # Add the previous column to the current column
    mat[, j] <- mat[, j] + mat[, j - 1]
  }

  # Return the updated matrix
  return(mat)
}


### functions to transfer full data to truncated triangular data
create_triangular_data <- function(Y_full, if_zero = FALSE) {
  N <- nrow(Y_full)     # number of days
  D <- ncol(Y_full) - 1 # Max Delay

  # # check if N <= D
  # if (N <= (D + 1)) {
  #   stop("The number of rows (N) cannot be smaller than D + 1.")
  # }

  # out matrix
  Y_triangular <- matrix(NA, nrow = N, ncol = D + 1)

  if(N > D){
    for (i in 1:(N-D)) {
      # keeps the full data
      Y_triangular[i, ] <- Y_full[i, ]
      # Y_triangular[N-i+1, 1:i] <- Y_full[N-i+1, 1:i]
    }

    for (j in 1:D) {
      # keeps
      Y_triangular[N-j+1, 1:j] <- Y_full[N-j+1, 1:j]
    }
  }else{
    for (i in seq_len(N)) {
      Y_triangular[N - i + 1, 1:i] <- Y_full[N - i + 1, 1:i]
    }
  }

  if(if_zero){
    Y_triangular[is.na(Y_triangular)] <- 0
  }

  return(Y_triangular)
}

extract_last_valid <- function(mat, D = ncol(mat)) {
  # Number of rows in the matrix
  N <- nrow(mat)

  # Initialize a vector to store the result
  last_valid <- numeric(N)

  # If sample size is less than D, search in all columns
  if (N < D) {
    for (i in 1:N) {
      # Traverse from right to left to find the first valid value
      for (j in ncol(mat):1) {
        if (!is.na(mat[i, j])) {
          last_valid[i] <- mat[i, j]
          break
        }
      }
    }
    return(last_valid)
  }

  # Handle rows greater than or equal to D
  # Process the first N-D complete rows
  for (i in 1:(N - D + 1)) {
    # Traverse from right to left to find the last valid value
    for (j in ncol(mat):1) {
      if (!is.na(mat[i, j])) {
        last_valid[i] <- mat[i, j]
        break
      }
    }
  }

  # Handle the last D rows with incomplete data
  for (i in (N - D + 1):N) {
    # Calculate offset
    offset <- N - i

    # Traverse from right to left across all columns
    for (j in ncol(mat):1) {
      if (!is.na(mat[i, j])) {
        last_valid[i] <- mat[i, j]
        break
      }
    }
  }

  return(last_valid)
}

