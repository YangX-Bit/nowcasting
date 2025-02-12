library(dplyr)
library(tidyr)

simulateData <- function(
    #------------------------------
    # 1) Input control
    #------------------------------
    params = list(
      # A. Data
      data = list(
        alpha_lamb = c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6)),
        beta_lamb  = 0.5,
        N_obs       = 30,
        date_start = as.Date("2024-01-01"),
        D = 20
      ),
      # C. Model selection
      q_model = list(
        method       = "fixed_b",  # "fixed_q", "fixed_b", "linear_b", "ou_b"
        # different params for diff methods
        method_params= list(
          # fixed_q: q_D, q_lambda
          # fixed_b: b, phi
          # linear_b: beta0, beta1, phi_init, phi_sigma
          # ou_b: alpha, mu, b_init, b_sigma, phi_init, phi_sigma
          b = 0.5, phi = 0.2
        )
      )
    )
){
  #---------------------------------------------------------
  # 2) Local variebles
  #---------------------------------------------------------
  # (A) data
  alpha_lamb  <- params$data$alpha_lamb
  beta_lamb   <- params$data$beta_lamb
  N_obs       <- params$data$N_obs
  date_start  <- params$data$date_start
  D           <- params$data$D

  # (B) model for qd
  method          <- params$q_model$method
  method_params   <- params$q_model$method_params

  # other parameters
  max_delay <- 100

  #---------------------------------------------------------
  # 3) Check the input
  #---------------------------------------------------------
  if (!(length(alpha_lamb) == N_obs)) {
    stop("The length of `alpha_lamb` should be equal to the length of `N_obs`！")
  }


  #---------------------------------------------------------
  # 5) use generateQ() to generate q(d) and b_t
  #---------------------------------------------------------
    print(params$q_model)
  simsQ_out <- generateQ(
    method        = method,
    method_params = method_params,
    N_obs          = N_obs,
    D             = D,
    max_delay = max_delay
  )
    print(simsQ_out)

  #---------------------------------------------------------
  # 6) simulation
  #---------------------------------------------------------
  simulation_result <- runSimulation(
    alpha_lamb = alpha_lamb,
    beta_lamb  = beta_lamb,
    N_obs      = N_obs,
    date_start = date_start,
    simsQ_out  = simsQ_out,
    max_delay  = max_delay,
    D = D
  )

  #---------------------------------------------------------
  # 7) output
  #---------------------------------------------------------
  return(simulation_result)
}

runSimulation <- function(
    alpha_lamb,
    beta_lamb,
    N_obs,
    date_start,
    simsQ_out,
    max_delay,
    D
) {
  # Generate the date sequence
  date_seq <- seq.Date(from = date_start, by = "day", length.out = N_obs)

  # Initialize variables
  lambda_t     <- numeric(N_obs)        # Disease intensity
  case_true    <- integer(N_obs)        # True number of cases
  case_reported<- matrix(0, nrow = N_obs, ncol = max_delay + 1)  # Reported cases
  rownames(case_reported) <- as.character(date_seq)

  # Simulation process
  for (tt in seq_len(N_obs)){

    # 1) λ_t ~ Gamma
    lambda_t[tt] <- rgamma(1, shape = alpha_lamb[tt], rate = beta_lamb)

    # 2) True number of cases
    case_true[tt] <- rpois(1, lambda = lambda_t[tt])

    # 3) Get the reporting proportion for the current time point
    if (is.vector(simsQ_out$qd)){
      # If qd is a vector (constant model)
      prob_temp <- simsQ_out$qd
    } else {
      # If qd is a matrix (time-varying model, etc.)
      prob_temp <- simsQ_out$qd[tt, ]
    }

    # 4) Calculate single-day reporting proportions
    p_temp <- c(prob_temp[1], diff(prob_temp))

    # 5) Distribute the true cases to each delay day
    case_reported[tt, ] <- rmultinom(n = 1, size = case_true[tt], prob = p_temp)

  }

  qd_out <- if(is.vector(simsQ_out$qd)){
    simsQ_out$qd[1:(D+1)]
  } else {
    simsQ_out$qd[, 1:(D+1)]
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
    lambda_t   = round(lambda_t),
    # b(t) parameter
    b        = round(simsQ_out$b, 4),
    # intercept
    phi = round(simsQ_out$phi, 4),
    # q(d) reporting proportion
    qd         = round(qd_out, 4),
    # Date sequence
    date_seq   = date_seq,
    # Used D value
    D     = D
  ))
}

generateQ <- function(method, method_params, N_obs, D, max_delay = 100) {
  #----------------------------------------------------------------
  # method:         "fixed_q", "fixed_b", "rw_b", "ou_b"
  # method_params:  a list of parameters required by each method
  # N_obs:          number of observations (time points)
  # D:              maximum delay used in certain methods
  # max_delay:      size for qd columns (default 100)
  #
  # Returns a list:
  #   $qd    : either a matrix (N_obs x (max_delay+1)) or a vector
  #   $b_t   : vector (length N_obs) of b(t)
  #   $phi   : vector (length N_obs) of phi(t)
  #----------------------------------------------------------------

  # Output containers
  b_out  <- NULL
  qd_out   <- NULL
  phi_out  <- NULL

  #----------------------------------------------------------------
  # 1) fixed_q: A simple method that returns an exponential distribution q
  #----------------------------------------------------------------
  if (method == "fixed_q") {
    # Expecting something like:
    # method_params$q_D (delay) and method_params$q_lambda (rate)
    if (!all(c("q_D", "b", "phi") %in% names(method_params))) {
      stop("method=fixed_q requires 'q_D' and 'q_lambda' in method_params!")
    }
    q_D      <- method_params$q_D
    q_lambda <- method_params$q_lambda
    b <- method_params$b
    phi <- method_params$phi

    # Generate a per-day exponential distribution, summing to 1
    # qd_out <- generate_exponential_q(D = q_D, lambda = q_lambda)  # your own function
    qd_out <- 1 - phi * exp(-b * (0:q_D))
    qd_out <- qd_out / max(qd_out)

    # Return placeholders for b, phi
    b_out <- b
    phi_out <- phi

    #----------------------------------------------------------------
    # 2) fixed_b: phi and b are constants
    #----------------------------------------------------------------
  } else if (method == "fixed_b") {
    # Expecting method_params$b, method_params$phi
    if (!all(c("b", "phi") %in% names(method_params))) {
      stop("method=fixed_b requires 'b' and 'phi' in method_params!")
    }
    b   <- method_params$b
    phi <- method_params$phi

    # A single vector of size (max_delay+1): q[d] = 1 - (1-phi)*exp(-b * d)
    qd_out <- 1 - phi * exp(-b * (0:max_delay))

    # b, phi are scalars
    b_out  <- b
    phi_out  <- phi

    #----------------------------------------------------------------
    # 3) rw_b: second-order random walk for both b(t) and phi(t)
    #----------------------------------------------------------------
  } else if (method == "rw_b") {
    if (!all(c("b_intercept", "b_sigma", "phi_intercept", "phi_sigma") %in% names(method_params))) {
      stop("method=rw_b requires b_intercept, b_sigma, phi_intercept, phi_sigma!")
    }

    # Extract parameters
    b_intercept    <- method_params$b_intercept
    b_sigma   <- method_params$b_sigma
    phi_intercept  <- method_params$phi_intercept
    phi_sigma <- method_params$phi_sigma

    b <- arima.sim(model = list(order = c(0, 1, 0)), n.start = 100, n = N_obs, sd = b_sigma)
    b_out <- exp(b_intercept + b - mean(b))

    phi <- arima.sim(model = list(order = c(0, 1, 0)), n.start = 100, n = N_obs, sd = phi_sigma)
    phi_out <- plogis(phi_intercept + phi - mean(phi))

    qd_out <- matrix(NA, nrow = N_obs, ncol = max_delay + 1)
    for (i in seq_len(N_obs)) {
      b_i <- b_out[i]
      phi_i <- phi_out[i]
      qd_out[i, ] <- 1 - phi_i * exp(-b_i * (0:max_delay))
    }

  } else if (method == "ou_b") {
    if (!all(c("b_init", "b_sigma", "phi_init", "phi_sigma",
               "b_theta", "b_mu", "phi_theta", "phi_mu") %in% names(method_params))) {
      stop("method=ou_b requires b_init, b_sigma, phi_init, phi_sigma, b_theta, b_mu, phi_theta, phi_mu!")
    }

    # Extract parameters
    b_init <- method_params$b_init
    b_mu <- method_params$b_mu
    b_sigma <- method_params$b_sigma
    b_theta <- method_params$b_theta

    phi_init <- method_params$phi_init
    phi_sigma <- method_params$phi_sigma
    phi_mu <- method_params$phi_mu
    phi_theta <- method_params$phi_theta

    # # Logistic transform
    # b_star_init <- inverse_logistic_transform(b_init, 0.05, 1)
    # b_star_sigma <- b_sigma / (b_init * (1 - b_init))
    # phi_star_init <- inverse_logistic_transform(phi_init, 0, 1)
    # phi_star_sigma <- phi_sigma / (phi_init * (1 - phi_init))

    # # Transform OU targets to logistic scale
    # mu_b_star <- inverse_logistic_transform(b_mu, 0.05, 1)
    # mu_phi_star <- inverse_logistic_transform(phi_mu, 0, 1)

    # Initialize arrays
    b <- phi <- numeric(N_obs)

    # Set initial values
    b[1] <- b_init
    phi[1] <- phi_init

    # OU updates in logistic-transformed space
    for (i in 2:N_obs) {
      drift_b_star <- b_theta * (b_mu - b[i-1])
      b[i] <- b[i-1] + drift_b_star + rnorm(1, 0, b_sigma)

      drift_phi_star <- phi_theta * (phi_mu - phi[i-1])
      phi[i] <- phi[i-1] + drift_phi_star + rnorm(1, 0, phi_sigma)
    }

    # Transform back to original scale (automatically constrained)
    b_out <- exp(b)
    phi_out <- plogis(phi)

    qd_out <- matrix(NA, nrow = N_obs, ncol = max_delay + 1)
    for (i in seq_len(N_obs)) {
      b_i <- b_out[i]
      phi_i <- phi_out[i]
      qd_out[i, ] <- 1 - phi_i * exp(-b_i * (0:max_delay))
    }
  } else if (method == "sin_b") {
    #-------------------------------------------------------------
    # sin_b: Sine + noise model (final-scale baseline) for b(t) and phi(t)
    #
    # Required parameters in 'method_params':
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

    if (!all(required_params %in% names(method_params))) {
      stop(
        "method=sin_b requires the following parameters in method_params:\n",
        paste(required_params, collapse = ", ")
      )
    }

    # Extract parameters
    b_min       <- method_params$b_min
    b_max       <- method_params$b_max
    phi_min     <- method_params$phi_min
    phi_max     <- method_params$phi_max

    b_baseline  <- method_params$b_baseline
    phi_baseline<- method_params$phi_baseline

    freq        <- method_params$freq
    amp_b       <- method_params$amp_b
    amp_phi     <- method_params$amp_phi

    sigma_b     <- method_params$sigma_b
    sigma_phi   <- method_params$sigma_phi

    # Initialize output containers
    b_out   <- numeric(N_obs)
    phi_out <- numeric(N_obs)
    qd_out  <- matrix(NA_real_, nrow = N_obs, ncol = max_delay + 1)

    # Loop over time
    for (t in seq_len(N_obs)) {
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
      qd_out[t, ] <- 1 - (1 - phi_out[t]) * exp(-b_out[t] * (0:max_delay))
    }

  } else {
    stop("method must be one of: 'fixed_q', 'fixed_b', 'rw_b', 'ou_b'!")
  }

  return(list(qd = qd_out, b = b_out, phi = phi_out))
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

