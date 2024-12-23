library(dplyr)
library(tidyr)

generate_alpha <- function(n, range = c(0, 1), dist = "uniform", seed = NULL) {
  # check the range 
  if (length(range) != 2 || range[1] >= range[2]) {
    stop("range must be a vector of length 2 with range[1] < range[2].")
  }
  
  # set random seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # produce random alpha value
  if (dist == "uniform") {
    # uniform
    alpha <- runif(n, min = range[1], max = range[2])
  } else if (dist == "normal") {
    # normal
    mean_val <- mean(range)  # mean is the middle value of the range
    sd_val <- (range[2] - range[1]) / 4  # SD be the 1/4 of the range
    alpha <- rnorm(n, mean = mean_val, sd = sd_val)
    
    # 
    alpha <- pmax(pmin(alpha, range[2]), range[1])
  } else {
    stop("Unsupported distribution. Choose 'uniform' or 'normal'.")
  }
  
  return(alpha)
}


simsQ <- function(
    method = c("constant", "time_varying", "random_walk"),
    D = NULL,
    t = NULL, 
    # Parameters for different methods:
    # For "constant": provide 'b'
    # For "time_varying": provide 'beta0', 'beta1'
    # For "random_walk": provide 'b' (initial), 'sigma_rw'
    b = NULL, 
    beta0 = NULL, beta1 = NULL,
    sigma_rw = NULL,
    ensure_qd_one = TRUE
) {
  # 'method' defines how Q is simulated:
  #   "constant": Q(d) = 1 - exp(-b * d)
  #   "time_varying": b_t = exp(beta0 + beta1 * t), Q_t(d) = 1 - exp(-b_t[i] * d)
  #   "random_walk": b_t evolves as a random walk in log-space, Q_t(d) = 1 - exp(-b_t[i] * d)
  #
  # 'ensure_qd_one': if TRUE, after computing qd, 
  # we normalize so that qd(..., D+1) = 1 exactly.
  #
  # Returns a list with:
  #   qd: 
  #     - If method="constant", a vector of length D+1
  #     - Otherwise, a matrix (t x (D+1))
  #   b_t:
  #     - The b_t values used. For "constant", b_t is just a single value.
  #     - For "time_varying" and "random_walk", b_t is a vector of length t.
  
  method <- match.arg(method)
  
  # A helper function to normalize qd so that qd at D+1 is exactly 1.
  normalize_qd <- function(qd) {
    if (is.vector(qd)) {
      # For a single vector
      final_val <- qd[D+1]
      if (final_val > 0) qd <- qd / final_val
    } else {
      # For a matrix, normalize each row
      final_vals <- qd[, D+1]
      # Avoid division by zero
      valid_rows <- final_vals > 0
      qd[valid_rows, ] <- qd[valid_rows, ] / final_vals[valid_rows]
    }
    return(qd)
  }
  
  if (method == "constant") {
    if (is.null(b)) stop("For 'constant' method, 'b' must be provided.")
    # No need for t
    qd <- 1 - exp(-b * (1:(D+1)))
    b_t <- b
    
  } else if (method == "time_varying") {
    if (is.null(t)) stop("For 'time_varying' method, 't' must be provided.")
    if (is.null(beta0) || is.null(beta1)) {
      stop("For 'time_varying' method, 'beta0' and 'beta1' must be provided.")
    }
    qd <- matrix(NA, nrow = t, ncol = D + 1)
    b_t <- exp(beta0 + beta1 * (1:t))
    for (i in 1:t) {
      qd[i, ] <- 1 - exp(-b_t[i] * (1:(D+1)))
    }
    
  } else { # method == "random_walk"
    if (is.null(t)) stop("For 'random_walk' method, 't' must be provided.")
    if (is.null(b) || is.null(sigma_rw)) {
      stop("For 'random_walk' method, 'b' (initial) and 'sigma_rw' must be provided.")
    }
    qd <- matrix(NA, nrow = t, ncol = D + 1)
    b_t <- numeric(t)
    b_t[1] <- b
    # Evolve b_t as a random walk in log-space
    for (i in 2:t) {
      b_t[i] <- exp(log(b_t[i - 1]) + rnorm(1, mean = 0, sd = sigma_rw))
    }
    # Compute qd
    for (i in 1:t) {
      qd[i, ] <- 1 - exp(-b_t[i] * (1:(D+1)))
    }
  }
  
  # Normalize qd if required
  if (ensure_qd_one) {
    qd <- normalize_qd(qd)
  }
  
  return(list(qd = qd, b_t = b_t))
}
# set.seed(1)
# simsQ(method = "random_walk", b = 0.3, D = 15, t = 30, sigma_rw = 0.1)
# simsQ("time_varying", t = 30, beta0 = 0.1, beta1 = -0.05)
# simsQ("constant", b = 0.3)

simsDataGenQ <- function(
    # Basic model parameters
  alpha_lamb = c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6)),
  beta_lamb = 0.5,
  days = 30,
  
  # Delays: 
  # D_trunc is the truncated delay window we usually observe or assume as "final".
  # D is the full delay window considered when if_fully_reported = TRUE.
  D_trunc = NULL,
  D_complete = NULL,
  if_fully_reported = FALSE,
  
  # Reporting dynamics method: "constant", "time_varying", or "random_walk".
  # Use a named list `method_params` to hold the parameters required by the chosen method.
  method = c("constant", "time_varying", "random_walk"),
  method_params = list(
    # For "constant": list(b = 3)
    # For "time_varying": list(beta0 = -2, beta1 = 0.1)
    # For "random_walk": list(b = 3, sigma_rw = 1)
  ),
  
  # Starting date
  date_start = as.Date("2024-01-01"),
  
  # Random seed
  seed = 123
) {
  # Check inputs
  method <- match.arg(method)
  if (length(alpha_lamb) < days) {
    stop("The length of alpha_lamb cannot be less than 'days'.")
  }
  
  if (D_trunc >= D_complete) {
    stop("D_trunc must be strictly less than D_complete in the 'not fully reported' scenario.")
  }
  
  if_fully_reported <- as.logical(if_fully_reported)
  if (if_fully_reported) {
    D_used <- D_trunc
    # Make sure we REALLY converge to 1 at D_trunc => ensure_qd_one=TRUE
    ensure_qd_one_flag <- TRUE
  } else {
    D_used <- D_complete
    # We do NOT want to force qd to 1 at D_complete => ensure_qd_one=FALSE
    ensure_qd_one_flag <- FALSE
  }
  
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate the date sequence
  date_end <- date_start + days - 1
  date_seq <- as.Date(date_start:date_end)
  
  # Prepare output arrays
  lambda_t <- numeric(days)
  case_true <- matrix(numeric(days), ncol = 1)
  rownames(case_true) <- as.character(date_seq)
  
  case_reported <- matrix(0, nrow = days, ncol = D_used + 1)
  rownames(case_reported) <- as.character(date_seq)
  
  # Generate reporting probabilities qd based on the chosen method.
  # The simsQ function should return a list with:
  # $qd: cumulative probabilities (either a vector for constant or a matrix for varying methods)
  # $b_t: corresponding b(t) (if applicable)
  
  # Extract parameters according to the chosen method
  if (method == "constant") {
    # Expected parameters in method_params: b
    if (!"b" %in% names(method_params)) stop("For method='constant', method_params must include 'b'.")
    simsQ_out <- simsQ(method = "constant", 
                       b = method_params$b, 
                       D = D_used, ensure_qd_one = ensure_qd_one_flag)
  } else if (method == "time_varying") {
    # Expected parameters in method_params: beta0, beta1
    if (!all(c("beta0", "beta1") %in% names(method_params))) {
      stop("For method='time_varying', method_params must include 'beta0' and 'beta1'.")
    }
    simsQ_out <- simsQ(method = "time_varying", 
                       t = days, 
                       D = D_used, 
                       beta0 = method_params$beta0, 
                       beta1 = method_params$beta1, 
                       ensure_qd_one = ensure_qd_one_flag)
  } else {
    # method == "random_walk"
    # Expected parameters in method_params: b, sigma_rw
    if (!all(c("b", "sigma_rw") %in% names(method_params))) {
      stop("For method='random_walk', method_params must include 'b' and 'sigma_rw'.")
    }
    simsQ_out <- simsQ(method = "random_walk", 
                       t = days, 
                       D = D_used,
                       b = method_params$b, 
                       sigma_rw = method_params$sigma_rw, 
                       ensure_qd_one = ensure_qd_one_flag)
  }
  
  # Simulate data
  for (tt in seq_len(days)) {
    # 1. Draw lambda_t from Gamma distribution
    lambda_t[tt] <- rgamma(1, shape = alpha_lamb[tt], rate = beta_lamb)
    
    # 2. Simulate the true number of cases for day t
    case_true[tt] <- rpois(1, lambda = lambda_t[tt])
    
    # 3. Extract cumulative probabilities for day t
    if (method == "constant") {
      # qd is a constant vector
      prob_temp <- simsQ_out$qd
    } else {
      # qd is a matrix of size (days x (D_used+1))
      prob_temp <- simsQ_out$qd[tt, ]
    }
    
    # Convert cumulative probabilities to incremental probabilities
    # p_temp[d] = qd[d] - qd[d-1] with qd[0] = 0
    p_temp <- c(prob_temp[1], diff(prob_temp))
    
    # 4. Distribute the total cases into delay intervals using rmultinom
    reported_temp <- rmultinom(1, size = case_true[tt], prob = p_temp)
    
    # 5. Convert from increments to cumulative reported cases
    case_reported[tt, ] <- cumsum(reported_temp)
  }
  
  # If we used a full delay window (if_fully_reported = TRUE),
  # we only return the truncated part of the matrix to mimic partial observation.
  if (!if_fully_reported) {
    case_reported <- case_reported[, 1:(D_trunc + 1), drop = FALSE]
    simsQ_out$qd <- simsQ_out$qd[, 1:(D_trunc + 1), drop = FALSE]
  }
  
  
  
  # Return a list of results
  return(list(
    case_reported = case_reported,
    case_true = case_true,
    lambda_t = lambda_t,
    qd = simsQ_out$qd,
    b = simsQ_out$b_t,
    # Save which D was actually used
    D_used = D_used,
    ensure_qd_one_flag = ensure_qd_one_flag
  ))
}
# spline basis
create_basis <- function(N_obs, n_knots = 5){
  time_points <- 1:N_obs
  
  knot_dist <- 1 / (n_knots + 1)
  probs <- seq(knot_dist, 1 - knot_dist, by = knot_dist)
  knots <- quantile(time_points, probs = probs)
  # spline basis matrix
  spline_basis <- bs(time_points, knots = knots, degree = 3, intercept = TRUE)
  X_spline <- cbind(1, as.matrix(spline_basis))
  
  return(X_spline)
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

### functions to get coordinates of non-NAs

find_non_na_coords <- function(mat) {
  # dimension
  N <- nrow(mat)
  D <- ncol(mat)
  
  # indices for non NAs
  non_na_indices <- which(!is.na(mat), arr.ind = TRUE)
  
  coords_df <- as.matrix(non_na_indices)
  
  # rename cols
  colnames(coords_df) <- c("row", "col")
  
  return(coords_df)
}
