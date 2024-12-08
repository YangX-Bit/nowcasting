library(dplyr)
library(tidyr)

### functions to generate Generalized Dirichlet distribution
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###
rGeneralizedDirichlet <- function(n = 1, alpha, beta) {
  # n is the number of samples to generate
  # alpha and beta are the parameter vectors for the Generalized Dirichlet distribution
  k <- length(alpha)  # The dimension of the vector
  X <- matrix(0, n, k)  # Initialize a matrix to store the results

  # Generate Beta-distributed random numbers sequentially according to the definition
  for (i in 1:n) {
    x <- numeric(k)  # Store the k components of each sample
    x[1] <- rbeta(1, alpha[1], beta[1])  # Generate X_1
    prod_term <- 1 - x[1]  # The remaining product term
    for (j in 2:(k - 1)) {
      x[j] <- rbeta(1, alpha[j], beta[j]) * prod_term  # Generate X_j
      prod_term <- prod_term - x[j]  # Update the remaining product term
    }
    x[k] <- prod_term  # The last component
    X[i, ] <- x  # Store the generated sample
  }

  return(X)
}


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


### functions to generate p(delay probability)
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###
simsP <- function(Type = "GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 30, seed = 123) {
  set.seed(seed)

  Type <- match.arg(Type, c("GD", "Multi_GD", "basin"))

  if (Type == "GD") {
    # fixed prob
    p <- as.numeric(rGeneralizedDirichlet(1, gd_alpha, gd_beta))

  } else if (Type == "Multi_GD") {
    # multi
    p <- rGeneralizedDirichlet(days, gd_alpha, gd_beta)

  } else {  # Type == "basin"
    # spit into three parts
    n <- split_into_three(days)
    p <- matrix(0, nrow = days, ncol = D + 1)

    for (i in 1:days) {
      if (i <= n[1]) {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 1, beta = 1000, order = "decreasing")

      } else if (i <= (n[1] + n[2])) {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 10, beta = 1)

      } else {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 1, beta = 1000, order = "increasing")
      }
    }
  }

  return(p)
}
# 

simsQ <- function(method = "constant", 
                  b = 3, D = 15, t = 30,
                  beta0 = -2, beta1 = 0.1, sigma_rw = 0.05) {
  # 'method' defines how Q is simulated
  # 'b' can be a constant or a time-varying function
  # 'd' is the delay time, and 't' is the current time point (optional)
  # 'sigma_rw' is the standard deviation of the random walk increment
  
  if (method == "constant") {
    # Q with a constant b
    qd <- 1 - exp(-b * 1:(D + 1))
    
    return(list(qd = qd, b_t = b))
    
  } else if (method == "time_varying") {
    # Time-varying Q, b_t increases over time
    if (is.null(t)) stop("Time 't' must be provided for time-varying method.")
    
    qd <- matrix(NA, nrow = t, ncol = D + 1)
    b_t <-  exp(beta0 + beta1 * 1:t)   
    
    for (i in 1:t) {
      qd[i,] <- 1 - exp(-b_t[i] * 1:(D + 1))
    }
    
    return(list(qd = qd, b_t = b_t)) 
    
  } else if (method == "random_walk") {
    # Random walk Q, b_t fluctuates over time
    if (is.null(t)) stop("Time 't' must be provided for random_walk method.")
    qd <- matrix(NA, nrow = t, ncol = D + 1)  # Initialize matrix for each time step
    b_t <- numeric(t)                         # Vector to store b_t values
    b_t[1] <- b  # Initialize the first b_t value
    
    for (i in 2:t) {
      b_t[i] <- exp(log(b_t[i - 1]) + rnorm(1, mean = 0, sd = sigma_rw))  # Random walk step
    }
    
    # Compute q(d) for each delay at each time step
    for (i in 1:t) {
      qd[i, ] <- 1 - exp(-b_t[i] * 1:(D + 1))
    }
    
    return(list(qd = qd, b_t = b_t))  # Return both Q and the time-varying b_t
    
  } else {
    stop("Unknown method for Q simulation.")
  }
}

set.seed(1)
simsQ(method = "random_walk", b = 0.3, D = 15, t = 30, sigma_rw = 0.3)
simsQ("time_varying", t = 30)
simsQ("constant", b = 3)

simsDataGenP <- function(alpha =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta = 0.5,
                         p = NULL,
                         days = length(alpha), D = 15, seed = 123){
  if(length(alpha) < days){
    stop("Error! The length of alpha cannot be less than days!")
  }
  
  set.seed(seed)
  
  # Cut the probability according to the Maximal delay that we concern
  p_cut <- ifelse(is.null(nrow(p)), p[c(1:(D+1))], p[, c(1:(D+1))])
  
  # Arrays to accumulate the results across multiple simulations
  lambda_t <- numeric(days)  # Lambda for each day 
  case_true <- numeric(days) # Actual number of cases per day

  # Simulate the true number of cases per day
  for (t in 1:days) {
    # Draw the Poisson intensity parameter lambda_t from a Gamma distribution
    lambda_t[t] <- rgamma(1, shape = alpha[t], rate = beta)
    
    # Draw the actual number of cases N(t, ∞) from a Poisson distribution
    case_true[t] <- rpois(1, lambda = lambda_t[t])
  }
  
  case_reported <- matrix(0, nrow = days, ncol = D + 1)
  # Reported cases based on delay distribution
  if (is.null(nrow(p_cut))) {
    p_cut <- p[c(1:(D+1))] # cut by maximal delay that we concern
    for (i in 1:days) {
      case_reported[i,] = rmultinom(1, size = case_true[i], prob = p_cut)
    }
  } else {
    p_cut <-  p[, c(1:(D+1))]
    for (i in 1:days) {
      case_reported[i,] = rmultinom(1, size = case_true[i], prob = p_cut[i, ])
    }
  }
  out <- list(case_reported = t(apply(case_reported, 1, cumsum)), 
              case_true = case_true, lambda_t = lambda_t,
              case_reported_non_accum = case_reported
              )
  
  return(out)
}

simsDataGenQ <- function(alpha_lamb =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta_lamb = 0.5,
                         b = 3, method = "constant",
                         beta0 = -2, beta1 = 0.1,
                         sigma_rw = 1,
                         date_start = as.Date("2024-01-01"), days = 30, D = 15, seed = 123){
  if(length(alpha_lamb) < days){
    stop("Error! The length of alpha cannot be less than days!")
  }
  
  # date
  date_end = date_start + days - 1
  date = as.Date(date_start:date_end)
  
  set.seed(seed)
  
  # Arrays to accumulate the results
  lambda_t <- numeric(days)  # Lambda for each day 
  case_true <- numeric(days) # Actual number of cases per day
  Q <- matrix(0, nrow = days, ncol = D + 1)
  case_reported <-  matrix(0, nrow = days, ncol = D + 1)
  
  # report proportion
  if(method == "constant"){
    simsQ_out <- simsQ("constant", b = b, D = D)
    prob_temp <- simsQ_out$qd
  }else if(method == "time_varying"){
    simsQ_out <- simsQ("time_varying", D = D, t = days, beta0 = beta0, beta1 = beta1)
  }else{
    simsQ_out <- simsQ(method = "random_walk", b = b, D = D, t = days, sigma_rw = sigma_rw)
  }
  
  # Simulate the true number of cases per day
  for (t in 1:days) {
    # Draw the Poisson intensity parameter lambda_t from a Gamma distribution
    lambda_t[t] <- rgamma(1, shape = alpha_lamb[t], rate = beta_lamb)
    
    # Draw the actual number of cases from a Poisson distribution
    case_true[t] <- rpois(1, lambda = lambda_t[t])
    
    # D times from binom
    reported_temp <- numeric(D+1)
    
    if(method != "constant"){ prob_temp <- simsQ_out$qd[t,] }
    
    for(d in 1:(D+1)){
      reported_temp[d] <- rbinom(1, size = case_true[t], prob = prob_temp[d])
    }
    case_reported[t, ] <- sort(reported_temp)
  }
  rownames(case_reported) <- as.character(date)

  out <- list(case_reported = case_reported, case_true = case_true,
              lambda_t = lambda_t, qd = simsQ_out$qd, b = simsQ_out$b_t)
  return(out)
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
