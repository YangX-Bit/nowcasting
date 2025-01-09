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
        n_obs       = 30,            
        date_start = as.Date("2024-01-01"),
        D = 20,
        seed       = 123
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
  n_obs       <- params$data$n_obs
  date_start  <- params$data$date_start
  D           <- params$data$D
  seed        <- params$data$seed
  
  # (B) model for qd
  method          <- params$q_model$method
  method_params   <- params$q_model$method_params
  
  # other parameters
  max_delay <- 100
  
  #---------------------------------------------------------
  # 3) Check the input
  #---------------------------------------------------------
  if (!(length(alpha_lamb) == n_obs)) {
    stop("The length of `alpha_lamb` should be equal to the length of `n_obs`！")
  }
 
  
  #---------------------------------------------------------
  # 4) Set Random Seed
  #---------------------------------------------------------
  set.seed(seed)
  
  #---------------------------------------------------------
  # 5) use generateQ() to generate q(d) and b_t
  #---------------------------------------------------------
  simsQ_out <- generateQ(
    method        = method,
    method_params = method_params,
    n_obs          = n_obs,
    D             = D,
    max_delay = max_delay
  )

  #---------------------------------------------------------
  # 6) simulation
  #---------------------------------------------------------
  simulation_result <- runSimulation(
    alpha_lamb = alpha_lamb,
    beta_lamb  = beta_lamb,
    n_obs      = n_obs,
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
    n_obs,
    date_start,
    simsQ_out,
    max_delay,
    D
) {
  # Generate the date sequence
  date_seq <- seq.Date(from = date_start, by = "day", length.out = n_obs)
  
  # Initialize variables
  lambda_t     <- numeric(n_obs)        # Disease intensity
  case_true    <- integer(n_obs)        # True number of cases
  case_reported<- matrix(0, nrow = n_obs, ncol = max_delay + 1)  # Reported cases
  rownames(case_reported) <- as.character(date_seq)

  # Simulation process
  for (tt in seq_len(n_obs)){
    
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
    reported_temp <- rmultinom(n = 1, size = case_true[tt], prob = p_temp)
    
    # 6) Cumulative reported cases
    case_reported[tt, ] <- reported_temp
  }

  qd_out <- if(is.vector(simsQ_out$qd)){
    simsQ_out$qd[c(1:(D+1))]
  } else{
    simsQ_out$qd[, c(1:(D-+1))]
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
    b_t        = round(simsQ_out$b_t, 4),
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

generateQ <- function(method, method_params, n_obs, D, max_delay = 100){
  # input the method we use use, and relavant parameters
  # set number of observations, Delay D
  
  # max_delay indicates in practice, the max delay that we consider
  # since bt we set is from (0.05,1), max_col to be 100 is enough
  
  # output
  b_t_out <- NULL
  qd_out  <- NULL
  phi_out <- NULL
  
  if (method == "fixed_q"){
    if (!all(c("q_D", "q_lambda") %in% names(method_params))) {
      stop("method=constant requires b in method_params!")
    }
    q_D   <- method_params$q_D
    q_lambda <- method_params$q_lambda
    
    qd_out <- generate_exponential_q(D = q_D, lambda = q_lambda)
    b_t_out <- NA_integer_
    phi_out <- NA_integer_
    
  } else if (method == "fixed_b"){
    if (!all(c("b", "phi") %in% names(method_params))) {
      stop("method=constant requires b in method_params!")
    }
    b    <- method_params$b
    phi <- method_params$phi
    
    qd_out <- 1 - (1 - phi)*exp(-b * c(0:max_delay))
    b_t_out <- b
    phi_out <- phi
    
  } else if (method == "linear_b"){
    if (!all(c("beta0","beta1","phi_init","phi_sigma") %in% names(method_params))) {
      stop("method=time_varying requires beta0 and beta1!")
    }
    beta0 <- method_params$beta0
    beta1 <- method_params$beta1
    phi_sigma <- method_params$phi_sigma
    
    b_t_out <- exp(beta0 + beta1*(1:n_obs))
    qd_out  <- matrix(NA, nrow=n_obs, ncol=max_delay+1)
    phi_out <- numeric(n_obs); 
    phi_out[1] <- method_params$phi_init
    
    qd_out[1, ] <- 1 - (1 - phi_out[1]) * exp(-b_t_out[1] * (0:max_delay))
    for (i in c(2:n_obs)){
      phi_out[i] <- max(0, min(rnorm(1, phi_out[i-1], phi_sigma) ) )
      qd_out[i,] <- 1 - (1 - phi_out[i]) * exp(-b_t_out[i] * (0:max_delay))
    }
  } else if (method == "ou_b"){
    # OU (Ornstein-Uhlenbeck)
    #  alpha, mu, b_init, b_sigma, phi_init, phi_sigma
    if (!all(c("alpha","mu","b_init","b_sigma",
               "phi_init","phi_sigma") %in% names(method_params))) {
      stop("method=ou requires alpha, mu, b_init, sigma_ou!")
    }
    alpha   <- method_params$alpha
    mu      <- method_params$mu
    phi_sigma <- method_params$phi_sigma
    b_sigma <- method_params$b_sigma
    
    b_t_out <- numeric(n_obs); phi_out <- numeric(n_obs)
    b_t_out[1] <- method_params$b_init; phi_out[1] <- method_params$phi_init
    
    # Iterate to generate log_b_t using the OU process
    for (i in 2:n_obs){
      drift <- alpha * (mu - b_t_out[i - 1])
      b_proposal <- b_t_out[i - 1] + drift + rnorm(1, 0, b_sigma)
      phi_proposal <- rnorm(1, phi_out[i - 1], phi_sigma)
      
      b_t_out[i] <- max(min(b_proposal, 1), 0.05)
      phi_out[i] <- max(min(phi_proposal, 1), 0)
    }
    
    qd_out <- matrix(NA, nrow=n_obs, ncol=max_delay+1)
    for (i in seq_len(n_obs)){
      qd_out[i,] <- 1 - (1 - phi_out[i]) * exp(-b_t_out[i] * (0:max_delay))
    }
    
  } else {
    stop("method must be one of fixed_q, fixed_b, linear_b, ou_b!")
  }
  
  return(list(qd=qd_out, b_t=b_t_out, phi = phi_out))
}

## generate exponential decay q. A quick easy way
generate_exponential_q <- function(D, lambda = 0.3) {
  q <- exp(-lambda * (1:(D+1)))
  q <- q / sum(q)  # Normalize to sum to 1
  q_out <- cumsum(q)
  return(q_out)
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
