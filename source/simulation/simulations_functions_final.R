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
        days       = 30,            
        date_start = as.Date("2024-01-01"),
        seed       = 123
      ),
      # B. Delay report structure
      reporting = list(
        D_trunc          = 5,
        D_complete       = 10,
        if_fully_reported= FALSE   # FALSE => NFR; TRUE => FR
      ),
      # C. Model selection
      q_model = list(
        method       = "constant",  # "constant","time_varying","random_walk","ou"
        # different params for diff methods
        method_params= list(
          # constant: b
          # time_varying: beta0,beta1
          # random_walk: b_init, sigma_rw
          # ou: alpha, mu, b_init, sigma_ou
          b = 0.5
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
  days        <- params$data$days
  date_start  <- params$data$date_start
  seed        <- params$data$seed
  
  # (B) structure
  D_trunc          <- params$reporting$D_trunc
  D_complete       <- params$reporting$D_complete
  if_fully_reported<- params$reporting$if_fully_reported
  
  # (C) model for qd
  method          <- params$q_model$method
  method_params   <- params$q_model$method_params
  
  #---------------------------------------------------------
  # 3) Check the input
  #---------------------------------------------------------
  if (length(alpha_lamb) < days) {
    stop("The length of `alpha_lamb` should be equal to the length of `days`！")
  }
  
  if (D_trunc >= D_complete) {
    stop("Ensure `D_trunc < D_complete` (for NFR scenario). For Fully Reported Case, set if_fully_reported=TRUE.")
  }
  
  if_fully_reported <- as.logical(if_fully_reported)
  # In FR scenario, use D_trunc to generate Q(d) and ensure q_d(D_trunc) = 1
  # In NFR scenario, use D_complete to generate Q(d) without forcing q_d to 1
  if (if_fully_reported) {
    D_used            <- D_trunc
    ensure_qd_one <- TRUE
  } else {
    D_used            <- D_complete
    ensure_qd_one <- FALSE
  }
  
  #---------------------------------------------------------
  # 4) Set Random Seed
  #---------------------------------------------------------
  set.seed(seed)
  
  #---------------------------------------------------------
  # 5) Internal Function to Generate b_t and q(d)
  #---------------------------------------------------------
  makeQ <- function(method, days, D, ensure_qd, method_params){
    
    # normalize_qd => Force the D+1 point to be 1
    normalize_qd <- function(qd_mat){
      if (is.vector(qd_mat)) {
        fval <- qd_mat[D+1]
        if (fval > 0) qd_mat <- qd_mat / fval
      } else {
        # matrix
        fvals <- qd_mat[, D+1]
        valid <- fvals > 0
        qd_mat[valid, ] <- qd_mat[valid, ] / fvals[valid]
      }
      return(qd_mat)
    }
    
    # outputs
    b_t_out <- NULL
    qd_out  <- NULL
    
    if (method == "constant"){
      if (!("b" %in% names(method_params))) {
        stop("method=constant requires b in method_params!")
      }
      b0    <- method_params$b
      qd_out<- 1 - exp(-b0 * (1:(D+1)))
      b_t_out <- b0
      
    } else if (method == "time_varying"){
      if (!all(c("beta0","beta1") %in% names(method_params))) {
        stop("method=time_varying requires beta0 and beta1!")
      }
      beta0 <- method_params$beta0
      beta1 <- method_params$beta1
      
      b_t_out <- exp(beta0 + beta1*(1:days))
      qd_out  <- matrix(NA, nrow=days, ncol=D+1)
      for (i in seq_len(days)){
        qd_out[i,] <- 1 - exp(-b_t_out[i] * (1:(D+1)))
      }
      
    } else if (method == "random_walk"){
      if (!all(c("b","sigma_rw") %in% names(method_params))) {
        stop("method=random_walk requires b_init and sigma_rw!")
      }
      b_init   <- method_params$b_init
      sigma_rw <- method_params$sigma_rw
      
      b_t_out  <- numeric(days)
      b_t_out[1] <- b_init
      for (i in 2:days){
        proposal   <- b_t_out[i-1] + rnorm(1,0,sigma_rw)
        b_t_out[i] <- max(min(proposal, 1), 0.05)
      }
      
      qd_out <- matrix(NA, nrow=days, ncol=D+1)
      for (i in seq_len(days)){
        qd_out[i,] <- 1 - exp(-b_t_out[i] * (1:(D+1)))
      }
      
    } else if (method == "ou"){
      # OU (Ornstein-Uhlenbeck)
      #  alpha, mu, b_init, sigma_ou
      if (!all(c("alpha","mu","b_init","sigma_ou") %in% names(method_params))) {
        stop("method=ou requires alpha, mu, b_init, sigma_ou!")
      }
      alpha   <- method_params$alpha
      mu      <- method_params$mu
      b_init  <- method_params$b_init
      sigma_o <- method_params$sigma_ou
      
      b_t_out <- numeric(days)
      b_t_out[1] <- b_init
      
      # Iterate to generate log_b_t using the OU process
      for (i in 2:days){
        drift <- alpha * (mu - b_t_out[i - 1])
        proposal <- b_t_out[i - 1] + drift + rnorm(1, 0, sigma_o)
        
        # Ensure log_b_t does not go below log(0.05)
        b_t_out[i] <- max(min(proposal, 1), 0.05)
      }
      
      # log_b_t_out <- numeric(days)
      # log_b_t_out[1] <- log(b_init)
      # 
      # # Iterate to generate log_b_t using the OU process
      # for (i in 2:days){
      #   drift <- alpha * (log(mu) - log_b_t_out[i - 1])
      #   proposal <- log_b_t_out[i - 1] + drift + rnorm(1, 0, log(sigma_o))
      #   
      #   # Ensure log_b_t does not go below log(0.05)
      #   log_b_t_out[i] <- max(proposal, log(0.05))
      # }
      # # Transform back to the original scale
      # b_t_out <- exp(log_b_t_out)
      
      qd_out <- matrix(NA, nrow=days, ncol=D+1)
      for (i in seq_len(days)){
        qd_out[i,] <- 1 - exp(-b_t_out[i] * (1:(D+1)))
      }
      
    } else {
      stop("method must be one of constant, time_varying, random_walk, ou!")
    }
    
    # if ensure_qd = T, then do normalization
    if (ensure_qd) {
      qd_out <- normalize_qd(qd_out)
    }
    
    return(list(qd=qd_out, b_t=b_t_out))
  }
  
  #---------------------------------------------------------
  # 6) use makeQ() to generate q(d) and b_t
  #---------------------------------------------------------
  simsQ_out <- makeQ(
    method        = method,
    days          = days,
    D             = D_used,
    ensure_qd     = ensure_qd_one,
    method_params = method_params
  )
  
  #---------------------------------------------------------
  # 7) start simulation
  #---------------------------------------------------------
  # generate the date sequence
  date_seq <- seq.Date(from=date_start, by="day", length.out = days)
  
  #
  lambda_t     <- numeric(days)        # disease intensity
  case_true    <- integer(days)        
  case_reported<- matrix(0, nrow=days, ncol=D_used+1)
  rownames(case_reported) <- as.character(date_seq)
  
  for (tt in seq_len(days)){
    
    # 1) λ_t ~ Gamma
    lambda_t[tt] <- rgamma(1, shape=alpha_lamb[tt], rate=beta_lamb)
    
    # 2) true num of cases
    case_true[tt] <- rpois(1, lambda=lambda_t[tt])
    
    # 3) Cumulated report proportion of a certain day
    if (is.vector(simsQ_out$qd)){
      # constant -> simsQ_out$qd 是一个 vector
      prob_temp <- simsQ_out$qd
    } else {
      # else -> simsQ_out$qd is a (days x (D+1)) matrix
      prob_temp <- simsQ_out$qd[tt, ]
    }
    
    # 4) cumu prop -> single prop
    p_temp <- c(prob_temp[1], diff(prob_temp))
    
    # 5) distribute case_true[tt] to each delay d = 0,1,2,...
    reported_temp <- rmultinom(n=1, size=case_true[tt], prob=p_temp)
    
    # 6) cumulated report case
    case_reported[tt, ] <- cumsum(reported_temp)
  }
  
  case_true = as.matrix(case_true)
  rownames(case_true) = as.character(date_seq)
  
  #---------------------------------------------------------
  # 8) if_fully_reported = FALSE (i.e. NFR), we only keep first D_trunc+1 cols
  #---------------------------------------------------------
  if (!if_fully_reported){
    case_reported <- case_reported[, 1:(D_trunc+1), drop=FALSE]
    # cut them by D_trunc
    #  constant => qd is vector; or qd is matrix
    if (is.vector(simsQ_out$qd)){
      simsQ_out$qd <- simsQ_out$qd[1:(D_trunc+1)]
    } else {
      simsQ_out$qd <- simsQ_out$qd[, 1:(D_trunc+1), drop=FALSE]
    }
  }
  
  #---------------------------------------------------------
  # 9) output
  #---------------------------------------------------------
  return(list(
    # reported cases
    case_reported = case_reported,
    # true cases
    case_true  = case_true,
    lambda_t   = round(lambda_t),
    # b(t) & q(d)
    b_t        = round(simsQ_out$b_t, 4),
    qd         = round(simsQ_out$qd, 4),
    # other inf
    date_seq   = date_seq,
    D_used     = D_used,
    # lable for FR/NFR
    if_fully_reported = if_fully_reported
  ))
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
