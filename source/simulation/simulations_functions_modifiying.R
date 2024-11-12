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

# sum(rGeneralizedDirichlet(n = 1, alpha = seq(3,3.3, by = 0.01), beta = rep(30,30))[1:15])


### functions to generate p(delay probability)
# Parameters:
#  n - number of samples
#  alpha - to control the rbeta(). If alpha is large, than the value is large.
#  beta - to control the rbeta(). If beta is large, than the value is large.
###
# generate_controlled_sum <- function(n, alpha, beta, order = "increasing") {
#   # n random numbers
#   values <- rbeta(n, alpha, beta)
#   
#   # normalize to  1
#   normalized_values <- values / sum(values)
#   
#   # Arrange order
#   if (order == "increasing") {
#     normalized_values <- sort(normalized_values)
#   } else if (order == "decreasing") {
#     normalized_values <- sort(normalized_values, decreasing = TRUE)
#   } 
#   
#   return(normalized_values)
# }

### functions to split n into three shares
# Parameters:
#  n - number of samples
###
# split_into_three <- function(n) {
#   base <- n %/% 3
#   remainder <- n %% 3
#   return(c(rep(base + 1, remainder), rep(base, 3 - remainder)))
# }

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
# p1 <- simsP(Type = "GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 58)
# p2 <- simsP(Type = "Multi_GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 58)
# p3 <- simsP(Type = "basin", D = 15, days = 58)



simsQ <- function(method = "constant", 
                  b = 3, D = 15, t = NULL,
                  beta0 = -2, beta1 = 0.1){
  # 'method' defines how Q is simulated
  # 'b' can be a constant or a time-varying function
  # 'd' is the delay time, and 't' is the current time point (optional)
  
  if (method == "constant") {
    qd <- 1 - exp(-b * 1:(D + 1))  # Q with constant b
  } else if (method == "time_varying") {
    if (is.null(t)) stop("Time 't' must be provided for time-varying method.")
    for (i in 1:t) {
      b_t <- exp(beta0 + beta1 * i)  
      
      qd <- 1 - exp(-b_t * 1:(D + 1))
    }
  } else {
    stop("Unknown method for Q simulation.")
  }
  
  return(qd)
}

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


# par(mfrow = c(5,6))
# for (i in 1:30) {
#   plot(p3[i,], type = "l",
#        xlab = "Delays", ylab = "")
# }
# #
# for (i in 1:30) {
#   plot(t(apply(p3, 1, cumsum))[i,], type = "l", 
#        xlab = "Delays", ylab = "")
# }

### functions to generate simulation data
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###

# old version with multiple sims
# simsDataGenP <- function(alpha =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta = 0.5,
#                          p = NULL,
#                         days = 30, D = 15, seed = 123, simulations = 1){
#   if(length(alpha) < days){
#     stop("Error! The length of alpha cannot be less than days!")
#   }
# 
#   set.seed(seed)
# 
#   # Cut the probability according to the Maximal delay that we concern
#   p_cut <- ifelse(is.null(nrow(p)), p[c(1:(D+1))], p[, c(1:(D+1))])
# 
# 
#   # Arrays to accumulate the results across multiple simulations
#   lambda_t <- matrix(0, nrow = days, ncol = simulations)  # Lambda for each day and each simulation
#   true_cases <- matrix(0, nrow = days, ncol = simulations)  # Actual number of cases per day for each simulation
#   reported_cases <- list()
# 
#   # Loop over the number of simulations
#   for (sim in 1:simulations) {
#     # Simulate the true number of cases per day
#     for (t in 1:days) {
#       # Draw the Poisson intensity parameter lambda_t from a Gamma distribution
#       lambda_t[t, sim] <- rgamma(1, shape = alpha[t], rate = beta)
# 
#       # Draw the actual number of cases N(t, ∞) from a Poisson distribution
#       true_cases[t, sim] <- rpois(1, lambda = lambda_t[t, sim])
#     }
# 
#     reported_cases_temp <- matrix(0, nrow = days, ncol = D + 1)
#     # Reported cases based on delay distribution
#     if (is.null(nrow(p_cut))) {
#       p_cut <- p[c(1:(D+1))] # cut by maximal delay that we concern
#       for (i in 1:days) {
#         reported_cases_temp[i, ] = rmultinom(1, size = true_cases[i, sim], prob = p_cut)
#       }
#     } else {
#       p_cut <-  p[, c(1:(D+1))]
#       for (i in 1:days) {
#         reported_cases_temp[i, ] = rmultinom(1, size = true_cases[i, sim], prob = p_cut[i, ])
#       }
#     }
#     reported_cases[[sim]] <- reported_cases_temp
#   }
#   # Take the average of lambda_t, true_cases, and reported_cases across simulations
#   lambda_t <- apply(lambda_t, 1, mean)
#   true_cases <- apply(true_cases, 1, mean)
#   reported_cases <- Reduce("+", reported_cases) / simulations
# 
#   out <- list(case_reported = t(apply(reported_cases, 1, cumsum)), 
#               case_true = true_cases, lambda = lambda_t,
#               case_reported_non_accum =reported_cases)
# 
#   return(out)
# }

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

# old version with multiple sims
# simsDataGenQ <- function(alpha_lamb =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta_lamb = 0.5,
#                         b = 3, method = "constant",
#                         beta0 = -2, beta1 = 0.1,
#                         days = 30, D = 15, seed = 123, simulations = 1){
#   if(length(alpha_lamb) < days){
#     stop("Error! The length of alpha cannot be less than days!")
#   }
# 
#   set.seed(seed)
# 
#   # Arrays to accumulate the results across multiple simulations
#   lambda_t <- matrix(0, nrow = days, ncol = simulations)  # Lambda for each day and each simulation
#   case_true <- matrix(0, nrow = days, ncol = simulations)  # Actual number of cases per day for each simulation
#   case_reported <- array(0, dim = c(days, D + 1, simulations))  # Reported cases matrix for delays 0 to D days for each simulation
# 
#   # Loop over the number of simulations
#   for (sim in 1:simulations) {
#     # Simulate the true number of cases per day
#     for (t in 1:days) {
#       # Draw the Poisson intensity parameter lambda_t from a Gamma distribution
#       lambda_t[t, sim] <- rgamma(1, shape = alpha_lamb[t], rate = beta_lamb)
# 
#       # Draw the actual number of cases from a Poisson distribution
#       case_true[t, sim] <- rpois(1, lambda = lambda_t[t, sim])
# 
#       # D times from binom
#       reported_temp <- numeric(D+1)
# 
#       # report proportion
#       if(method == "constant"){
#         prob_temp <- simsQ("constant", b = b)
#       }else if(method == "time_varying"){
#         prob_temp <- simsQ("time_varying", t = t, beta0 = beta0, beta1 = beta1)
#       }
# 
#       for(d in 1:(D+1)){
#         reported_temp[d] <- rbinom(1, size = case_true[t, sim], prob = prob_temp[d])
#       }
#       case_reported[t, , sim] <- sort(reported_temp)
#       #reported_cases[t, , sim] <- reported_temp
#     }
#   }
#   avg_reported_cases <- apply(case_reported, c(1, 2), mean)
#   out <- list(case_reported = avg_reported_cases, case_true = apply(case_true, 1, mean),
#               lambda = apply(lambda_t, 1, mean))
#   return(out)
# }

simsDataGenQ <- function(alpha_lamb =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta_lamb = 0.5,
                         b = 3, method = "constant",
                         beta0 = -2, beta1 = 0.1,
                         sigma_rw = 1,
                         date_start = as.Date("2024-01-01"),days = 30, D = 15, seed = 123){
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

  out <- list(case_reported = case_reported, case_true = case_true,
              lambda_t = lambda_t, qd = simsQ_out$qd, b = simsQ_out$b_t,
              date = date)
  return(out)
}




### functions to transfer the simulation data to the form of data in the paper
# Parameters:
#  data - Matrix of cases in each day with delays
###
dataTransform <- function(data, 
                          start_date = as.Date("2011-07-04")){
  
  # sequence of the start date
  admission_dates <- start_date - 0:(nrow(data) - 1)
  
  
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
  
  # 
  for (i in 1:(N-D)) {
    # keeps the full data
    Y_triangular[i, ] <- Y_full[i, ]
    # Y_triangular[N-i+1, 1:i] <- Y_full[N-i+1, 1:i]
  }
  
  for (j in 1:D) {
    # keeps
    Y_triangular[N-j+1, 1:j] <- Y_full[N-j+1, 1:j]
  }
  
  if(if_zero){
    Y_triangular[is.na(Y_triangular)] <- 0  
  }
  
  return(Y_triangular)
}

### functions to create basis

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


### functions to get coordinates of non-NAs
# find_non_na <- function(mat) {
#   # coordinates
#   non_na_indices <- which(!is.na(mat), arr.ind = TRUE)
#   
#   coords_df <- as.data.frame(non_na_indices)
#   
#   # extract the value for each element
#   coords_df$value <- mat[cbind(non_na_indices[, 1], non_na_indices[, 2])]
#   
#   # rename the col
#   colnames(coords_df) <- c("row", "col", "value")
#   
#   return(coords_df)
# }


find_non_na_coords <- function(mat) {
  # dimension
  N <- nrow(mat)
  D <- ncol(mat)

  # indices for non NAs
  non_na_indices <- which(!is.na(mat), arr.ind = TRUE)

  coords_df <- as.data.frame(non_na_indices)

  # rename cols
  colnames(coords_df) <- c("row", "col")

  return(coords_df)
}


