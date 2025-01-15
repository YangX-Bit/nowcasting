#################################
### beta prior ###

calculate_beta_params <- function(mu, sigma2) {
  if (mu <= 0 || mu >= 1) {
    stop("The mean (mu) must be in the range (0, 1).")
  }
  if (sigma2 <= 0 || sigma2 >= mu * (1 - mu)) {
    stop("The variance (sigma2) must be > 0 and < mu * (1 - mu).")
  }
  
  # Calculate alpha + beta
  alpha_plus_beta <- mu * (1 - mu) / sigma2 - 1
  
  # Calculate alpha and beta
  alpha <- mu * alpha_plus_beta
  beta <- (1 - mu) * alpha_plus_beta
  
  return(list(alpha = alpha, beta = beta))
}

beta_mu_beta <- function(a,b){
  mu <- a/(a+b)
  sigma <- a*b/((a+b)^2 * (a+b+1))
  return(list(mu,sigma))
}

# scenario 1
mu_1 = 0.3
sigma2_1 = 0.1

scenario_1 <- calculate_beta_params(mu_1, sigma2_1)

beta_mu_beta(scenario_1$alpha, scenario_1$beta)

# scenario 2
mu_2 = 0.15
sigma2_2 = 0.05

scenario_2 <- calculate_beta_params(mu_2, sigma2_2)

beta_mu_beta(scenario_2$alpha, scenario_2$beta)



beta_mu_beta(2, 1)



#################################
### plots for exp_grow func ###

