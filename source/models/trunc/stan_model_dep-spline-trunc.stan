data {
  int<lower=0> N_obs;         // Number of time points (T)
  int<lower=0> D;             // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;   // Reported cases (T x D matrix)
  // bsplines
  int<lower=1> J;                          // number of basis functions
  matrix[N_obs, J] X_spline;                   // basis function values(N_obs, J)
  
  int<lower=1> K;                   // number of obs in total
  array[K, 2] int<lower=1> obs_index; // coordinates
  real<lower=0>sigma_b;
}

parameters {
  real<lower=0> alpha_lambda;  // Gamma prior shape parameter for lambda
  real<lower=0> beta_lambda;   // Gamma prior rate parameter for lambda
  vector<lower=0>[N_obs] lambda_t;   // Poisson intensities (λ[t]) at each time point
  //
  vector<lower=0>[J] beta;               // global spline coefficients (no clusters)
}

transformed parameters {
  vector[N_obs] b_t;                 // pre-calc b_t
  for (t in 1:N_obs) {
    b_t[t] = dot_product(X_spline[t], beta); // b(t) = X_spline[t] * beta
  }
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  
  beta ~ normal(0, 10);  // prior for spline coef
  
  // Gamma prior on Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  for (t in 2:N_obs) {
    b_t[t] - b_t[t - 1] ~ normal(0, sigma_b);  // constraint for b_t
  }
  

  // Likelihood: Marginalized Poisson likelihood for N_t and Binomial for Y
  for (k in 1:K) {
    int t = obs_index[k, 1];         
    int d = obs_index[k, 2];
    
    real q_d = 1 - exp(- b_t[t] * d);
    Y[t, d] ~ poisson(lambda_t[t] * q_d); //
  }
}

generated quantities {
  vector<lower=0>[N_obs] b_t_est;
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    b_t_est[t] = b_t[t] ;  // posterior of  b(t) 
    N_t[t] = poisson_rng(lambda_t[t]); // Sample N_t from Poisson distribution
  }
}


