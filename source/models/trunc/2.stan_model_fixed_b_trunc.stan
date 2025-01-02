data {
  int<lower=0> N_obs;         // Number of time points (T)
  int<lower=0> D;             // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;   // Reported cases (T x D matrix)
  int<lower=1> K;                   // number of obs in total
  array[K, 2] int<lower=1> obs_index; // coordinates
}

parameters {
  real<lower=0> alpha_lambda;  // Gamma prior shape parameter for lambda
  real<lower=0> beta_lambda;   // Gamma prior rate parameter for lambda
  vector<lower=0>[N_obs] lambda_t;   // Poisson intensities (λ[t]) at each time point
  real<lower=0> b;
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  b ~ uniform(0.1, 1);
  
  // Gamma prior on Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // Likelihood: Marginalized Poisson likelihood for N_t and Binomial for Y
  for(k in 1:K){
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    
    real q_d = 1 - exp(- b * d);
    Y[t,d] ~ poisson(lambda_t[t] * q_d);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]); // Sample N_t from Poisson distribution
  }
}

