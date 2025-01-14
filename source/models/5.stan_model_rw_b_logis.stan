data {
  int<lower=0> N_obs;          // Number of time points (T)
  int<lower=0> D;              // Maximum delay
  array[N_obs, D] int<lower=0> Y; // Reported cases, dimension T x D
  int<lower=1> K;                 // Number of non-missing entries
  array[K, 2] int<lower=0> obs_index; // Coordinates (t, d) of non-missing Y
}

parameters {
  // A) Hyperparameters for lambda
  real<lower=0> alpha_lambda;         // Shape
  real<lower=0> beta_lambda;          // Rate
  // B) Poisson intensities
  vector<lower=0>[N_obs] lambda_t;    // lambda_t[t]
  // C) Random-walk std dev for transformed parameters
  real<lower=0> sigma_b_star;
  real<lower=0> sigma_phi_star;
  // D) Time-varying parameters on transformed scale
  vector[N_obs] b_star;    // Will be transformed to [0.05, 1]
  vector[N_obs] phi_star;  // Will be transformed to [0, 1]
}

transformed parameters {
  // Transform parameters back to original scale
  vector<lower=0.05,upper=1>[N_obs] b_t;
  vector<lower=0,upper=1>[N_obs] phi_t;
  matrix[N_obs, D+1] q_d_matrix;
  
  // Transform parameters back to original scale using logistic function
  for (t in 1:N_obs) {
    b_t[t] = 0.05 + 0.95 * inv_logit(b_star[t]);
    phi_t[t] = inv_logit(phi_star[t]);
  }
  
  // Compute q_d values using transformed parameters
  for (t in 1:N_obs) {
    for (d in 0:D) {
      q_d_matrix[t, (d+1)] = 1 - (1 - phi_t[t]) * exp(-b_t[t] * d);
    }
  }
}

model {
  // 1) Priors on alpha_lambda and beta_lambda (for Gamma)
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);
  
  // 2) Prior on random walk standard deviations
  sigma_b_star ~ normal(0, 2);
  sigma_phi_star ~ normal(0, 2);
  
  // 3) Priors on Poisson intensities
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // 4) Initialize first values for both parameters
  b_star[1] ~ normal(logit((0.2 - 0.05)/0.95), 0.5);
  phi_star[1] ~ normal(logit(0.2), 0.5);
  
  if (N_obs > 1) {
    b_star[2] ~ normal(b_star[1], sigma_b_star);
    phi_star[2] ~ normal(phi_star[1], sigma_phi_star);
  }
  
  // 5) Second-order random walks on transformed scale
  for (t in 3:N_obs) {
    b_star[t] ~ normal(2 * b_star[t-1] - b_star[t-2], sigma_b_star);
    phi_star[t] ~ normal(2 * phi_star[t-1] - phi_star[t-2], sigma_phi_star);
  }
  
  // 6) Likelihood
  for (k in 1:K) {
    int t = obs_index[k, 1];
    int d = obs_index[k, 2];
    Y[t, (d+1)] ~ poisson(lambda_t[t] * q_d_matrix[t, (d+1)]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
