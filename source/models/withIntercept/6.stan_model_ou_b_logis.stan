data {
  int<lower=0> N_obs;                    // Number of time points (T)
  int<lower=0> D;                        // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;        // Reported cases matrix (T x D)
  int<lower=1> K;                        // Total number of observations (K)
  array[K, 2] int<lower=1> obs_index;    // Indexes for each observation (t, d)
}

parameters {
  real<lower=0> alpha_lambda;            // Shape parameter for Gamma distribution
  real<lower=0> beta_lambda;             // Rate parameter for Gamma distribution
  vector<lower=0>[N_obs] lambda_t;       // Poisson intensities (λ[t]) for each time point
  
  // Ornstein-Uhlenbeck process parameters for transformed space
  real<lower=0> alpha_b;                 // Mean-reversion rate for b
  real<lower=0.05,upper=1> mu_b;         // Long-term mean for b (original scale)
  real<lower=0> sigma_b_star;            // Diffusion coefficient for transformed b
  
  real<lower=0> alpha_phi;               // Mean-reversion rate for phi
  real<lower=0,upper=1> mu_phi;          // Long-term mean for phi (original scale)
  real<lower=0> sigma_phi_star;          // Diffusion coefficient for transformed phi
  
  // Parameters on transformed (unconstrained) scale
  vector[N_obs] b_star;                  // Transformed b parameter
  vector[N_obs] phi_star;                // Transformed phi parameter
}

transformed parameters {
  // Transform parameters back to original scale
  vector<lower=0.05,upper=1>[N_obs] b_t;
  vector<lower=0,upper=1>[N_obs] phi;
  matrix[N_obs, D] q_d_matrix;
  
  // Transform mu to logistic scale for the OU process
  real mu_b_star = logit((mu_b - 0.05)/0.95);    // Target for b_star
  real mu_phi_star = logit(mu_phi);              // Target for phi_star
  
  // Transform parameters back to original scale using logistic function
  for (n in 1:N_obs) {
    b_t[n] = 0.05 + 0.95 * inv_logit(b_star[n]);
    phi[n] = inv_logit(phi_star[n]);
  }
  
  // Compute q_d values using transformed parameters
  for (n in 1:N_obs) {
    for (d in 1:D) {
      q_d_matrix[n, d] = 1 - (1 - phi[n]) * exp(-b_t[n] * d);
    }
  }
}

model {
  // Priors for Poisson process parameters
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  
  // Priors for OU parameters
  alpha_b ~ normal(0.5, 0.2);
  mu_b ~ normal(0.5, 0.2);
  sigma_b_star ~ normal(0.1, 0.05);
  
  alpha_phi ~ normal(0.5, 0.2);
  mu_phi ~ uniform(0, 1);
  sigma_phi_star ~ normal(0.1, 0.05);
  
  // Gamma prior for Poisson intensities
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // Ornstein-Uhlenbeck process on transformed scale
  // Initial values
  b_star[1] ~ normal(mu_b_star, sigma_b_star);
  phi_star[1] ~ normal(mu_phi_star, sigma_phi_star);
  
  // OU process evolution
  for (t in 2:N_obs) {
    // Mean reversion to transformed targets
    b_star[t] ~ normal(b_star[t-1] + alpha_b * (mu_b_star - b_star[t-1]), 
                      sigma_b_star);
    phi_star[t] ~ normal(phi_star[t-1] + alpha_phi * (mu_phi_star - phi_star[t-1]), 
                        sigma_phi_star);
  }
  
  // Likelihood
  for (k in 1:K) {
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    Y[t,d] ~ poisson(lambda_t[t] * q_d_matrix[t, d]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
