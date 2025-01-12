data {
  int<lower=0> N_obs;          // Number of time points (T)
  int<lower=0> D;              // Maximum delay
  array[N_obs, D] int<lower=0> Y; // Reported cases, dimension T x D
  int<lower=1> K;                 // Number of non-missing entries
  array[K, 2] int<lower=1> obs_index; // Coordinates (t, d) of non-missing Y
}

// -------------------------------------------------------------------
// PARAMETERS
// -------------------------------------------------------------------
parameters {
  // A) Hyperparameters for lambda
  real<lower=0> alpha_lambda;         // Shape
  real<lower=0> beta_lambda;          // Rate
  // B) Poisson intensities
  vector<lower=0>[N_obs] lambda_t;    // lambda_t[t]
  // C) Random-walk std dev for transformed b
  real<lower=0> sigma_b_star;
  // D) Time-varying b(t) on transformed scale
  //    Now b_star is unconstrained (real line) and will be transformed to [0.05, 1]
  vector[N_obs] b_star;
  // E) phi remains on [0,1] scale
  real<lower=0,upper=1> phi;
}

// -------------------------------------------------------------------
// TRANSFORMED PARAMETERS
// -------------------------------------------------------------------
transformed parameters {
  // Transform b_star back to original scale [0.05, 1]
  vector<lower=0.05,upper=1>[N_obs] b_t;
  // q[t] computed using transformed b_t
  vector[N_obs] q;
  
  // Logistic transform from real line to [0.05, 1]
  for (t in 1:N_obs) {
    b_t[t] = 0.05 + 0.95 * inv_logit(b_star[t]);
  }
  
  // Compute q using transformed b_t
  for (t in 1:N_obs) {
    q[t] = 1 - (1 - phi) * exp(-b_t[t] * D);
  }
}

// -------------------------------------------------------------------
// MODEL
// -------------------------------------------------------------------
model {
  // 1) Priors on alpha_lambda and beta_lambda (for Gamma)
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);
  
  // 2) Prior on random walk standard deviation (now for transformed scale)
  sigma_b_star ~ normal(0, 2);
  
  // 3) Priors on Poisson intensities
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // 4) SECOND-ORDER RANDOM WALK for b_star (on transformed scale)
  //    Initialize with values that will transform to reasonable b_t values
  b_star[1] ~ normal(logit((0.2 - 0.05)/0.95), 0.5);
  
  if (N_obs > 1) {
    b_star[2] ~ normal(b_star[1], sigma_b_star);
  }
  
  // Random walk on transformed scale
  for (t in 3:N_obs) {
    b_star[t] ~ normal(2 * b_star[t-1] - b_star[t-2], sigma_b_star);
  }
  
  // 5) Prior for phi (kept as is)
  phi ~ uniform(0, 1);
  
  // 6) Likelihood
  for (k in 1:K) {
    int t = obs_index[k, 1];
    int d = obs_index[k, 2];
    Y[t, d] ~ poisson(lambda_t[t] * q[t]);
  }
}

// -------------------------------------------------------------------
// GENERATED QUANTITIES
// -------------------------------------------------------------------
generated quantities {
  vector<lower=0>[N_obs] N_t;
  
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
