data {
  int<lower=0> N_obs;          // Number of time points (T)
  int<lower=0> D;              // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y; // Reported cases matrix (T x D)
  int<lower=1> K;              // Total number of observations (K)
  array[K, 2] int<lower=1> obs_index; // Indexes for each observation (t, d)
}

parameters {
  // --- Poisson parameters ---
  real<lower=0> alpha_lambda;  
  real<lower=0> beta_lambda;   
  vector<lower=0>[N_obs] lambda_t;  

  // --- OU process for b_t ---
  real<lower=0> alpha_b;  
  real mu_b;             
  real<lower=0> sigma_b; 
  vector<lower=0.05, upper=1>[N_obs] b_t;

  // --- Hierarchical Normal for phi ---
  //real mu_phi;                // global mean of phi
  //real<lower=0> sigma_phi;     // global SD controlling phi dispersion
  vector<lower=0, upper=1>[N_obs] phi; // each phi[t]
}

transformed parameters {
  matrix[N_obs, D] q_d_matrix;

  // Precompute q_d_matrix
  for (n in 1:N_obs) {
    for (d in 1:D) {
      q_d_matrix[n, d] = 1 - (1 - phi[n]) * exp(- b_t[n] * d);
    }
  }
}

model {
  // --- Priors on lambda ---
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);

  // --- Prior for OU process b_t ---
  alpha_b ~ normal(0.5, 0.2);
  mu_b    ~ normal(0.5, 0.2);
  sigma_b ~ normal(0.1, 0.05);

  // --- Hierarchical Normal priors for phi ---
  //mu_phi   ~ uniform(0, 1);      // or normal(0.5, 0.2), etc.
  //sigma_phi ~ normal(0.1, 0.05);  // or some other reasonable prior

  // --- Gamma prior for Poisson intensities ---
  lambda_t ~ gamma(alpha_lambda, beta_lambda);

  // --- OU process for b_t ---
  b_t[1] ~ normal(mu_b, sigma_b);
  for (t in 2:N_obs) {
    b_t[t] ~ normal(b_t[t-1] + alpha_b * (mu_b - b_t[t-1]), sigma_b);
  }

  // --- Hierarchical Normal for phi ---
  //     (truncated automatically by <lower=0, upper=1>)
  for (t in 1:N_obs) {
    phi[t] ~ normal(0.1, 0.05);
  }

  // --- Likelihood ---
  for (k in 1:K) {
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    Y[t, d] ~ poisson(lambda_t[t] * q_d_matrix[t, d]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]); 
  }
}
