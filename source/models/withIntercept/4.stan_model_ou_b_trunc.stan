data {
  int<lower=0> N_obs;                    // Number of time points (T)
  int<lower=0> D;                        // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;        // Reported cases matrix (T x D)
  int<lower=1> K;                        // Total number of observations (K)
  array[K, 2] int<lower=1> obs_index;     // Indexes for each observation (t, d)
}

parameters {
  real<lower=0> alpha_lambda;            // Shape parameter for Gamma distribution
  real<lower=0> beta_lambda;             // Rate parameter for Gamma distribution
  vector<lower=0>[N_obs] lambda_t;       // Poisson intensities (λ[t]) for each time point

  // Ornstein-Uhlenbeck process parameters
  real<lower=0> alpha;                    // Mean-reversion rate
  real mu;                                 // Long-term mean
  real<lower=0> sigma_ou;                 // Diffusion coefficient

  vector<lower=0.05, upper=1>[N_obs] b_t;                  // Logarithm of b_t parameter for stability
  //vector<lower=0>[N_obs] phi;
  real<lower=0> phi;
}

transformed parameters {
  // vector<lower=0.05, upper=1>[N_obs] b_t;
  matrix[N_obs, D] q_d_matrix;

  // Calculate b_t[t] = exp(log_b_t[t]) to ensure b_t[t] is positive
  // for (n in 1:N_obs){
  //   b_t[n] = exp(log_b_t[n]);
  // }

  // Precompute all q_d values
  for (n in 1:N_obs){
    for (d in 1:D){
      q_d_matrix[n, d] = phi + exp(- b_t[n] * d);
    }
  }
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);

  alpha ~ normal(0.5, 0.2);       // Prior for alpha
  mu    ~ normal(0.5, 0.2);       // Prior for mu
  sigma_ou ~ normal(0.1, 0.05);   // Prior for sigma_ou
  phi ~ uniform(0,0.5);

  // Gamma prior for Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);

  // Ornstein-Uhlenbeck process for b_t
  b_t[1] ~ normal(mu, sigma_ou);
  for (t in 2:N_obs){
    b_t[t] ~ normal(b_t[t-1] + alpha * (mu - b_t[t-1]), sigma_ou);
  }

  // Likelihood: Vectorized
  for (k in 1:K){
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    Y[t,d] ~ poisson(lambda_t[t] * q_d_matrix[t, d]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs){
    N_t[t] = poisson_rng(lambda_t[t]); // Generate final number of cases
  }
}
