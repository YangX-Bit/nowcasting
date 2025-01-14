data {
  int<lower=0> N_obs;          // Number of time points (T)
  int<lower=0> D;              // Maximum delay
  array[N_obs, D] int<lower=0> Y; // Reported cases, dimension T x D
  int<lower=1> K;                 // Number of non-missing entries
  array[K, 2] int<lower=0> obs_index; // Coordinates (t, d) of non-missing Y
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

  // C) Random-walk std dev for b and phi
  real<lower=0> sigma_b;
  //real<lower=0> sigma_phi;

  // D) Time-varying b(t) and phi(t) as raw parameters
  //    We'll give them second-order RW priors below.
  //    b_t: we allow b_t[t] >= 0. Optionally clamp to [0.05, 1], etc.
  vector<lower=0>[N_obs] b_t;

  //    phi(t) is fraction in [0, 1].
  real<lower=0,upper=1> phi;
}

// -------------------------------------------------------------------
// TRANSFORMED PARAMETERS
// -------------------------------------------------------------------
transformed parameters {
  // We compute q[t] for each time t = 1..N_obs
  // Using: q[t] = 1 - (1 - phi[t]) * exp(-b_t[t] * D)
  vector[N_obs] q;

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

  // 2) Priors on random walk standard deviations
  //    You can use half-normal, half-Cauchy, exponential, etc.:
  sigma_b   ~ normal(0, 2);
  //sigma_phi ~ normal(0, 2);

  // 3) Priors on Poisson intensities
  //    Each lambda_t[t] ~ Gamma(alpha_lambda, beta_lambda)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);

  // 4) SECOND-ORDER RANDOM WALK for b_t
  //    b_t[1], b_t[2] have simpler priors; then for t>=3: second-order step
  //    Adjust center or scale to taste:
  b_t[1] ~ normal(0.2, 0.1);   
  if (N_obs > 1) {
    b_t[2] ~ normal(b_t[1], sigma_b);
  }
  for (t in 3:N_obs) {
    // b_t[t] ~ Normal( 2*b_t[t-1] - b_t[t-2], sigma_b )
    b_t[t] ~ normal(2 * b_t[t - 1] - b_t[t - 2], sigma_b);
  }

  // 5) SECOND-ORDER RANDOM WALK for phi
  //    phi in [0,1], so we do a random walk on that restricted scale.
  //    The first two get simpler priors, then second-order for t>=3.
  //    You might prefer Beta() for phi[1], or uniform(0,1):
  phi ~ uniform(0,1); 
  // if (N_obs > 1) {
  //   phi[2] ~ normal(phi[1], sigma_phi);
  // }
  // for (t in 3:N_obs) {
  //   phi[t] ~ normal(2 * phi[t - 1] - phi[t - 2], sigma_phi);
  // }

  // 6) Likelihood
  //    For each (t, d) in obs_index, Y[t,d] ~ Poisson(lambda_t[t] * q[t])
  for (k in 1:K) {
    int t = obs_index[k, 1];
    int d = obs_index[k, 2];
    Y[t, (d+1)] ~ poisson(lambda_t[t] * q[t]);
  }
}

// -------------------------------------------------------------------
// GENERATED QUANTITIES
// -------------------------------------------------------------------
generated quantities {
  // Optionally generate the latent total N_t for each time point
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    // Could also do: N_t[t] = binomial_rng( ??? ) ...
    // but typically we say: N_t ~ Poisson(lambda_t[t])
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
