data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  // Hyperparameters
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
    // Hyperparameters for alpha
  real mean_logit_alpha;
  real<lower=0> sd_logit_alpha;
}

transformed data {
  real<lower=0> factor = (T-1) * (2*T-1) / (6.0 * T);
}

parameters {
  vector<lower=0>[T] lambda;            // expected number of cases
  vector[T] log_b;                      // log rate of accumulated reporting probability
  vector[T] logit_phi;                  // logit of delayed reporting probability
  vector[T] logit_alpha;
  real<lower=0> sigma_log_b;
  real<lower=0> sigma_logit_phi;
  real<lower=0> sigma_logit_alpha;
}

transformed parameters {
  vector<lower=0>[T] b = exp(log_b);                     // rate of accumulated reporting probability
  vector<lower=0,upper=1>[T] phi = inv_logit(logit_phi); // delayed reporting probability
  vector<lower=0,upper=1>[T] alpha   = inv_logit(logit_alpha);
  matrix[T, D+1] q;                                      // accumulated reporting probability
  for (d in 0:D)
    for (t in 1:(T-d))
      q[t, d+1] = 1 - (1 - phi[t]) * exp(-b[t] * d);
}

model {
  // Priors
  lambda ~ lognormal(0, 2.5);
  sigma_log_b ~ lognormal(-2, 1);
  sigma_logit_phi ~ lognormal(-2, 1);
  sigma_logit_alpha ~ lognormal(-2, 1);

  // First-order random walks
  log_b[1] ~ normal(mean_log_b, sqrt(sd_log_b^2 + sigma_log_b^2 * factor));
  logit_phi[1] ~ normal(mean_logit_phi, sqrt(sd_logit_phi^2 + sigma_logit_phi^2 * factor));
  logit_alpha[1] ~ normal(mean_logit_alpha, sqrt(sd_logit_alpha^2 + sigma_logit_alpha^2 * factor));
  for (t in 2:T) {
    log_b[t] ~ normal(log_b[t-1], sigma_log_b);
    logit_phi[t] ~ normal(logit_phi[t-1], sigma_logit_phi);
    logit_alpha[t] ~ normal(logit_alpha[t-1], sigma_logit_alpha);
  }

  // Likelihood
  for (d in 0:D)
    for (t in 1:(T-d))
      Y[t, d+1] ~ poisson(lambda[t] * q[t, d+1]);
}

generated quantities {
  vector<lower=0>[T] N;                 // number of cases
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}
