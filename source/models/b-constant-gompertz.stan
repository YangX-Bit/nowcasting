data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  // Hyperparameters
  real<lower=0> alpha_phi;
  real<lower=0> beta_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
}

parameters {
  vector<lower=0>[T] lambda;            // expected number of cases
  real<lower=0> b;                      // rate of accumulated reporting probability
  real<lower=0, upper=1> phi;           // delayed reporting probability
}

transformed parameters {
  vector<lower=0, upper=1>[D+1] q;      // accumulated reporting probability
  q[1] = phi;
  for (d in 1:D)
    q[d+1] = exp(log(phi) * exp(-b * d));
}

model {
  // Priors
  lambda ~ lognormal(0, 2.5);
  phi ~ beta(alpha_phi, beta_phi);
  b ~ lognormal(mean_log_b, sd_log_b);

  // Likelihood
  for (d in 0:D)
    for (t in 1:(T-d))
      Y[t, d+1] ~ poisson(lambda[t] * q[d+1]);
}

generated quantities {
  vector<lower=0>[T] N;                 // number of cases
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}

