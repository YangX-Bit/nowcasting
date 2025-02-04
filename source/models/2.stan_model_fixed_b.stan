data {
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
}

parameters {
  vector<lower=0>[T] lambda;            // expected number of cases
  real<lower=0> b;                      // rate of accumulated reporting probability
  real<lower=0, upper=1> phi;           // delayed reporting probability
}

transformed parameters {
  vector<lower=0, upper=1>[D+1] q;      // accumulated reporting probability
  q[1] = 1 - phi;
  for (d in 1:D)
    q[d+1] = 1 - phi * exp(- b * d);
}

model {
  // Priors
  b ~ lognormal(0, 0.3);
  phi ~ uniform(0, 1);
  lambda ~ lognormal(0, 3);

  // Likelihood
  for (t in 1:T) {
    int Dmax = min(T-t, D);
    for (d in 0:Dmax)
      Y[t, d+1] ~ poisson(lambda[t] * q[d+1]);
  }
}

generated quantities {
  vector<lower=0>[T] N;                 // number of cases
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}

