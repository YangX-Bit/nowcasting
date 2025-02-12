data {
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
}

transformed data {
  vector[D+1] p_alpha = rep_vector(1.0, D+1);
}

parameters {
  vector<lower=0>[T] lambda;            // expected number of cases
  simplex[D+1] p;                       // reporting probability
}

transformed parameters {
  vector[D+1] q = cumulative_sum(p);    // accumulated reporting probability
}

model {
  // Priors
  lambda ~ lognormal(0, 3);
  p ~ dirichlet(p_alpha);               // simplex uniform

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

