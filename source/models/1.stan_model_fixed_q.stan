data {
  // Dimensions
  int<lower=0> T;                   // number of time points
  int<lower=0> D;                   // maximum delay
  int<lower=1> nobs;                // number of obs in total
  // Data
  array[T, D+1] int<lower=0> Y;     // reported cases (t x d matrix)
  array[nobs, 2] int<lower=0> obs_index; // coordinates
  vector<lower=0>[T] alpha;         // gamma prior shape parameter for lambda
}

parameters {
  vector<lower=0>[T] lambda;        // expected number of cases
  simplex[D+1] p;                   // reporting probability
}

model {
  // Priors
  // p ~ uniform over a simplex (flat Dirichlet distribution)
  lambda[t] ~ lognormal(0, 3);

  // Gamma prior on Poisson intensities (lambda)
  for(t in 1:T){
    lambda[t] ~ lognormal(0, 3);
  }

  // Likelihood: Marginalized Poisson likelihood for N and Binomial for Y
  for(i in 1:nobs){
    int t = obs_index[i,1];
    int d = obs_index[i,2];

    real q_d = sum(p[1:(d+1)]);
    Y[t, (d+1)] ~ poisson(lambda[t] * q_d);
  }
}

generated quantities {
  vector<lower=0>[T] N;
  for (t in 1:T) {
    N[t] = poisson_rng(lambda[t]); // Sample N from Poisson distribution
  }
}

