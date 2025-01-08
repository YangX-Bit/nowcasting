data {
  int<lower=0> N_obs;         // Number of time points (T)
  int<lower=0> D;             // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;   // Reported cases (T x D matrix)
  int<lower=1> K;                   // number of obs in total
  array[K, 2] int<lower=1> obs_index; // coordinates
}

parameters {
  real<lower=0> alpha_lambda;  // Gamma prior shape parameter for lambda
  real<lower=0> beta_lambda;   // Gamma prior rate parameter for lambda
  vector<lower=0>[N_obs] lambda_t;   // Poisson intensities (λ[t]) at each time point
  real b0;// Parameter for delay function q(d)
  real b1;
  real<lower=0> phi;
  //vector<lower=0>[N_obs] phi;
  //real b2;
}

transformed parameters {
  vector[N_obs] b_vec;   // b(t) = exp(...)
  vector[N_obs] q_d; 
  
  for (t in 1:N_obs) {
    //real x = b0 + b1 * n + b2 * n^2;   
    real x = b0 + b1 * t;
    b_vec[t] = exp(x);
    q_d[t] =phi + exp(- b_vec[t] * (D)) ;
  }
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  b0 ~ uniform(-2, 2);  // Prior for parameter b
  b1 ~ normal(0, 0.1); 
  //b2 ~ normal(0, 5);
  phi ~ uniform(0,0.5);
  
  // Gamma prior on Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // Likelihood: Marginalized Poisson likelihood for N_t and Binomial for Y
  for(k in 1:K){
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    Y[t,d] ~ poisson(lambda_t[t] * q_d[t]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]); // Sample N_t from Poisson distribution
  }
}


