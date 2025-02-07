data {
  int<lower=0> N_obs;         // Number of time points (T)
  int<lower=0> D;             // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;   // Reported cases (T x D matrix)
  // bsplines
  int<lower=1> J;                          // number of basis functions
  matrix[N_obs, J] X_spline;                   // basis function values(N_obs, J)
}

parameters {
  real<lower=0> alpha_lambda;  // Gamma prior shape parameter for lambda
  real<lower=0> beta_lambda;   // Gamma prior rate parameter for lambda
  vector<lower=0>[N_obs] lambda_t;   // Poisson intensities (λ[t]) at each time point
  //
  vector<lower=0>[J] beta;               // global spline coefficients (no clusters)
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  
  beta ~ uniform(0, 100);  // prior for spline coef
  
  // Gamma prior on Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);


  // Likelihood: Marginalized Poisson likelihood for N_t and Binomial for Y
  for (t in 1:N_obs){
    real b_t = dot_product(X_spline[t], beta); 
    for (d in 1:D){
      real q_d = 1 - exp(- b_t * d);  // q(d) as a function of delay d
      Y[t,d] ~ poisson(lambda_t[t] * q_d);
    }
  }
}

generated quantities {
  vector[N_obs] b_t_est;
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs) {
    b_t_est[t] = dot_product(X_spline[t], beta);  // posterior of  b(t) 
    N_t[t] = poisson_rng(lambda_t[t]); // Sample N_t from Poisson distribution
  }
}


