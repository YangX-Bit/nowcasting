data {
  // Data
  int<lower=0> T;                       // number of time points
  int<lower=0> D;                       // maximum delay
  array[T, D+1] int<lower=0> Y;         // reported cases (t x (d+1) matrix)
  // Hyperparameters for phi and b (as before)
  real mean_logit_phi;
  real<lower=0> sd_logit_phi;
  real mean_log_b;
  real<lower=0> sd_log_b;
  // Hyperparameters for alpha
  real mean_logit_alpha;
  real<lower=0> sd_logit_alpha;
}

parameters {
  vector<lower=0>[T] lambda;            // expected number of cases

  // b and phi as OU processes
  vector[T] log_b;                      // log rate of accumulated reporting probability
  vector[T] logit_phi;                  // instantaneous reporting prob (logit)

  // alpha as OU process
  vector[T] logit_alpha;                // long-term reporting completeness (logit)

  // OU hyperparameters
  real mu_log_b;                        
  real mu_logit_phi;                    
  real mu_logit_alpha;                  

  real<lower=0> theta_log_b;            
  real<lower=0> theta_logit_phi; 
  real<lower=0> theta_logit_alpha;

  real<lower=0> sigma_log_b;            
  real<lower=0> sigma_logit_phi;        
  real<lower=0> sigma_logit_alpha;      
}

transformed parameters {
  vector<lower=0>[T] b       = exp(log_b);
  vector<lower=0,upper=1>[T] phi     = inv_logit(logit_phi);
  vector<lower=0,upper=1>[T] alpha   = inv_logit(logit_alpha);
  matrix[T, D+1] q;                    // accumulated reporting probability with alpha

  for (d in 0:D)
    for (t in 1:(T-d))
      q[t, d+1] = alpha[t] - (alpha[t] - phi[t]) * exp(-b[t] * d);

}

model {
  // Priors
  lambda ~ lognormal(0, 2.5);

  mu_log_b        ~ normal(mean_log_b, sd_log_b);
  mu_logit_phi    ~ normal(mean_logit_phi, sd_logit_phi);
  mu_logit_alpha  ~ normal(mean_logit_alpha, sd_logit_alpha);

  theta_log_b       ~ lognormal(0, 1);
  theta_logit_phi   ~ lognormal(0, 1);
  theta_logit_alpha ~ lognormal(0, 1);

  sigma_log_b       ~ lognormal(-2, 1);
  sigma_logit_phi   ~ lognormal(-2, 1);
  sigma_logit_alpha ~ lognormal(-2, 1);

  // OU processes
  log_b[1]        ~ normal(mu_log_b, sigma_log_b);
  logit_phi[1]    ~ normal(mu_logit_phi, sigma_logit_phi);
  logit_alpha[1]  ~ normal(mu_logit_alpha, sigma_logit_alpha);

  for (t in 2:T) {
    log_b[t] ~ normal(
      log_b[t-1] + theta_log_b * (mu_log_b - log_b[t-1]),
      sigma_log_b
    );
    logit_phi[t] ~ normal(
      logit_phi[t-1] + theta_logit_phi * (mu_logit_phi - logit_phi[t-1]),
      sigma_logit_phi
    );
    logit_alpha[t] ~ normal(
      logit_alpha[t-1] + theta_logit_alpha * (mu_logit_alpha - logit_alpha[t-1]),
      sigma_logit_alpha
    );
  }

  // Likelihood
  for (d in 0:D)
    for (t in 1:(T-d))
      Y[t, d+1] ~ poisson(lambda[t] * q[t, d+1]);
}

generated quantities {
  vector<lower=0>[T] N;
  for (t in 1:T)
    N[t] = poisson_rng(lambda[t]);
}
