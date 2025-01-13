data {
  int<lower=0> N_obs;          
  int<lower=1> D;              
  array[N_obs, D] int<lower=0> Y;   
  int<lower=1> K;              
  array[K, 2] int<lower=1> obs_index;  
}

transformed data {
  // 预计算延迟时间序列
  vector[D-1] delays;
  for (d in 1:(D-1)) {
    delays[d] = d - 1.0;  // 从0开始的延迟时间
  }
}

parameters {
  real<lower=0> alpha_lambda; 
  real<lower=0> beta_lambda;  
  vector<lower=0>[N_obs] lambda_t;
  real<lower=0> alpha_b; 
  real<lower=0.05, upper=1> mu_b; 
  real<lower=0> sigma_b_star;
  real<lower=0> alpha_phi; 
  real<lower=0, upper=1> mu_phi;
  real<lower=0> sigma_phi_star;
  vector[N_obs] b_star;
  vector[N_obs] phi_star;
}

transformed parameters {
  vector<lower=0.05,upper=1>[N_obs] b_t = 0.05 + 0.95 * inv_logit(b_star);
  vector<lower=0,upper=1>[N_obs] phi_t = inv_logit(phi_star);
  real mu_b_star = logit((mu_b - 0.05) / 0.95);
  real mu_phi_star = logit(mu_phi);
  matrix[N_obs, D] p_d_matrix;
  
  // 向量化计算p_d_matrix
  p_d_matrix[,1] = phi_t;  // day 0
  for(n in 1:N_obs) {
    p_d_matrix[n,2:D] = (1 - phi_t[n]) * 
      (exp(-b_t[n] * delays[1:(D-1)]) - exp(-b_t[n] * (delays[1:(D-1)] + 1)))';
  }
}

model {
  // 先验
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // OU过程先验 
  alpha_b ~ normal(0.5, 0.2);
  mu_b ~ normal(0.5, 0.2);
  sigma_b_star ~ normal(0.1, 0.05);
  alpha_phi ~ normal(0.5, 0.2);
  mu_phi ~ normal(0.2, 0.05) T[0,1];
  sigma_phi_star ~ normal(0.1, 0.05);
  
  // OU初始化
  b_star[1] ~ normal(mu_b_star, sigma_b_star);
  phi_star[1] ~ normal(mu_phi_star, sigma_phi_star);
  
  // OU演化 - 向量化
  b_star[2:N_obs] ~ normal(b_star[1:(N_obs-1)] + 
    alpha_b * (mu_b_star - b_star[1:(N_obs-1)]), sigma_b_star);
  phi_star[2:N_obs] ~ normal(phi_star[1:(N_obs-1)] + 
    alpha_phi * (mu_phi_star - phi_star[1:(N_obs-1)]), sigma_phi_star);
    
  // 似然
  for(k in 1:K) {
    Y[obs_index[k,1], obs_index[k,2]] ~ 
      poisson(lambda_t[obs_index[k,1]] * p_d_matrix[obs_index[k,1], obs_index[k,2]]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for(t in 1:N_obs) {
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
