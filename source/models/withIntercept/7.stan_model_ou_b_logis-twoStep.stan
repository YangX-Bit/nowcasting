data {
  int<lower=0> N_obs;                  // Number of time points
  int<lower=0> D;                      // Maximum delay in the data (用于构建 Y[t,d] 的列数)
  array[N_obs, D] int<lower=0> Y;      // Reported cases matrix (T x D)
  int<lower=1> K;                      // Total number of observations
  array[K, 2] int<lower=1> obs_index;  // Indexes for each observation (t, d)
}

parameters {
  // Gamma 先验用到的两个参数
  real<lower=0> alpha_lambda;
  real<lower=0> beta_lambda;

  // 各时间点的 Poisson 强度 λ[t]
  vector<lower=0>[N_obs] lambda_t;

  // Ornstein-Uhlenbeck 对于 phi 的参数
  real<lower=0> alpha_phi;
  real<lower=0> mu_phi;
  real<lower=0> sigma_phi;

  // 每个时间点的 phi
  // phi[n] \in (0, 1) 用于控制 "0延迟" 报告率
  vector<lower=0, upper=1>[N_obs] phi;
}

transformed parameters {
  // 在这里用 phi 来计算 b_t
  // 采用 (1 - phi[n]) * exp(-b_t[n]*max_delay) = epsilon 的思路
  // 若不想引入太多超参, 可以把 epsilon 直接写死成一个很小的值, 如 1e-6
  real max_delay = 100;
  real epsilon = 1e-6;
  vector[N_obs] b_t;           
  matrix[N_obs, D] q_d_matrix; 

  // 计算 b_t[n]
  for (n in 1:N_obs) {
    // 注意: 当 phi[n] 很接近1时, (1 - phi[n]) 非常小, b_t[n] 可能相对较大
    //       需保证 1 - phi[n] > 0, 所以 phi[n] 不能取到1的上界
    b_t[n] = - (1 / max_delay) * log( epsilon / (1 - phi[n]) );

    // 预先算好 q_d_matrix
    for (d in 1:D){
      q_d_matrix[n, d] = 1 - (1 - phi[n]) * exp(- b_t[n] * d);
    }
  }
}

model {
  //------------------------------------------------
  // priors
  //------------------------------------------------
  // 对 alpha_lambda, beta_lambda 给予一些较宽松的先验
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);

  // phi 的 OU 过程参数先验，可再酌情调整
  alpha_phi ~ normal(0.5, 0.2);
  mu_phi    ~ uniform(0, 1);
  sigma_phi ~ normal(0.1, 0.05);

  // Gamma prior for Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);

  // Ornstein-Uhlenbeck process for phi
  phi[1] ~ normal(mu_phi, sigma_phi);
  for (t in 2:N_obs){
    phi[t] ~ normal(phi[t-1] + alpha_phi * (mu_phi - phi[t-1]), sigma_phi);
  }

  //------------------------------------------------
  // likelihood
  //------------------------------------------------
  // 各观测点 (t, d) 的泊松似然： Y[t,d] ~ Poisson( lambda[t] * q_d_matrix[t,d] )
  for (k in 1:K){
    int t = obs_index[k, 1];
    int d = obs_index[k, 2];
    Y[t, d] ~ poisson(lambda_t[t] * q_d_matrix[t, d]);
  }
}

generated quantities {
  // 这里可根据自己需要, 生成一些感兴趣的量
  // 比如每个 t 的 "最终规模" (不区分延迟) 的一个预测
  vector[N_obs] N_t;
  for (t in 1:N_obs){
    N_t[t] = poisson_rng(lambda_t[t]);
  }
}
