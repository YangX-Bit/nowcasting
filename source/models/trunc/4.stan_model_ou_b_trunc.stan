data {
  int<lower=0> N_obs;                    // 时间点数 T
  int<lower=0> D;                        // 最大延迟 D
  array[N_obs, D] int<lower=0> Y;        // 报告病例数矩阵 (T x D)
  int<lower=1> K;                        // 总观测数 K
  array[K, 2] int<lower=1> obs_index;     // 每个观测的 (t, d) 索引
}

parameters {
  real<lower=0> alpha_lambda;            // Gamma分布 shape参数
  real<lower=0> beta_lambda;             // Gamma分布 rate参数
  vector<lower=0>[N_obs] lambda_t;       // 每个时间点的Poisson强度 λ[t]

  // Ornstein-Uhlenbeck过程参数
  real<lower=0> alpha;                    // 均值回归速率
  real mu;                                 // 长期均值
  real<lower=0> sigma_ou;                 // 扩散系数

  vector[N_obs] log_b_t;                  // 对数b_t参数，用于稳定性
}

transformed parameters {
  vector[N_obs] b_t;
  matrix[N_obs, D] q_d_matrix;

  // 计算b_t[t] = exp(log_b_t[t])，确保b_t[t]为正
  for (n in 1:N_obs){
    b_t[n] = exp(log_b_t[n]);
  }

  // 预计算所有q_d
  for (n in 1:N_obs){
    for (d in 1:D){
      q_d_matrix[n, d] = 1 - exp(- b_t[n] * d);
    }
  }
}

model {
  // 先验分布
  alpha_lambda ~ uniform(0, 10);
  beta_lambda  ~ uniform(0, 10);

  alpha ~ normal(0.5, 0.2);       // 根据业务调整
  mu    ~ normal(0.5, 0.2);       // 根据业务调整
  sigma_ou ~ normal(0.1, 0.05);   // 根据业务调整

  // Poisson强度的Gamma先验
  lambda_t ~ gamma(alpha_lambda, beta_lambda);

  // OU过程：对数空间的b_t
  log_b_t[1] ~ normal(log(mu), sigma_ou);
  for (t in 2:N_obs){
    log_b_t[t] ~ normal(log_b_t[t-1] + alpha * (log(mu) - log_b_t[t-1]), sigma_ou);
  }

  // 似然函数：向量化
  for (k in 1:K){
    int t = obs_index[k,1];
    int d = obs_index[k,2];
    Y[t,d] ~ poisson(lambda_t[t] * q_d_matrix[t, d]);
  }
}

generated quantities {
  vector<lower=0>[N_obs] N_t;
  for (t in 1:N_obs){
    N_t[t] = poisson_rng(lambda_t[t]); // 生成最终病例数
  }
}
