test <- out_2$b_spline[[4]]
test_coef <- test$summary()


test_coef[test_coef$variable == "N_t[1]", ]


Nts <- test$draws(variables = "N_t")
Nts

library(posterior) 
posterior_samples <- test$draws()
posterior_df <- as_draws_df(posterior_samples)
head(posterior_df)

install.packages("bayesplot")
library(bayesplot)

mcmc_hist(posterior_samples, pars = c("N_t[30]"))


test_1 <- nowcasting_moving_window(data_2$case_reported, scoreRange = as.Date("2024-02-26"), 
                         case_true = data_2$case_true,
                         start_date = as.Date("2024-01-01"),
                         D = D_trunc, sigma_b = 0.1, seeds = 666,
                         #models_to_run = c("fixed_q", "fixed_b", "b_poly","b_spline"),
                         models_to_run = c("b_spline"),
                         compiled_models = compiled_models,
                         iter_sampling = 2000, iter_warmup = 1000, refresh = 0,
                         num_chains = 3, suppress_output = T)

post <- test_1$b_spline[[1]]$draws()
post_samples <- as_draws_df(post)

mcmc_hist(post_samples, pars = c("N_t[30]"))

N_array <- Z_array <- b_array <-  numeric(57)
for (i in 1:57) {
  N_array[i] <- round(mean(as.numeric(post_samples[, c(paste0("N_t[", i, "]"))][[1]])),2)
  Z_array[i] <- round(mean(as.numeric(post_samples[, c(paste0("Z_t[", i, "]"))][[1]])),2)
  b_array[i] <- round(mean(as.numeric(post_samples[, c(paste0("b_t[", i, "]"))][[1]])),2)
}
qd <-  1 - exp(- b_array * D)

as.numeric(test_1$case_true[[1]]) - Z_array - N_array

as.numeric(test_1$case_true[[1]])[27]
Z_array[27]
N_array[27]

plot(data_2$b - b_array)
data_2$case_reported[c(25:30),]


###########


model1 <-out_ou_b_FR$ou_b[[1]]



# 提取 b_t 和 phi 参数的采样数据
model1_bt <- model1$draws(variables = c("b_t"), format = "draws_df")
model1_phi <- model1$draws(variables = c("phi"), format = "draws_df")

# 提取前 30 列（假设是感兴趣的参数）
bt_samples <- model1_bt[, c(1:30)]
phi_samples <- model1_phi[, c(1:30)]

# 计算均值
mean_values_bt <- apply(bt_samples, 2, mean)
mean_values_phi <- apply(phi_samples, 2, mean)

# 计算上下限 (95% CI)
ci_lower_bt <- apply(bt_samples, 2, function(x) quantile(x, 0.025)) # 2.5% 分位数
ci_upper_bt <- apply(bt_samples, 2, function(x) quantile(x, 0.975)) # 97.5% 分位数

ci_lower_phi <- apply(phi_samples, 2, function(x) quantile(x, 0.025))
ci_upper_phi <- apply(phi_samples, 2, function(x) quantile(x, 0.975))

# 打印结果
print(mean_values_bt)
print(ci_lower_bt)
print(ci_upper_bt)

print(mean_values_phi)
print(ci_lower_phi)
print(ci_upper_phi)

library(ggplot2)

library(ggplot2)
library(dplyr)

# 假设 exp_grow_func, model1_bt, model1_phi 已定义
# 数据准备
x <- 1:15  # x 轴范围

# 提取 b_t 和 phi 参数的均值和区间
mean_values_bt <- apply(model1_bt[, c(1:30)], 2, mean)
ci_lower_bt <- apply(model1_bt[, c(1:30)], 2, function(x) quantile(x, 0.025))
ci_upper_bt <- apply(model1_bt[, c(1:30)], 2, function(x) quantile(x, 0.975))

mean_values_phi <- apply(model1_phi[, c(1:30)], 2, mean)
ci_lower_phi <- apply(model1_phi[, c(1:30)], 2, function(x) quantile(x, 0.025))
ci_upper_phi <- apply(model1_phi[, c(1:30)], 2, function(x) quantile(x, 0.975))

library(ggplot2)
library(dplyr)

# 定义绘图函数
plot_multiple_curves <- function(model1, x_range, num_rows, num_cols, exp_grow_func, ou_b_FR) {
  # 提取 b_t 和 phi 参数的采样数据
  model1_bt <- model1$draws(variables = c("b_t"), format = "draws_df")
  model1_phi <- model1$draws(variables = c("phi_t"), format = "draws_df")
  
  # 提取前 30 列（假设是感兴趣的参数）
  bt_samples <- model1_bt[, c(1:30)]
  phi_samples <- model1_phi[, c(1:30)]
  
  # 计算均值和上下限 (95% CI)
  mean_values_bt <- apply(bt_samples, 2, mean)
  ci_lower_bt <- apply(bt_samples, 2, function(x) quantile(x, 0.025))
  ci_upper_bt <- apply(bt_samples, 2, function(x) quantile(x, 0.975))
  
  mean_values_phi <- apply(phi_samples, 2, mean)
  ci_lower_phi <- apply(phi_samples, 2, function(x) quantile(x, 0.025))
  ci_upper_phi <- apply(phi_samples, 2, function(x) quantile(x, 0.975))
  
  # 数据生成
  plot_data <- lapply(1:(num_rows * num_cols), function(i) {
    # 黑色线数据
    y_black <- exp_grow_func(x_range, ou_b_FR$b_t[c(1:30)][i], ou_b_FR$phi[c(1:30)][i])
    
    # 红色线数据及区间
    y_red <- exp_grow_func(x_range, mean_values_bt[i], mean_values_phi[i])
    y_red_lower <- exp_grow_func(x_range, ci_lower_bt[i], ci_lower_phi[i])
    y_red_upper <- exp_grow_func(x_range, ci_upper_bt[i], ci_upper_phi[i])
    
    # 数据框
    data.frame(
      x = rep(x_range, 2),  # 重复 x
      y = c(y_black, y_red),
      group = c(rep("Black", length(x_range)), rep("Red", length(x_range))),
      lower = c(rep(NA, length(x_range)), y_red_lower),  # 下限仅对红色有效
      upper = c(rep(NA, length(x_range)), y_red_upper),  # 上限仅对红色有效
      plot_id = i  # 用于区分不同子图
    )
  }) %>% bind_rows()  # 合并所有子图数据
  
  # 绘图
  ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_line(size = 1) +  # 绘制线
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Black" = "black", "Red" = "red")) +
    scale_fill_manual(values = c("Black" = NA, "Red" = "red")) +
    facet_wrap(~plot_id, nrow = num_rows, ncol = num_cols) +  # 网格布局
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Comparison of Curves with Uncertainty Intervals",
         x = "X", y = "Y")
}

# 调用函数示例
plot_multiple_curves(
  model1 = model1,                      # cmdstanr 的模型对象
  x_range = 1:15,                       # X 轴范围
  num_rows = 5,                         # 图的行数
  num_cols = 5,                         # 图的列数
  exp_grow_func = exp_grow_func,        # 用户定义的增长函数
  ou_b_FR = ou_b_FR                     # 包含 b_t 和 phi 数据的对象
)

