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
