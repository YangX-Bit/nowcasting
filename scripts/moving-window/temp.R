#nowcast 

data$date
#
 

# length for nowcast
N_obs <- 60
D <- 15
k_back <- 10  

now <- as.Date("2024-02-01")
when <- seq(now-k_back+1, now, by="1 day")

case_data <- data$case_reported





########### loop to fit ##########
scoreRange <- seq(as.Date("2024-01-20"),as.Date("2024-01-21"),by="1 day")

scoreRange <- seq(as.Date("2024-01-20"),as.Date("2024-02-28"),by="1 day")



plot_list <- list()

for (i in 40:length(scoreRange)) {
  #What's "today"
  now <- scoreRange[i]
  # show the status
  cat(paste("====================\nnow=",now,
            " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
  
  #when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
  
  # Nowcast #
  # get data
  data_use <- head(data$case_reported, N_obs - length(scoreRange) + i)
  
  # truncated version of data
  data_trunc <- create_triangular_data(data_use, if_zero = T)
  N_obs_local <- nrow(data_trunc)
  
  # coordinates for non-NAs
  indices_data_trunc <- find_non_na(data_trunc)
  
  X_spline <- create_basis(N_obs_local, n_knots = 5)

  # input list
  stan_data_trunc <- list(N_obs = N_obs_local, D = D + 1, Y = data_trunc,
                          K = nrow(indices_data_trunc), obs_index = indices_data_trunc,
                          J = ncol(X_spline),
                          X_spline = X_spline,
                          sigma_b = 0.1)
  
  fit_trunc <- stan(
    file = file.path(path_proj, "source", "models",
                     "trunc", "stan_model_dep-spline-trunc.stan"),  
    data = stan_data_trunc, 
    iter = 2000, chains = 3, seed = 123
  )
  
  # extract parameters
  samples_nt <- rstan::extract(fit_trunc, pars = "N_t")$N_t
  
  nowcasts <- data.frame(mean = apply(samples_nt, 2, mean),
                         lower = apply(samples_nt, 2, quantile, probs = 0.025),
                         upper = apply(samples_nt, 2, quantile, probs = 0.975),
                         date =  data$date[1:N_obs_local],
                         case_true = data$case_true[1:N_obs_local],
                         case_reported = apply(data_trunc, 1, max))
  
  p <- ggplot(nowcasts, aes(x = date)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.5) +
    geom_line(aes(y = mean), color = "blue") +
    
    geom_col(aes(y = case_true), fill = "red", alpha = 0.6) +
    
    geom_col(aes(y = case_reported), fill = "blue", alpha = 0.4) +
    
    labs(title = "Nowcast with True and Reported Cases",
         x = "Date",
         y = "Cases / Nowcast") +
    theme_minimal()
  
  plot_list[[i]] <- p
  
  print(p)
}

plot_list[[1]]


N_obs_local <- 60
find_non_na_coords(matrix(c(1,2,NA,4,5,6), nrow = 2))
