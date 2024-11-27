nowcasting_moving_window <- function(data, scoreRange,
                                     N_obs = 60, D = 15, sigma_b = 0.1, seeds = 123,
                                     path_p_change, path_p_fixed,
                                     iter = 2000, warmup = 1000, refresh = 500){
  # plot list
    plot_list <- list()
  # fit list
    model_p_fixed_list <- list()
    model_p_change_list <- list()
  i <- 1
  for (i in 1:length(scoreRange)) {
    #What's "today"
    now <- scoreRange[i]
    # show the status
    cat(paste("====================\nnow=",now,
              " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
    
    #when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
    
    # Nowcast #
    # cut the length of data
    data_use <- head(data, N_obs - length(scoreRange) + i)
    
    # truncated version of data
    data_trunc <- create_triangular_data(data_use, if_zero = F)
    N_obs_local <- nrow(data_trunc)
    
    # coordinates for non-NAs
    indices_data_trunc <- find_non_na_coords(data_trunc)
    
    data_trunc[is.na(data_trunc)] <- 0 # to avoid NAs in data
    
    X_spline <- create_basis(N_obs_local, n_knots = 5) # functions to create basis
    
    # input list
    stan_data_trunc <- list(N_obs = N_obs_local, D = D + 1, Y = data_trunc,
                            K = nrow(indices_data_trunc), obs_index = indices_data_trunc,
                            J = ncol(X_spline),
                            X_spline = X_spline,
                            sigma_b = sigma_b)
    
    fit_trunc <- stan(
      file = path_p_change,  
      data = stan_data_trunc, 
      iter = iter, warmup = warmup, chains = 3, seed = seeds,
      #control = list(adapt_delta = 0.96, max_treedepth = 15),
      refresh = refresh
    )

    
    fit_trunc_fixped_q <- stan(
      file = path_p_fixed,  
      data = stan_data_trunc, 
      iter = iter, warmup = warmup, chains = 3, seed = seeds,
      refresh = refresh, 
    )
    
    # extract parameters
    samples_nt <- rstan::extract(fit_trunc, pars = "N_t")$N_t
    samples_nt_fixped_q <- rstan::extract(fit_trunc_fixped_q, pars = "N_t")$N_t
    
    p <- nowcasts_plot(samples_nt, samples_nt_fixped_q, N_obs = N_obs_local)
    
    plot_list[[i]] <- p
    model_p_fixed_list[[i]] <- fit_trunc
    model_p_change_list[[i]] <- fit_trunc_fixped_q
  }
  # output 
  return(list(
    plots = plot_list[!sapply(plot_list, is.null)],
    models = list(
      fixed = model_p_fixed_list[!sapply(model_p_fixed_list, is.null)],
      change = model_p_change_list[!sapply(model_p_change_list, is.null)]
    )
  ))
}

