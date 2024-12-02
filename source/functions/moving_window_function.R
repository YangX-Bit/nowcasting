slice_data <- function(data, scoreRange, 
                       start_date = NULL, window_day_length = NULL) {
  # Validate input arguments
  if (is.null(start_date) && is.null(window_day_length)) {
    stop("Either start_date or window_length must be provided.")
  }
  
  dates <- as.Date(rownames(data))
  
  # Initialize the result list
  result <- list()
  
  if (!is.null(start_date)) {
    # Slice based on start_date and scoreRange
    for (score in scoreRange) {
      end_date <- score
      slice <- data[dates >= as.Date(start_date) & dates <= end_date, , drop = FALSE]
      if (nrow(slice) > 0) {
        result[[paste("Range_from", start_date, "to", as.Date(end_date), sep = "_")]] <- as.matrix(slice)
      }
    }
  } else if (!is.null(window_day_length)) {
    # Slice based on window_length and scoreRange
    for (score in scoreRange) {
      start_window <- as.Date(score) - days(window_day_length)
      slice <- data[dates >= start_window & dates <= as.Date(score), , drop = FALSE]
      if (nrow(slice) > 0) {
        result[[paste(window_day_length,"days_window_up_to", last(rownames(slice)), sep = "_")]] <- as.matrix(slice)
      }
    }
  }
  return(result)
}


a <- slice_data(data, scoreRange, start_date = as.Date("2009-07-01"), window_day_length = NULL)
length(scoreRange)
last(rownames(a[[6]]))

nowcasting_moving_window <- function(data, scoreRange, case_true = NULL,
                                     start_date = NULL, predict_length = NULL,
                                     D = 20, sigma_b = 0.1, seeds = 123,
                                     path_p_change, path_p_fixed,
                                     iter = 2000, warmup = 1000, refresh = 500,
                                     num_chains = 3){
    if(is.null(case_true)){
      stop("You must input true cases.")
    }
  # get the date
    if(is.null(start_date)){ start_date = rownames(data)[1] 
    }else {
      data <- data[rownames(data) >= start_date,]
    } 
    
    data <- as.matrix(data)
    scoreRange <- as.Date(scoreRange)
    
    # create used data
    data_list <- slice_data(data, scoreRange, 
                            start_date = start_date, window_day_length = predict_length)
    scoreRange <- tail(scoreRange, length(data_list)) #remove invalid scoring date
    
  # plot list
    plot_list <- list()
  # fit list
    model_p_fixed_list <- list()
    model_p_change_list <- list()
  for (i in 1:length(scoreRange)) {
    #What's "today"
    now <- scoreRange[i]
    # show the status
    cat(paste("====================\nnow=",now,
              " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
    #when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
    
    # Nowcast #
    # cut the length of data
    data_use <- data_list[[i]]
    
    # cut the true cases
    case_true_temp <- case_true[rownames(case_true) %in% rownames(data_use), , drop = FALSE]
    
    # truncated version of data
    data_trunc <- create_triangular_data(data_use, if_zero = F)
    # For real data the last valid column is the last non-NA column
    case_reported <- extract_last_valid(data_trunc)
    N_obs_local <- nrow(data_trunc) # num of obs
    # coordinates for non-NAs
    indices_data_trunc <- find_non_na_coords(data_trunc)
    data_trunc[is.na(data_trunc)] <- 0 # to avoid NAs in data
    
    X_spline <- create_basis(N_obs_local, n_knots = 5) # functions to create basis
    
    if(nrow(data_trunc) <= D + 1){
      warning("The number of rows of the input data is smaller than number of max delay D, which might cause inaccuracy." )
    }
    
    # input list
    stan_data_trunc <- list(N_obs = N_obs_local, D = D + 1, Y = data_trunc,
                            K = nrow(indices_data_trunc), obs_index = indices_data_trunc,
                            J = ncol(X_spline),
                            X_spline = X_spline,
                            sigma_b = sigma_b)
    
    fit_trunc <- stan(
      file = path_p_change,
      data = stan_data_trunc,
      iter = iter, warmup = warmup, chains = num_chains, seed = seeds,
      #control = list(adapt_delta = 0.96, max_treedepth = 15),
      refresh = refresh
    )

    fit_trunc_fixped_q <- stan(
      file = path_p_fixed,
      data = stan_data_trunc,
      iter = iter, warmup = warmup, chains = num_chains, seed = seeds,
      refresh = refresh,
    )

    # extract parameters
    samples_nt <- rstan::extract(fit_trunc, pars = "N_t")$N_t
    samples_nt_fixped_q <- rstan::extract(fit_trunc_fixped_q, pars = "N_t")$N_t
    
    p <- nowcasts_plot(samples_nt, samples_nt_fixped_q, N_obs = N_obs_local,
                       dates = as.Date(rownames(data_use)),
                       case_true = case_true_temp, case_reported = case_reported)
    
    plot_list[[i]] <- p
    model_p_fixed_list[[i]] <- fit_trunc
    model_p_change_list[[i]] <- fit_trunc_fixped_q
  }
  # output 
  return(list(
    plots = plot_list[!sapply(plot_list, is.null)],
    fixed = model_p_fixed_list[!sapply(model_p_fixed_list, is.null)],
    change = model_p_change_list[!sapply(model_p_change_list, is.null)]
  ))
}

