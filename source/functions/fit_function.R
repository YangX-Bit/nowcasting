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


nowcasting_moving_window <- function(data, scoreRange, case_true = NULL,
                                     start_date = NULL, predict_length = NULL,
                                     D = 20, sigma_b = 0.1, seeds = 123,
                                     models_to_run = c("fixed_q", "fixed_b", "b_poly", "b_spline"),
                                     compiled_models,
                                     iter_sampling = 2000, iter_warmup = 1000, refresh = 500,
                                     num_chains = 3, suppress_output = TRUE){
  if(is.null(case_true)){
    stop("You must input true cases.")
  }
  
  if (is.null(compiled_models) || !all(models_to_run %in% names(compiled_models))) {
    stop("You must provide compiled models matching 'models_to_run'.")
  }
  
  # get the date
  if(is.null(start_date)){ start_date = rownames(data)[1] 
  }else {
    data <- data[rownames(data) >= start_date,]
  } 
  
  # prepare data
  data <- as.matrix(data)
  scoreRange <- as.Date(scoreRange)
  data_list <- slice_data(data, scoreRange, 
                          start_date = start_date, window_day_length = predict_length)
  scoreRange <- tail(scoreRange, length(data_list)) #remove invalid scoring date
  # result list
  model_fits <- list()
  for (i in 1:length(scoreRange)) {
    #What's "today"
    now <- scoreRange[i]
    # show the status
    cat(paste("====================\nnow=",now,
              " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
    
    # prepare the data for Stan
    data_use <- data_list[[i]]
    data_trunc <- create_triangular_data(data_use, if_zero = F)
   
    # information for plot
    model_fits[["case_true"]][[i]] <- case_true[rownames(case_true) 
                                                %in% rownames(data_use), , drop = FALSE]
    model_fits[["case_reported"]][[i]] <- extract_last_valid(data_trunc)
    model_fits[["dates"]][[i]] <- as.Date(rownames(data_use))
    
    N_obs_local <- nrow(data_trunc) # num of obs
    indices_data_trunc <- find_non_na_coords(data_trunc) # coordinates for non-NAs
    data_trunc[is.na(data_trunc)] <- 0 # to avoid NAs in data
    X_spline <- create_basis(N_obs_local, n_knots = 5) # functions to create basis
    if(nrow(data_trunc) <= D + 1){
      warning("The number of rows of the input data is smaller than number of max delay D, which might cause inaccuracy." )
    }
    
    # input list
    stan_data_trunc <- list(N_obs = N_obs_local, D = D + 1, Y = data_trunc,
                            K = nrow(indices_data_trunc), obs_index = indices_data_trunc,
                            J = ncol(X_spline), X_spline = X_spline, sigma_b = sigma_b)
    # return(stan_data_trunc)
    # Fit models based on what is selected 
    for (model_name in models_to_run) {
      compiled_model <- compiled_models[[model_name]]
      if (is.null(compiled_model)) {
        stop(paste("Model path for", model_name, "is not specified in model_paths."))
      }

      # Fit the Stan model
      sampling_code <- function() {
        compiled_model$sample(
          data = stan_data_trunc,
          seed = seeds,
          iter_sampling = iter_sampling,
          iter_warmup = iter_warmup,
          chains = num_chains,
          refresh = refresh
        )
      }
      
      if (suppress_output) {
        fit <- suppressWarnings(suppressMessages(sampling_code()))
      } else {
        fit <- sampling_code()
      }

      # Store the result
      model_fits[[model_name]][[i]] <- fit
    }
  }
  return(model_fits)
}

normalize_matrix_columns <- function(matrix_data) {
  # result matrix
  result_matrix <- matrix_data
  
  # num of col
  ncols <- ncol(matrix_data)
  
  for (i in 1:nrow(matrix_data)) {
    # avoid this if last col is zero
    if (matrix_data[i, ncols] != 0) {
      result_matrix[i, ] <- matrix_data[i, ] / matrix_data[i, ncols]
    }
  }
  
  return(result_matrix)
}


# measurement metrics
calculate_metrics <- function(df, methods = c("fixed_q", "fixed_b", "b_poly", "b_spline")) {
  # Ensure required columns exist
  required_columns <- c("date", "case_true")
  for (method in methods) {
    required_columns <- c(required_columns, 
                          paste0("mean_", method), 
                          paste0("lower_", method), 
                          paste0("upper_", method))
  }
  if (!all(required_columns %in% names(df))) {
    stop("The data frame is missing required columns for the specified methods.")
  }
  
  # Initialize result list
  results <- list()
  
  # Calculate metrics for each method
  for (method in methods) {
    mean_col <- paste0("mean_", method)
    lower_col <- paste0("lower_", method)
    upper_col <- paste0("upper_", method)
    
    # Calculate errors
    error <- df[[mean_col]] - df$case_true
    abs_error <- abs(error)
    rel_error <- error / df$case_true
    abs_rel_error <- abs(rel_error)
    
    # Metrics
    RMSE <- sqrt(mean(error^2, na.rm = TRUE))
    RMSPE <- sqrt(mean(rel_error^2, na.rm = TRUE)) * 100
    MAE <- mean(abs_error, na.rm = TRUE)
    MAPE <- mean(abs_rel_error, na.rm = TRUE) * 100
    
    # Prediction interval metrics
    interval_width <- mean(df[[upper_col]] - df[[lower_col]], na.rm = TRUE)
    coverage <- mean((df$case_true >= df[[lower_col]]) & (df$case_true <= df[[upper_col]]), na.rm = TRUE)
    
    # Store results (force to numeric to avoid list nesting)
    results[[method]] <- c(
      RMSE = RMSE,
      RMSPE = RMSPE,
      MAE = MAE,
      MAPE = MAPE,
      Interval_Width = interval_width,
      Coverage_Rate = coverage
    )
  }
  
  # Convert to a data frame
  results_df <- do.call(rbind, results)
  results_df <- as.data.frame(results_df)  # Ensure it's a data frame
  results_df$Method <- rownames(results_df)  # Add method names as a column
  rownames(results_df) <- NULL
  
  return(results_df)
}




