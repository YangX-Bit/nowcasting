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
      start_window <- as.Date(score) - window_day_length
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
                                     D = 20,
                                     methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
                                     compiled_models,
                                     iter_sampling = 2000, iter_warmup = 1000, refresh = 500,
                                     num_chains = 3, thin = 2,suppress_output = TRUE,
                                    posterior_draws_path = file.path(path_proj, "source", "models",
                                                                     "posterior_draws"), hypers = NULL
                                     ){
  if(is.null(case_true)){
    stop("You must input true cases.")
  }
  
  if (is.null(compiled_models) || !all(methods %in% names(compiled_models))) {
    stop("You must provide compiled models matching 'methods'.")
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
    #X_spline <- create_basis(N_obs_local, n_knots = 5) # functions to create basis
    if(nrow(data_trunc) <= D + 1){
      warning("The number of rows of the input data is smaller than number of max delay D, which might cause inaccuracy." )
    }
    
    stan_data_trunc <- c(list(T = N_obs_local, D = D, Y = data_trunc), hypers)

    # return(stan_data_trunc)
    # Fit models based on what is selected 
    for (model_name in methods) {
      compiled_model <- compiled_models[[model_name]]
      if (is.null(compiled_model)) {
        stop(paste("Model path for", model_name, "is not specified in model_paths."))
      }

      # Fit the Stan model
      sampling_code <- function() {
        compiled_model$sample(
          data = stan_data_trunc,
          iter_sampling = iter_sampling,
          iter_warmup = iter_warmup,
          chains = num_chains,
          refresh = refresh,
          thin = thin,
          output_dir = posterior_draws_path
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

### functions to get coordinates of non-NAs
find_non_na_coords <- function(mat) {
  # dimension
  N <- nrow(mat)
  D <- ncol(mat)
  
  # indices for non NAs
  non_na_indices <- which(!is.na(mat), arr.ind = TRUE)
  
  coords_df <- as.matrix(non_na_indices)
  
  # rename cols
  colnames(coords_df) <- c("row", "col")
  
  # make col consist to d = 0,...,D
  coords_df[,2] = coords_df[,2] - 1
  
  return(coords_df)
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



