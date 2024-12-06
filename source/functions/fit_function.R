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

    # fit_trunc_fixped_q <- stan(
    #   file = path_p_fixed,
    #   data = stan_data_trunc,
    #   iter = iter, warmup = warmup, chains = num_chains, seed = seeds,
    #   refresh = refresh,
    # )

    # extract parameters
    samples_nt <- rstan::extract(fit_trunc, pars = "N_t")$N_t
    #samples_nt_fixped_q <- rstan::extract(fit_trunc_fixped_q, pars = "N_t")$N_t
    
    # p <- nowcasts_plot(samples_nt, samples_nt_fixped_q, N_obs = N_obs_local,
    #                    dates = as.Date(rownames(data_use)),
    #                    case_true = case_true_temp, case_reported = case_reported)
    # 
    # plot_list[[i]] <- p
    model_p_change_list[[i]] <- fit_trunc
    #model_p_fixed_list[[i]] <- fit_trunc_fixped_q
  }
  # output 
  return(list(
    #plots = plot_list[!sapply(plot_list, is.null)],
    #fixed = model_p_fixed_list[!sapply(model_p_fixed_list, is.null)],
    change = model_p_change_list[!sapply(model_p_change_list, is.null)]
  ))
}

matrix_data <- data_41002[, c(1:21)]
# check the q shape
fit_and_plot <- function(matrix_data) {
  if (!is.matrix(matrix_data)) stop("Input must be a matrix.")
  
  # 获取矩阵的行数
  n_rows <- nrow(matrix_data)
  
  # 创建一个绘图窗口
  par(mfrow = c(ceiling(sqrt(n_rows)), ceiling(sqrt(n_rows)))) # 设置布局
  
  # 定义拟合函数
  model_function <- function(b, d) {
    1 - exp(-b * d)
  }
  i = 1
  # 对每一行进行拟合
  for (i in 1:n_rows) {
    # 获取当前行的数据
    d <- seq_along(matrix_data[i, ]) # d 的索引
    y <- matrix_data[i, ]            # 实际数据
    
    # 初始值设定
    start_b <- 1
    
    # 非线性最小二乘拟合
    fit <- try(nls(y ~ 1 - exp(-b * d), start = list(b = start_b)), silent = TRUE)
    
    # 如果拟合失败，跳过
    if (inherits(fit, "try-error")) {
      warning(paste("Fitting failed for row", i))
      next
    }
    
    # 获取拟合参数
    fitted_b <- coef(fit)["b"]
    
    # 计算拟合值
    fitted_y <- model_function(fitted_b, d)
    
    # 绘图
    plot(d, y, type = "p", pch = 16, col = "black", 
         main = paste("Row", i), xlab = "d", ylab = "y")
    lines(d, fitted_y, col = "red", lwd = 2) # 绘制拟合曲线
  }
  
  # 恢复绘图布局
  par(mfrow = c(1, 1))
}

# 示例数据
set.seed(123)
example_matrix <- matrix(runif(100, 0, 1), nrow = 10)

# 调用函数
fit_and_plot(as.matrix(matrix_data[c(20:40),]))

plot_exponential_fits <- function(data_matrix) {
  # Loop through each row of the input matrix
  for (i in 1:nrow(data_matrix)) {
    # Extract the data for the current row
    row_data <- data_matrix[i, ]
    
    # Define the exponential function
    exp_func <- function(d, b) {
      1 - exp(-b * d)
    }
    # Perform non-linear least squares fit
    fit_result <- nls(row_data ~ exp_func(seq_along(row_data), b), start = list(b = 0.1))
    
    # Extract the fitted parameter
    b_hat <- coef(fit_result)
    
    # Create a plot for the current row
    plot(seq_along(row_data), row_data, type = "b", col = "black", 
         main = paste("Row", i), xlab = "d", ylab = "1 - exp(-b * d)")
    
    # Plot the fitted exponential curve in red
    curve(exp_func(x, b_hat), add = TRUE, col = "red")
    
    # Add a legend
    legend("bottomright", legend = c("Actual Data", "Fitted Curve"), 
           col = c("black", "red"), lty = 1)
  }
}
data_matrix <- as.matrix(matrix_data[c(30:40),])
row_data <- data_matrix[1,]
plot_exponential_fits(as.matrix(matrix_data[c(30:40),]))
