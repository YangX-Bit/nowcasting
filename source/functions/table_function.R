nowcasts_table <- function(results_list, 
                           D = NULL,
                           report_unit = "week",
                           models_to_run = c("fixed_q", "fixed_b", "b_poly", "b_spline")) {
  # Basic checks
  if (is.null(D)) {
    stop("Parameter 'D' must be provided.")
  }
  if (!report_unit %in% c("week", "day")) {
    stop("report_unit must be 'week' or 'day'.")
  }
  
  library(lubridate)
  library(dplyr)
  
  # Decide factor for date shifting
  factor_loc <- if (report_unit == "week") 7 else 1
  
  # Prepare output list of data frames
  nowcasts_out <- list()
  
  n_runs <- length(results_list$case_true)  # how many sets of data we have
  
  for (i in seq_len(n_runs)) {
    # Extract data
    case_true     <- results_list[["case_true"]][[i]]
    case_reported <- results_list[["case_reported"]][[i]]
    dates         <- results_list[["dates"]][[i]]
    
    # Basic date references
    now      <- as.Date(dplyr::last(dates))
    earliest <- as.Date(dplyr::first(dates))
    last_date_for_delay <- now - days(D * factor_loc)
    
    # Initialize a data frame for storing nowcasts
    nowcasts_df <- data.frame(
      date          = dates,
      case_true     = case_true,
      case_reported = case_reported,
      row.names     = seq_along(dates)
    )
    
    # Store these references as columns (the same for all rows in this i)
    nowcasts_df$now                = now
    nowcasts_df$earliest           = earliest
    nowcasts_df$last_date_for_delay = last_date_for_delay
    
    # Dynamically add model results
    for (model_name in models_to_run) {
      # model_name might be "fixed_q", "fixed_b", ...
      samples <- results_list[[model_name]][[i]]$draws(variables = "N_t", format = "draws_matrix")
      nowcasts_df[[paste0("mean_", model_name)]]  <- apply(samples, 2, mean)
      nowcasts_df[[paste0("lower_", model_name)]] <- apply(samples, 2, quantile, probs = 0.025)
      nowcasts_df[[paste0("upper_", model_name)]] <- apply(samples, 2, quantile, probs = 0.975)
    }
    
    # Save to output list
    nowcasts_out[[i]] <- nowcasts_df
  }
  
  return(nowcasts_out)
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


highlight_metrics <- function(tables, method_names = NULL, date_labels = NULL) {
  
  # If no date labels are provided, use default numeric labels
  if (is.null(date_labels)) {
    date_labels <- paste0("Scenario ", 1:4)
  }
  
  # If no method names are provided, use default names
  if (is.null(method_names)) {
    method_names <- c("Method 1", "Method 2", "Method 3", "Method 4", "Method 5")
  }
  
  # Combine tables and add scenario column and custom method names
  combined_df <- do.call(rbind, lapply(seq_along(tables), function(i) {
    df <- tables[[i]]
    df$Scenario <- i
    # Replace the Method column with custom method names
    df$Method <- method_names[1:nrow(df)]
    return(df)
  }))
  
  # Define columns to be highlighted
  cols_to_process <- c("RMSE", "RMSPE", "MAE", "MAPE", "Interval_Width")
  
  # Process each column
  highlighted_df <- combined_df %>%
    group_by(Scenario) %>%
    mutate(across(all_of(cols_to_process), 
                  ~ ifelse(. == min(.), 
                           paste0("\\textcolor{red}{", formatC(., format = "f", digits = 5), "}"),
                           formatC(., format = "f", digits = 5)))) %>%
    # Special handling for Coverage Rate (highlight the maximum value)
    mutate(`Coverage_Rate` = 
             ifelse(`Coverage_Rate` == max(`Coverage_Rate`), 
                    paste0("\\textcolor{red}{", formatC(`Coverage_Rate`, format = "f", digits = 5), "}"),
                    formatC(`Coverage_Rate`, format = "f", digits = 5)))
  
  # Generate LaTeX table
  latex_table <- "\\begin{table}[htbp]\n\\centering\n"
  latex_table <- paste0(latex_table, "\\caption{Metrics Comparison}\n")
  latex_table <- paste0(latex_table, "\\begin{tabular}{c|c|c|c|c|c|c|c}\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  latex_table <- paste0(latex_table, "Scenario & RMSE & RMSPE & MAE & MAPE & Interval Width & Coverage Rate & Method \\\\\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  
  # Add data for each scenario
  for(i in seq_along(unique(highlighted_df$Scenario))) {
    scenario <- unique(highlighted_df$Scenario)[i]
    scenario_data <- highlighted_df[highlighted_df$Scenario == scenario,]
    
    # Add date row
    latex_table <- paste0(latex_table, "\\multicolumn{8}{l}{\\textbf{Now is ", date_labels[i], "}} \\\\\n")
    latex_table <- paste0(latex_table, "\\hline\n")
    
    # Add data rows for the current scenario
    for(j in 1:nrow(scenario_data)) {
      row <- scenario_data[j,]
      row_str <- paste(
        "", # First column left blank
        row$RMSE,
        row$RMSPE,
        row$MAE,
        row$MAPE,
        row$`Interval_Width`,
        row$`Coverage_Rate`,
        row$Method,
        sep = " & "
      )
      latex_table <- paste0(latex_table, row_str, " \\\\\n")
    }
    
    # Add a separator line between scenarios
    latex_table <- paste0(latex_table, "\\hline\n")
  }
  
  latex_table <- paste0(latex_table, "\\end{tabular}\n")
  latex_table <- paste0(latex_table, "\\end{table}")
  
  return(latex_table)
}