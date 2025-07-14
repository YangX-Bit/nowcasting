nowcasts_table <- function(results_list, 
                           D = NULL,
                           report_unit = "week",
                           methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
                           replicate_id = NA_integer_, alpha = 0.05) {
  # Basic checks
  if (is.null(D)) stop("Parameter 'D' must be provided.")
  if (!report_unit %in% c("week", "day")) {
    stop("report_unit must be 'week' or 'day'.")
  }
  
  library(lubridate)
  library(dplyr)
  
  # Decide factor for date shifting
  factor_loc <- if (report_unit == "week") 7 else 1
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
    nowcasts_df$now                <- now
    nowcasts_df$earliest           <- earliest
    nowcasts_df$last_date_for_delay <- last_date_for_delay
    if (!is.na(replicate_id)) {
      nowcasts_df$replicate_id       <- replicate_id # for replicate
    }
    
    lower_p <- alpha / 2
    upper_p <- 1 - alpha / 2
    
    # Dynamically add model results
    for (model_name in methods) {
      # model_name might be "fixed_q", "fixed_b", ...
      samples <- results_list[[model_name]][[i]]$draws(variables = "N", format = "draws_matrix")
      nowcasts_df[[paste0("mean_", model_name)]]  <- apply(samples, 2, mean)
      nowcasts_df[[paste0("lower_", model_name)]] <- apply(samples, 2, quantile, probs = lower_p, na.rm = TRUE)# if this na.rm = TRUE is okay?
      nowcasts_df[[paste0("upper_", model_name)]] <- apply(samples, 2, quantile, probs = upper_p, na.rm = TRUE)
    }
    
    # Save to output list
    nowcasts_out[[i]] <- nowcasts_df
  }
  
  return(nowcasts_out)
}


# measurement metrics
calculate_metrics <- function(df, methods = c("fixed_q", "fixed_b", "linear_b", "ou_b")) {
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
  
  for (method in methods) {
    mean_col  <- paste0("mean_", method)
    lower_col <- paste0("lower_", method)
    upper_col <- paste0("upper_", method)
    
    error <- df[[mean_col]] - df$case_true
    abs_error <- abs(error)
    
    # For RMSPE & MAPE, use a small epsilon where case_true == 0
    nonzero_idx <- df$case_true != 0
    
    rel_error <- error[nonzero_idx] / df$case_true[nonzero_idx]
    abs_rel_error <- abs(rel_error)
    
    # error <- df[[mean_col]] - df$case_true
    # abs_error <- abs(error)
    # rel_error <- error / df$case_true
    # abs_rel_error <- abs(rel_error)
    
    
    
    # Metrics
    RMSE  <- sqrt(mean(error^2, na.rm = TRUE))
    RMSPE <- sqrt(mean(rel_error^2, na.rm = TRUE)) * 100
    MAE   <- mean(abs_error, na.rm = TRUE)
    MAPE  <- mean(abs_rel_error, na.rm = TRUE) * 100
    
    # Prediction interval metrics
    interval_width <- mean(df[[upper_col]] - df[[lower_col]], na.rm = TRUE)
    coverage <- mean((df$case_true >= df[[lower_col]]) & (df$case_true <= df[[upper_col]]), na.rm = TRUE)
    
    results[[method]] <- c(
      RMSE = RMSE,
      RMSPE = RMSPE,
      MAE = MAE,
      MAPE = MAPE,
      Interval_Width = interval_width,
      Coverage_Rate = coverage
    )
  }
  
  # Convert to data frame
  results_df <- do.call(rbind, results)
  results_df <- as.data.frame(round(results_df,2))
  results_df$Method <- rownames(results_df)
  rownames(results_df) <- NULL
  
  return(results_df)
}



# --------------------------------------------------
# compute_all_nowcasts_tables()
# --------------------------------------------------
# out_list:  a list of length num_sims
#            e.g. out_list_1_NFR, where out_list[[i]]
#            has the data needed for nowcasts_table()
# D:         max delay
# report_unit: "day" or "week"
# methods: c("fixed_q","fixed_b","linear_b","ou_b"), etc.
# replicate_ids: a vector of replicate IDs, e.g. 1:num_sims
#
# Return: a list of length num_sims, 
#   each element is the *list* returned by nowcasts_table()
#   i.e. results_list[[i]] is itself a list of length T
#   (one data frame per time window).
# --------------------------------------------------

compute_all_nowcasts_tables <- function(
    out_list,
    D,
    report_unit = "day",
    methods = c("fixed_q", "fixed_b", "rw_b", "ou_b"),
    replicate_ids = 1:length(out_list)
){
  num_sims <- length(out_list)
  results_all <- vector("list", length = num_sims)
  
  for (i in seq_len(num_sims)) {
    # Call nowcasts_table for each replicate
    results_all[[i]] <- nowcasts_table(
      results_list  = out_list[[i]],  # this replicate's data
      D             = D,
      report_unit   = report_unit,
      methods = methods,
      replicate_id  = replicate_ids[i]   # store replicate ID
    )
  }
  return(results_all)
}


# --------------------------------------------------
# average_nowcasts_tables()
# --------------------------------------------------
# results_all: the object returned by compute_all_nowcasts_tables().
#   A list of length num_sims, each element is a list of T data frames.
# numeric_cols: optional vector of columns you want to average.
#   If NULL, we auto-detect columns that start with mean_, lower_, upper_,
#   plus case_true, case_reported, etc.
#
# Return: a list of length T, each element is a single data frame
#   with averaged columns across replicates.
# --------------------------------------------------

average_nowcasts_tables <- function(results_all, numeric_cols = NULL) {
  num_sims <- length(results_all)
  # Assume each replicate has the same number of windows T
  T <- length(results_all[[1]])  
  
  # Auto-detect numeric columns if not provided
  if (is.null(numeric_cols)) {
    df_example <- results_all[[1]][[1]]
    typical_prefixes <- c("mean_","lower_","upper_","case_true","case_reported")
    all_cols <- names(df_example)
    numeric_cols <- all_cols[sapply(all_cols, function(x) {
      startsWith(x, "mean_") ||
        startsWith(x, "lower_") ||
        startsWith(x, "upper_") ||
        x %in% c("case_true", "case_reported")
    })]
  }
  
  # Prepare output: a list of length T
  results_avg <- vector("list", length = T)
  
  for (t in seq_len(T)) {
    # Gather data frames from each replicate for window t
    df_list_t <- lapply(seq_len(num_sims), function(i) {
      results_all[[i]][[t]]
    })
    # Use first replicate as a template
    df_avg_t <- df_list_t[[1]]
    
    # For each numeric column, compute row-wise mean
    for (colname in numeric_cols) {
      vals <- sapply(df_list_t, function(df) df[[colname]])
      df_avg_t[[colname]] <- rowMeans(vals, na.rm = TRUE)
    }
    
    # Remove replicate_id if present
    if ("replicate_id" %in% names(df_avg_t)) {
      df_avg_t$replicate_id <- NULL
    }
    
    results_avg[[t]] <- df_avg_t
  }
  
  return(results_avg)
}


# --------------------------------------------------
# average_nowcasts_metrics()
# --------------------------------------------------
# results_all: the same object returned by compute_all_nowcasts_tables().
# methods: which models to compute metrics for, e.g. c("fixed_q","ou_b")
#
# Returns a list of length T, each is a data frame of averaged metrics.
# e.g. metrics_avg_t has columns: Method, RMSE, MAPE, Coverage_Rate, etc.
# --------------------------------------------------

average_nowcasts_metrics <- function(
    results_all,
    methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
    filter_length = NULL  # e.g., D, to use last D rows
) {
  library(dplyr)
  
  num_sims <- length(results_all)
  T <- length(results_all[[1]])  # number of windows
  
  # We'll store the final average metrics for each window t
  metrics_t_averaged <- vector("list", length = T)
  
  for (t in seq_len(T)) {
    # Collect the metrics from all replicates for this window
    metrics_for_t_list <- list()
    
    for (r in seq_len(num_sims)) {
      df_rt <- results_all[[r]][[t]]  # the data frame for replicate r, window t
      
      # If we want to filter the last 'filter_length' rows
      if (!is.null(filter_length) && filter_length > 0) {
        n_rows <- nrow(df_rt)
        start_row <- max(1, n_rows - filter_length + 1)
        df_rt <- df_rt[start_row:n_rows, , drop = FALSE]
      }
      
      # Compute metrics on the filtered or full data
      metrics_rt <- calculate_metrics(df_rt, methods = methods)
      metrics_for_t_list[[r]] <- metrics_rt
    }
    
    # Combine row-wise (all replicates for this t)
    metrics_for_t <- bind_rows(metrics_for_t_list)
    
    # Convert 'Method' into a factor with the specified order:
    metrics_for_t <- metrics_for_t %>%
      mutate(Method = factor(Method, levels = methods))
    
    # Average across replicates, keeping the factor order
    metrics_avg_t <- metrics_for_t %>%
      group_by(Method) %>%
      summarize(
        RMSE           = mean(RMSE, na.rm=TRUE),
        RMSPE          = mean(RMSPE, na.rm=TRUE),
        MAE            = mean(MAE, na.rm=TRUE),
        MAPE           = mean(MAPE, na.rm=TRUE),
        Interval_Width = mean(Interval_Width, na.rm=TRUE),
        Coverage_Rate  = mean(Coverage_Rate, na.rm=TRUE),
        .groups = "drop"
      ) %>%
      mutate_if(is.numeric, round, 2) %>%
      arrange(Method)  # ensure final ordering follows factor levels
    
    metrics_t_averaged[[t]] <- metrics_avg_t
  }
  
  return(metrics_t_averaged)
}


highlight_metrics <- function(tables, 
                              model_names = NULL, 
                              date_labels = NULL, 
                              scenario_names = NULL,
                              digits = 2, 
                              table_caption = "Metrics comparison", 
                              report_unit = "day", 
                              D = NULL, 
                              first_date = NULL,
                              report_scenario = "fr",       # e.g., "fr" for fully reported, "nfr" for not fully reported
                              model_origin = "constant"     # e.g., "constant", "random_walks", "ou_processes"
) {
  # Check required parameters
  if (is.null(date_labels)) {
    stop("Please provide date_labels (a vector of dates).")
  }
  if (is.null(model_names)) {
    stop("Please provide model_names (a vector of model names) with length matching the number of rows in each table.")
  }
  if (is.null(scenario_names)) {
    scenario_names <- paste0("Scenario ", 1:length(tables))
  }
  if (is.null(D) || is.null(first_date)) {
    stop("Please provide both 'D' (report delay) and 'first_date'.")
  }
  
  # Convert first_date and date_labels to Date objects
  first_date <- as.Date(first_date)
  date_labels <- as.Date(date_labels)
  
  # Determine conversion factor based on report_unit (default is day)
  factor_loc <- if (report_unit == "week") 7 else 1
  
  # Define unit label and header text for time column
  unit_label <- if (report_unit == "week") "week" else "day"
  time_header <- if (report_unit == "week") "Time (in weeks)" else "Time (in days)"
  
  # Combine data from each table and add scenario numbers and model names
  combined_df <- do.call(rbind, lapply(seq_along(tables), function(i) {
    df <- tables[[i]]
    df$Scenario <- i
    df$Method <- model_names[1:nrow(df)]
    return(df)
  }))
  
  # Columns to be highlighted
  cols_to_process <- c("RMSE", "RMSPE", "MAE", "MAPE", "Interval_Width")
  
  highlighted_df <- combined_df %>%
    group_by(Scenario) %>%
    mutate(across(all_of(cols_to_process), 
                  ~ ifelse(. == min(.), 
                           paste0("\\textcolor{red}{", formatC(., format = "f", digits = digits), "}"),
                           formatC(., format = "f", digits = digits)))) %>%
    mutate(Coverage_Rate = 
             ifelse(Coverage_Rate == max(Coverage_Rate), 
                    paste0("\\textcolor{red}{", formatC(Coverage_Rate, format = "f", digits = digits), "}"),
                    formatC(Coverage_Rate, format = "f", digits = digits))) %>%
    ungroup()
  
  # Construct LaTeX table string
  latex_table <- "\\begin{table}[H]\n\\footnotesize\n\\centering\n"
  # Use auto-adjust width (left aligned) for the first column and update header with time unit
  latex_table <- paste0(latex_table, "\\begin{tabular}{l|c|c|c|c|c|c|l}\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  latex_table <- paste0(latex_table, time_header, " & RMSE & RMSPE & MAE & MAPE & MCIW & MCR & Method \\\\\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  
  # Output 4 rows of data for each scenario
  for (i in seq_along(unique(highlighted_df$Scenario))) {
    scenario <- unique(highlighted_df$Scenario)[i]
    scenario_data <- highlighted_df[highlighted_df$Scenario == scenario, ]
    
    # Calculate total time units (days or weeks)
    total_units <- as.numeric((date_labels[i] - first_date) / factor_loc) + 1
    
    # Calculate completed and incomplete time units
    if (total_units <= D) {
      a <- 0
      b <- total_units
    } else {
      b <- D
      a <- total_units - b
    }
    total_units <- floor(total_units)
    a <- floor(a); b <- floor(b)
    
    # First row: Use scenario_names with numbering (e.g., "(1) Initial")
    row1 <- scenario_data[1, ]
    label1 <- paste0("(", i, ") ", scenario_names[i])
    row1_str <- paste(label1,
                      row1$RMSE,
                      row1$RMSPE,
                      row1$MAE,
                      row1$MAPE,
                      row1$Interval_Width,
                      row1$Coverage_Rate,
                      row1$Method,
                      sep = " & ")
    latex_table <- paste0(latex_table, row1_str, " \\\\\n")
    
    # Second row: Now time (total units)
    row2 <- scenario_data[2, ]
    label2 <- paste0("Now: $t=", total_units, "$")
    row2_str <- paste(label2,
                      row2$RMSE,
                      row2$RMSPE,
                      row2$MAE,
                      row2$MAPE,
                      row2$Interval_Width,
                      row2$Coverage_Rate,
                      row2$Method,
                      sep = " & ")
    latex_table <- paste0(latex_table, row2_str, " \\\\\n")
    
    # Third row: Completed units
    row3 <- scenario_data[3, ]
    label3 <- paste0("Completed ", unit_label, "s: ", a)
    row3_str <- paste(label3,
                      row3$RMSE,
                      row3$RMSPE,
                      row3$MAE,
                      row3$MAPE,
                      row3$Interval_Width,
                      row3$Coverage_Rate,
                      row3$Method,
                      sep = " & ")
    latex_table <- paste0(latex_table, row3_str, " \\\\\n")
    
    # Fourth row: Incomplete units
    row4 <- scenario_data[4, ]
    label4 <- paste0("Incomplete ", unit_label, "s: ", b)
    row4_str <- paste(label4,
                      row4$RMSE,
                      row4$RMSPE,
                      row4$MAE,
                      row4$MAPE,
                      row4$Interval_Width,
                      row4$Coverage_Rate,
                      row4$Method,
                      sep = " & ")
    latex_table <- paste0(latex_table, row4_str, " \\\\\n")
    
    latex_table <- paste0(latex_table, "\\hline\n")
  }
  
  latex_table <- paste0(latex_table, "\\end{tabular}\n")
  latex_table <- paste0(latex_table, "\\caption{", table_caption, "}\n")
  # Create label based on report_scenario and model_origin parameters
  latex_table <- paste0(latex_table, "\\label{tab:", report_scenario, "_", model_origin, "}\n")
  latex_table <- paste0(latex_table, "\\end{table}")
  
  return(cat(latex_table))
}


# highlight_metrics <- function(tables, model_names = NULL, date_labels = NULL, digits = 2, 
#                               table_caption = "Metrics comparison", 
#                               report_unit = "day", D = NULL, first_date = NULL) {
#   # If no date labels are provided, use default numeric labels
#   if (is.null(date_labels)) {
#     date_labels <- paste0("Scenario ", 1:4)
#   }
#   
#   # If no method names are provided, use default names
#   if (is.null(model_names)) {
#     model_names <- c("Method 1", "Method 2", "Method 3", "Method 4", "Method 5")
#   }
#   
#   # Ensure D and first_date are provided
#   if (is.null(D) || is.null(first_date)) {
#     stop("Both 'D' (report delay) and 'first_date' must be provided.")
#   }
#   
#   # Convert first_date and date_labels to Date objects
#   first_date <- as.Date(first_date)
#   date_labels <- as.Date(date_labels)
#   
#   # Factor to convert report_unit to days
#   factor_loc <- if (report_unit == "week") 7 else 1
#   
#   # Combine tables and add scenario column and custom method names
#   combined_df <- do.call(rbind, lapply(seq_along(tables), function(i) {
#     df <- tables[[i]]
#     df$Scenario <- i
#     df$Method <- model_names[1:nrow(df)]
#     return(df)
#   }))
#   
#   # Define columns to be highlighted
#   cols_to_process <- c("RMSE", "RMSPE", "MAE", "MAPE", "Interval_Width")
#   
#   # Process each column
#   highlighted_df <- combined_df %>%
#     group_by(Scenario) %>%
#     mutate(across(all_of(cols_to_process), 
#                   ~ ifelse(. == min(.), 
#                            paste0("\\textcolor{red}{", formatC(., format = "f", digits = digits), "}"),
#                            formatC(., format = "f", digits = digits)))) %>%
#     mutate(`Coverage_Rate` = 
#              ifelse(`Coverage_Rate` == max(`Coverage_Rate`), 
#                     paste0("\\textcolor{red}{", formatC(`Coverage_Rate`, format = "f", digits = digits), "}"),
#                     formatC(`Coverage_Rate`, format = "f", digits = digits)))
#   
#   # Generate LaTeX table
#   latex_table <- "\\begin{table}[htbp]\n\\centering\n"
#   latex_table <- paste0(latex_table, "\\begin{tabular}{c|c|c|c|c|c|c|c}\n")
#   latex_table <- paste0(latex_table, "\\hline\n")
#   latex_table <- paste0(latex_table, "Scenario & RMSE & RMSPE & MAE & MAPE & Interval Width & Coverage Rate & Method \\\\\n")
#   latex_table <- paste0(latex_table, "\\hline\n")
#   
#   # Add data for each scenario
#   for (i in seq_along(unique(highlighted_df$Scenario))) {
#     scenario <- unique(highlighted_df$Scenario)[i]
#     scenario_data <- highlighted_df[highlighted_df$Scenario == scenario,]
#     
#     # Calculate total days of data used
#     total_days <- as.numeric((date_labels[i] - first_date) / factor_loc) + 1 # Convert to correct unit
#     
#     # Compute a (completed) and b (incomplete) days
#     if (total_days <= D) {
#       a <- 0
#       b <- total_days
#     } else {
#       b <- D
#       a <- total_days - b
#     }
#     
#     # Add date row with additional reporting information
#     latex_table <- paste0(latex_table, "\\multicolumn{8}{l}{\\textbf{Now is ", date_labels[i], 
#                           ", with ", a, " ", report_unit, " completed, ", b, " ", report_unit, " incomplete.}} \\\\\n")
#     latex_table <- paste0(latex_table, "\\hline\n")
#     
#     # Add data rows for the current scenario
#     for (j in 1:nrow(scenario_data)) {
#       row <- scenario_data[j,]
#       row_str <- paste(
#         "", # First column left blank
#         row$RMSE,
#         row$RMSPE,
#         row$MAE,
#         row$MAPE,
#         row$`Interval_Width`,
#         row$`Coverage_Rate`,
#         row$Method,
#         sep = " & "
#       )
#       latex_table <- paste0(latex_table, row_str, " \\\\\n")
#     }
#     
#     # Add a separator line between scenarios
#     latex_table <- paste0(latex_table, "\\hline\n")
#   }
#   
#   latex_table <- paste0(latex_table, "\\end{tabular}\n")
#   latex_table <- paste0(latex_table, "\\caption{", table_caption,"}\n")
#   latex_table <- paste0(latex_table, "\\end{table}")
#   
#   return(cat(latex_table))
# }