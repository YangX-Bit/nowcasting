library(ggplot2)
library(gridExtra)
library(ggforce)

# nowcasts_plot <- function(results_list, models_to_run = c("fixed_q", "fixed_b", "b_poly", "b_spline"),
#                           title = NULL, x_lab = NULL, y_lab = "Cases / Nowcast") {
#   #
#   case_true <- results_list[["case_true"]][[i]]
#   case_reported <- results_list[["case_reported"]][[i]]
#   dates <- results_list[["dates"]][[i]]
#   
#   # Initialize an empty data frame for nowcast data
#   nowcasts <- data.frame(date = dates,
#                          case_true = case_true,
#                          case_reported = case_reported)
#   
#   # Dynamically add model results
#   for (model_name in names(models_to_run)) {
#     samples <- results_list[[model_name]]
#     nowcasts[[paste0("mean_", model_name)]] <- apply(samples, 2, mean)
#     nowcasts[[paste0("lower_", model_name)]] <- apply(samples, 2, quantile, probs = 0.025)
#     nowcasts[[paste0("upper_", model_name)]] <- apply(samples, 2, quantile, probs = 0.975)
#   }
#   
#   # Create the base ggplot object
#   p <- ggplot(nowcasts, aes(x = date)) +
#     geom_line(aes(y = case_true, color = "Real Cases")) +
#     geom_line(aes(y = case_reported, color = "Reported Cases"))
#   
#   # Dynamically add ribbons and lines for each model
#   for (model_name in names(models_to_run)) {
#     p <- p +
#       geom_ribbon(aes_string(ymin = paste0("lower_", model_name),
#                              ymax = paste0("upper_", model_name)),
#                   fill = model_name, alpha = 0.3) +  # Use model name as fill color
#       geom_line(aes_string(y = paste0("mean_", model_name), color = model_name))
#   }
#   
#   # Add manual color scale
#   model_colors <- c("Real Cases" = "red", "Reported Cases" = "black",
#                     setNames(rainbow(length(models_to_run)), names(models_to_run)))  # Dynamic colors
#   p <- p +
#     scale_color_manual(values = model_colors) +
#     labs(title = title,
#          x = x_lab,
#          y = y_lab,
#          color = NULL) +
#     theme_minimal() +
#     theme(
#       legend.position = c(0.1, 0.9),  # Legend position
#       legend.justification = c(0, 1),  # Top-left alignment
#       legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"), # Legend border
#       legend.key = element_rect(fill = "white", color = NA),
#       legend.text = element_text(size = 16),
#       legend.title = element_text(size = 16),
#       axis.text = element_text(size = 16),
#       axis.title = element_text(size = 16)
#     )
#   p
#   return(p)
# }

# check the q shape and output plots of fit
fit_exp_plot <- function(matrix_data, ncol = 3, nrow = 3, pages = 1) {
  if (!is.matrix(matrix_data)) stop("Input must be a matrix.")
  
  # Normalize matrix columns
  matrix_data <- normalize_matrix_columns(matrix_data)
  
  # get date for the row
  row_names <- rownames(matrix_data)
  if (is.null(row_names)) {
    row_names <- as.character(1:nrow(matrix_data))
  }
  
  n_rows <- nrow(matrix_data); D <- ncol(matrix_data) - 1
  coef_saved <- data.frame(b = as.numeric(rep(0, n_rows)))
  
  for (i in 1:n_rows) {
    data_fit <- data.frame(
      x = c(0:D),
      y = as.numeric(matrix_data[i, ])
    )
    tryCatch({
      model_fit <- nls(y ~ (1 - exp(-b * x)), data = data_fit, start = list(b = 0.5))
      coef_saved[i, 1] <- coef(model_fit)["b"]
    }, error = function(e) {
      warning(paste("Fitting failed for row", i, ":", e$message))
    })
  }
  
  # Prepare data for ggplot
  x_vals <- c(0:D)
  plot_data <- data.frame()
  for (i in 1:n_rows) {
    y_vals <- (1 - exp(-coef_saved$b[i] * x_vals))
    temp_data <- data.frame(
      x = x_vals,
      y = as.numeric(matrix_data[i, ]),
      fit = y_vals,
      Row = factor(rep(row_names[i], length(x_vals))) 
    )
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Output multiple pages
  plots <- list()
  for (page in pages) {
    p <- ggplot(plot_data, aes(x = x)) +
      geom_line(aes(y = y), color = "black") +
      geom_line(aes(y = fit), color = "red", linetype = "dashed", size = 1) +
      facet_wrap_paginate(~ Row, ncol = ncol, nrow = nrow, page = page) +
      labs(title = paste("Fitted Plots (Page", page, ")"), x = NULL, y = NULL) +
      theme_minimal()
    plots[[page]] <- p
  }
  
  return(plots)
}



nowcasts_plot <- function(results_list, D = NULL, report_unit = "week",
                          models_to_run = c("fixed_q", "fixed_b", "b_poly", "b_spline"),
                          title = NULL, x_lab = NULL, y_lab = "Cases / Nowcast") {
  
  p_out <- list(); nowcasts_out <- list();
  
  if(report_unit == "week"){ factor_loc = 7
  }else if(report_unit == "day"){ factor_loc = 1}
  else{ stop("Wrong input of parameter report_unit. It has to be 'week' or 'day'.")}
  for (i in 1:length(results_list$case_true)) {
    # Extract data for the specific time window
    case_true <- results_list[["case_true"]][[i]]
    case_reported <- results_list[["case_reported"]][[i]]
    dates <- results_list[["dates"]][[i]]
    
    now <- as.Date(last(dates))
    last_date_for_delay <- now - days(D) * factor_loc
    
    # Initialize an empty data frame for nowcast data
    nowcasts <- data.frame(date = dates,
                           case_true = case_true,
                           case_reported = case_reported,
                           row.names = c(1:length(dates)))
    
    # Dynamically add model results
    for (model_name in models_to_run) {
      # Extract samples
      samples <- results_list[[model_name]][[i]]$draws(variables = "lambda_t", format = "draws_matrix") 
      nowcasts[[paste0("mean_", model_name)]] <- apply(samples, 2, mean)
      nowcasts[[paste0("lower_", model_name)]] <- apply(samples, 2, quantile, probs = 0.025)
      nowcasts[[paste0("upper_", model_name)]] <- apply(samples, 2, quantile, probs = 0.975)
    }
    #now <- max(nowcasts$dates)
    nowcasts_out[[i]] <- nowcasts
    
    model_colors <- c(
      "Real Cases" = "red", 
      "Reported Cases" = "black",
      "fixed_q" = "#228B22",   # 
      "fixed_b" = "#ffff66",  # 
      "b_poly" = "#8446c6",  # 
      "b_spline" = "#4682B4" # 
    )
    # Create the base ggplot object
    p <- ggplot(nowcasts, aes(x = date)) +
      geom_line(aes(y = case_true, color = "Real Cases")) +
      geom_line(aes(y = case_reported, color = "Reported Cases")) +
      geom_vline(xintercept = last_date_for_delay, color = "orange", linetype = "dashed", size = 1) +
      annotate("text", x = last_date_for_delay, y = -1, 
               label = last_date_for_delay, vjust = 2, color = "orange") +
      annotate("text", x = now, y = -1, 
               label = paste0("now: ", now), hjust = 1,vjust = 2, color = "red") +
      geom_point(data = data.frame(x = now, y = 0), aes(x = x, y = y, shape = "Today"), 
                 size = 2, color = "red")+
      scale_shape_manual(values = c("Today" = 17)) 
    
    # Dynamically add ribbons and lines for each model
    for (model_name in models_to_run) {
      p <- p +
        geom_ribbon(aes_string(ymin = paste0("lower_", model_name),
                               ymax = paste0("upper_", model_name)),
                    fill = model_colors[model_name], alpha = 0.3) +  # Use model name as fill color
        geom_line(
          aes_string(y = paste0("mean_", model_name)), 
          color = model_colors[model_name]  # pre-defined color
        )
    }
    # Add manual color scale
    model_colors <- c("Real Cases" = "red", "Reported Cases" = "black",
                      setNames(rainbow(length(models_to_run)), models_to_run))  # Dynamic colors
    # legend_position <- adjust_legend_position(nowcasts, "date", "case_true")
    
    p <- p +
      scale_color_manual(values = model_colors) +
      labs(title = title,
           x = x_lab,
           y = y_lab,
           color = NULL) +
      theme_minimal() +
      theme(
        legend.position = c(0.7, 0.9),  # Legend position
        legend.justification = c(0, 1),  # Top-left alignment
        legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"), # Legend border
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
      )
    p_out[[i]] <- p
  }
  return(list(plots = p_out,
              nowcasts = nowcasts_out))
}

# Function to dynamically determine legend position
# adjust_legend_position <- function(data, x_var, y_var) {
#   x_mid <- mean(range(data[[x_var]], na.rm = TRUE))  # mid value
#   y_max <- max(data[[y_var]], na.rm = TRUE)          # max value
#   return(c(0.9, y_max * 0.9))                        # right-top
# }
