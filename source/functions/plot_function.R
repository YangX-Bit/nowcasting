library(ggplot2)
library(gridExtra)
library(ggforce)
library(lubridate)
library(dplyr)

# for showing the relastionship between bt and sigma
simulate_plot_q <- function(b_inits,
                            sigma_rws,
                            T = 30,
                            D = 5,
                            delta = 0,
                            n_sims = 100) {
  
  # Basic checks
  if (length(b_inits) != length(sigma_rws)) {
    stop("Length of 'b_inits' and 'sigma_rws' must be the same.")
  }
  if (delta < 0) {
    stop("'delta' must be a non-negative number.")
  }
  
  # Prepare a list to store averaged q_{(t,d)} for each parameter set
  # Each element is a T x (D+1) matrix
  avg_q_list <- vector("list", length(b_inits))
  
  # Color palette for multiple curves
  my_cols <- rainbow(length(b_inits))
  
  # Set up an empty plot
  plot.new()
  plot.window(xlim = c(0, D), ylim = c(0, 1))
  axis(1); axis(2)
  box()
  # title(main = "Average q(t,d) Curves", xlab = "Delay d", ylab = expression(q[(t,d)]))
  
  # Add a horizontal line y=1 in red
  abline(h = 1, col = "red", lty = 2, lwd = 2)
  
  for (i in seq_along(b_inits)) {
    # Temporary array for storing q(t,d) from n_sims simulations
    # Dimensions: T x (D+1) x n_sims
    q_mat_sims <- array(0, dim = c(T, D + 1, n_sims))
    
    for (s in seq_len(n_sims)) {
      
      # 1. Generate the random-walk trajectory of b_t
      b_vec <- numeric(T)
      b_vec[1] <- b_inits[i]
      for (t in 2:T) {
        b_vec[t] <- max(b_vec[t - 1] + rnorm(1, 0, sigma_rws[i]), 0.01)
      }
      
      # 2. Compute q_{(t,d)} = 1 - exp(- b_t (d + delta)) for d = 0..D
      q_mat <- matrix(0, nrow = T, ncol = D + 1)
      for (t_idx in 1:T) {
        for (d_idx in 0:D) {
          q_mat[t_idx, d_idx + 1] <- 1 - exp(-b_vec[t_idx] * (d_idx + delta))
        }
      }
      
      q_mat_sims[, , s] <- q_mat
    }
    
    # 3. Average over n_sims to get a smoother curve
    avg_q <- apply(q_mat_sims, c(1, 2), mean)
    avg_q_list[[i]] <- avg_q
    
    # 4. For illustration, plot q_{(T,d)} against d
    lines(0:D, avg_q[T, ], col = my_cols[i], lwd = 2)
  }
  
  # Add legend
  legend(
    "bottomright",
    legend = paste0("b_init=", b_inits, ", sigma_rw=", sigma_rws),
    col = my_cols, lwd = 2
  )
  
  invisible(avg_q_list)
}



#check the q shape and output plots of fit
fit_exp_plot <- function(matrix_data, ncol = 3, nrow = 3, pages = 1, if_fit = T) {
  if (!is.matrix(matrix_data)) stop("Input must be a matrix.")

  # Normalize matrix columns
  matrix_data <- normalize_matrix_columns(matrix_data)

  # get date for the row
  row_names <- rownames(matrix_data)
  if (is.null(row_names)) {
    row_names <- as.character(1:nrow(matrix_data))
  }

  n_rows <- nrow(matrix_data); D <- ncol(matrix_data) - 1
  coef_saved <- data.frame(b = as.numeric(rep(0, n_rows)),
                           phi = as.numeric(rep(0, n_rows)))

  for (i in 1:n_rows) {
    data_fit <- data.frame(
      x = c(0:D),
      y = as.numeric(matrix_data[i, ])
    )
    tryCatch({
      model_fit <- nls(y ~ (1 - (1 - phi)*exp(-b * x)), data = data_fit, start = list(b = 0.2, phi = 0.2))
      coef_saved[i, 1] <- coef(model_fit)["b"]
      coef_saved[i, 2] <- coef(model_fit)["phi"]
    }, error = function(e) {
      warning(paste("Fitting failed for row", i, ":", e$message))
    })
  }

  # Prepare data for ggplot
  x_vals <- c(0:D)
  plot_data <- data.frame()
  for (i in 1:n_rows) {
    y_vals <- (1 - (1 - coef_saved$phi[i]) * exp(-coef_saved$b[i] * x_vals))
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
      geom_line(aes(y = y), color = "black")

    # Conditionally add the fit line based on if_fit parameter
    if (if_fit) {
      p <- p + geom_line(aes(y = fit), color = "red", linetype = "dashed", size = 1)
    }

    p <- p +
      facet_wrap_paginate(~ Row, ncol = ncol, nrow = nrow, page = page) +
      labs(title = paste("Fitted Plots (Page", page, ")"), x = NULL, y = NULL) +
      theme_minimal()

    plots[[page]] <- p
  }
  list_out <- list(plots = plots,
                   coefficients = coef_saved)
  return(list_out)
}


# fit_exp_plot <- function(matrix_data, ncol = 3, nrow = 3, pages = 1, if_fit = T) {
#   if (!is.matrix(matrix_data)) stop("Input must be a matrix.")
# 
#   # Normalize matrix columns
#   matrix_data <- normalize_matrix_columns(matrix_data)
# 
#   # get date for the row
#   row_names <- rownames(matrix_data)
#   if (is.null(row_names)) {
#     row_names <- as.character(1:nrow(matrix_data))
#   }
# 
#   n_rows <- nrow(matrix_data); D <- ncol(matrix_data) - 1
#   #coef_saved <- data.frame(b = as.numeric(rep(0, n_rows)), delta = as.numeric(rep(0, n_rows)))
#   coef_saved <- as.numeric(rep(0, n_rows))
# 
#   for (i in 1:n_rows) {
#     data_fit <- data.frame(
#       x = c(0:D),
#       y = as.numeric(matrix_data[i, ])
#     )
#     tryCatch({
#       #model_fit <- nls(y ~ (1 - exp(-b * (x + delta))), data = data_fit, start = list(b = 0.5, delta = 0.2))
#       model_fit <- nls(y ~ (1 - exp(-b * (x))), data = data_fit, start = list(b = 0.5))
#       coef_saved[i] <- as.numeric(coef(model_fit)["b"])
#       #coef_saved[i, 2] <- coef(model_fit)["delta"]
#     }, error = function(e) {
#       warning(paste("Fitting failed for row", i, ":", e$message))
#     })
#   }
# 
#   # Prepare data for ggplot
#   x_vals <- c(0:D)
#   plot_data <- data.frame()
#   for (i in 1:n_rows) {
#     #y_vals <- (1 - exp(-coef_saved$b[i] * (x_vals + coef_saved$delta[i])))
#     y_vals <- (1 - exp(-coef_saved[i] * (x_vals)))
#     temp_data <- data.frame(
#       x = x_vals,
#       y = as.numeric(matrix_data[i, ]),
#       fit = y_vals,
#       Row = factor(rep(row_names[i], length(x_vals)))
#     )
#     plot_data <- rbind(plot_data, temp_data)
#   }
# 
#   # Output multiple pages
#   plots <- list()
#   for (page in pages) {
#     p <- ggplot(plot_data, aes(x = x)) +
#       geom_line(aes(y = y), color = "black")
# 
#     # Conditionally add the fit line based on if_fit parameter
#     if (if_fit) {
#       p <- p + geom_line(aes(y = fit), color = "red", linetype = "dashed", size = 1)
#     }
# 
#     p <- p +
#       facet_wrap_paginate(~ Row, ncol = ncol, nrow = nrow, page = page) +
#       labs(title = paste("Fitted Plots (Page", page, ")"), x = NULL, y = NULL) +
#       theme_minimal()
# 
#     plots[[page]] <- p
#   }
#   list_out <- list(plots = plots,
#                    coefficients = coef_saved)
#   return(list_out)
# }

nowcasts_plot <- function(nowcasts_list,
                          D = NULL,
                          report_unit = "week",
                          methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
                          title = NULL,
                          x_lab = NULL,
                          y_lab = "Cases / Nowcast",
                          legend_position = NULL,
                          combine_plots = FALSE,
                          ncol = NULL,
                          nrow = NULL) {
  library(ggplot2)
  library(lubridate)
  library(dplyr)
  library(patchwork)
  
  # Basic checks
  if (is.null(D)) {
    stop("Parameter 'D' must be provided.")
  }
  if (!report_unit %in% c("week", "day")) {
    stop("report_unit must be 'week' or 'day'.")
  }
  
  factor_loc <- if (report_unit == "week") 7 else 1
  
  # -- 1) Define a 4-color palette for your methods (in a fixed order)
  #    Even if the 'methods' vector changes names, 
  #    the 1st method in 'methods' gets color 1, 
  #    2nd method gets color 2, etc.
  default_colors <- c("#228B22", "#ffff66", "#8446c6", "#4682B4")
  
  # Subset/assign colors in the exact order that 'methods' appears
  # e.g., if methods = c("A", "B", "C", "D") => method_colors[A] = "#228B22", etc.
  method_colors <- setNames(default_colors[seq_along(methods)], methods)
  
  # -- 2) Labels for the methods (in your preferred textual form)
  #    If you want them always the same, you can define them in a named vector,
  #    then subset by 'methods'. Or simply use the method names themselves.
  #    Example: Named vector with your descriptive labels:
  default_labels <- c("fixed_q"  = "Fixed q",
                      "fixed_b"  = "Fixed b",
                      "linear_b" = "Linear b",
                      "ou_b"     = "OU b")
  
  # Subset these labels in the order of 'methods' 
  # (if user changes the method name to something else, you may need to adapt here).
  method_labels <- default_labels[methods]
  
  # -- 3) Colors for Real Cases and Reported Cases remain the same
  #    We'll combine them with the method colors
  model_colors <- c("Real Cases"     = "red",
                    "Reported Cases" = "black",
                    method_colors)
  
  # Collect the breaks (legend order) and labels in the exact order we want
  legend_breaks <- c("Real Cases", "Reported Cases", methods)
  legend_labels <- c("Real Cases", "Reported Cases", method_labels)
  
  p_out <- list()
  n_runs <- length(nowcasts_list)
  
  for (i in seq_len(n_runs)) {
    nowcasts_df <- nowcasts_list[[i]]
    now <- unique(nowcasts_df$now)
    earliest <- unique(nowcasts_df$earliest)
    last_date_for_delay <- unique(nowcasts_df$last_date_for_delay)
    
    # Build combined data for models
    model_data <- lapply(methods, function(model_name) {
      data.frame(
        date  = nowcasts_df$date,
        mean  = nowcasts_df[[paste0("mean_", model_name)]],
        lower = nowcasts_df[[paste0("lower_", model_name)]],
        upper = nowcasts_df[[paste0("upper_", model_name)]],
        model = model_name
      )
    }) %>% 
      do.call(rbind, .)
    
    # -- Start plot
    p <- ggplot() +
      # 3) Make Real/Reported Cases double width
      geom_line(data = nowcasts_df,
                aes(x = date, y = case_true, color = "Real Cases"),
                linewidth = 2) +
      geom_line(data = nowcasts_df,
                aes(x = date, y = case_reported, color = "Reported Cases"),
                linewidth = 2) +
      
      # Mark 'now'
      annotate("text", x = now, y = -1, 
               label = paste0("now: ", now), 
               hjust = 1, vjust = 2, color = "red") +
      geom_point(data = data.frame(x = now, y = 0), 
                 aes(x = x, y = y, shape = "Today"), 
                 size = 2, color = "red") +
      scale_shape_manual(values = c("Today" = 17), guide = "none") +
      
      # Ribbon for model intervals
      geom_ribbon(data = model_data, 
                  aes(x = date, ymin = lower, ymax = upper, fill = model),
                  alpha = 0.3) +
      geom_line(data = model_data,
                aes(x = date, y = mean, color = model),
                linewidth = 1) +
      
      # Use scale_color_manual with breaks & labels for exact legend order
      scale_color_manual(
        values = model_colors, 
        breaks = legend_breaks,
        labels = legend_labels,
        name   = "Legend"
      ) +
      # Fill for ribbons
      scale_fill_manual(values = method_colors, guide = "none") +
      
      labs(
        title = title,
        x     = x_lab,
        y     = y_lab
      ) +
      theme_minimal() +
      theme(
        legend.position      = legend_position,   
        legend.justification = c(0, 1),
        legend.background    = element_rect(fill = "white", color = "black", 
                                            size = 0.5, linetype = "solid"),
        legend.key           = element_rect(fill = "white", color = NA),
        legend.text          = element_text(size = 12),
        legend.title         = element_text(size = 12),
        axis.text            = element_text(size = 12),
        axis.title           = element_text(size = 12)
      )
    
    # Show vertical line if 'last_date_for_delay' is valid
    if (last_date_for_delay >= earliest) {
      p <- p + 
        geom_vline(xintercept = last_date_for_delay, color = "orange", 
                   linetype = "dashed", size = 1) +
        annotate("text", x = last_date_for_delay, y = -1, 
                 label = last_date_for_delay, vjust = 2, color = "orange")
    }
    
    p_out[[i]] <- p
  }
  
  # Combine plots in grid if requested
  if (combine_plots) {
    if (is.null(ncol) && is.null(nrow)) {
      ncol <- ceiling(sqrt(n_runs))
      nrow <- ceiling(n_runs / ncol)
    }
    
    patchwork_plot <- wrap_plots(p_out, guides = "collect") +
      plot_layout(ncol = ncol, nrow = nrow) +
      plot_annotation(title = title) &
      theme(legend.position = "bottom")
    
    return(patchwork_plot)
    
  } else {
    return(p_out)
  }
}


