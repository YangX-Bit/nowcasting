library(ggplot2)
library(gridExtra)
library(ggforce)
library(lubridate)
library(dplyr)
library(ggpubr)

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
      model_fit <- nls(y ~ (1 - (1 - phi) * exp(-b * x)), data = data_fit, start = list(b = 0.2, phi = 0.9))
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



# nowcasts_plot <- function(nowcasts_list,
#                           D = NULL,
#                           report_unit = "week",
#                           methods = c("fixed_q", "fixed_b", "linear_b", "ou_b"),
#                           title = NULL,
#                           x_lab = NULL,
#                           y_lab = "Cases / Nowcast",
#                           legend_position = NULL,
#                           combine_plots = FALSE,
#                           ncol = NULL,
#                           nrow = NULL) {
#   library(ggplot2)
#   library(lubridate)
#   library(dplyr)
#   library(patchwork)
#   
#   # Basic checks
#   if (is.null(D)) {
#     stop("Parameter 'D' must be provided.")
#   }
#   if (!report_unit %in% c("week", "day")) {
#     stop("report_unit must be 'week' or 'day'.")
#   }
#   
#   factor_loc <- if (report_unit == "week") 7 else 1
#   
#   # -- 1) Define a 4-color palette for your methods (in a fixed order)
#   #    Even if the 'methods' vector changes names, 
#   #    the 1st method in 'methods' gets color 1, 
#   #    2nd method gets color 2, etc.
#   default_colors <- c("#228B22", "#ffff66", "#8446c6", "#4682B4")
#   
#   # Subset/assign colors in the exact order that 'methods' appears
#   # e.g., if methods = c("A", "B", "C", "D") => method_colors[A] = "#228B22", etc.
#   method_colors <- setNames(default_colors[seq_along(methods)], methods)
#   
#   # -- 2) Labels for the methods (in your preferred textual form)
#   #    If you want them always the same, you can define them in a named vector,
#   #    then subset by 'methods'. Or simply use the method names themselves.
#   #    Example: Named vector with your descriptive labels:
#   default_labels <- c("fixed_q"  = "Fixed q",
#                       "fixed_b"  = "Fixed b",
#                       "linear_b" = "Linear b",
#                       "ou_b"     = "OU b")
#   
#   # Subset these labels in the order of 'methods' 
#   # (if user changes the method name to something else, you may need to adapt here).
#   method_labels <- default_labels[methods]
#   
#   # -- 3) Colors for Actual Cases and Reported Cases remain the same
#   #    We'll combine them with the method colors
#   model_colors <- c("Actual Cases"     = "red",
#                     "Reported Cases" = "black",
#                     method_colors)
#   
#   # Collect the breaks (legend order) and labels in the exact order we want
#   legend_breaks <- c("Actual Cases", "Reported Cases", methods)
#   legend_labels <- c("Actual Cases", "Reported Cases", method_labels)
#   
#   p_out <- list()
#   n_runs <- length(nowcasts_list)
#   
#   for (i in seq_len(n_runs)) {
#     nowcasts_df <- nowcasts_list[[i]]
#     now <- unique(nowcasts_df$now)
#     earliest <- unique(nowcasts_df$earliest)
#     last_date_for_delay <- unique(nowcasts_df$last_date_for_delay)
#     
#     # Build combined data for models
#     model_data <- lapply(methods, function(model_name) {
#       data.frame(
#         date  = nowcasts_df$date,
#         mean  = nowcasts_df[[paste0("mean_", model_name)]],
#         lower = nowcasts_df[[paste0("lower_", model_name)]],
#         upper = nowcasts_df[[paste0("upper_", model_name)]],
#         model = model_name
#       )
#     }) %>% 
#       do.call(rbind, .)
#     
#     # -- Start plot
#     p <- ggplot() +
#       # 3) Make Real/Reported Cases double width
#       geom_line(data = nowcasts_df,
#                 aes(x = date, y = case_true, color = "Actual Cases"),
#                 linewidth = 2) +
#       geom_line(data = nowcasts_df,
#                 aes(x = date, y = case_reported, color = "Reported Cases"),
#                 linewidth = 2) +
#       
#       # Mark 'now'
#       annotate("text", x = now, y = -1, 
#                label = paste0("now: ", now), 
#                hjust = 1, vjust = 2, color = "red") +
#       geom_point(data = data.frame(x = now, y = 0), 
#                  aes(x = x, y = y, shape = "Today"), 
#                  size = 2, color = "red") +
#       scale_shape_manual(values = c("Today" = 17), guide = "none") +
#       
#       # Ribbon for model intervals
#       geom_ribbon(data = model_data, 
#                   aes(x = date, ymin = lower, ymax = upper, fill = model),
#                   alpha = 0.3) +
#       geom_line(data = model_data,
#                 aes(x = date, y = mean, color = model),
#                 linewidth = 1) +
#       
#       # Use scale_color_manual with breaks & labels for exact legend order
#       scale_color_manual(
#         values = model_colors, 
#         breaks = legend_breaks,
#         labels = legend_labels,
#         name   = "Legend"
#       ) +
#       # Fill for ribbons
#       scale_fill_manual(values = method_colors, guide = "none") +
#       
#       labs(
#         title = title,
#         x     = x_lab,
#         y     = y_lab
#       ) +
#       theme_minimal() +
#       theme(
#         legend.position      = legend_position,   
#         legend.justification = c(0, 1),
#         legend.background    = element_rect(fill = "white", color = "black", 
#                                             size = 0.5, linetype = "solid"),
#         legend.key           = element_rect(fill = "white", color = NA),
#         legend.text          = element_text(size = 12),
#         legend.title         = element_text(size = 12),
#         axis.text            = element_text(size = 12),
#         axis.title           = element_text(size = 12)
#       )
#     
#     # Show vertical line if 'last_date_for_delay' is valid
#     if (last_date_for_delay >= earliest) {
#       p <- p + 
#         geom_vline(xintercept = last_date_for_delay, color = "orange", 
#                    linetype = "dashed", size = 1) +
#         annotate("text", x = last_date_for_delay, y = -1, 
#                  label = last_date_for_delay, vjust = 2, color = "orange")
#     }
#     
#     p_out[[i]] <- p
#   }
#   
#   # Combine plots in grid if requested
#   if (combine_plots) {
#     if (is.null(ncol) && is.null(nrow)) {
#       ncol <- ceiling(sqrt(n_runs))
#       nrow <- ceiling(n_runs / ncol)
#     }
#     
#     patchwork_plot <- wrap_plots(p_out, guides = "collect") +
#       plot_layout(ncol = ncol, nrow = nrow) +
#       plot_annotation(title = title) &
#       theme(legend.position = "bottom")
#     
#     return(patchwork_plot)
#     
#   } else {
#     return(p_out)
#   }
# }

# nowcasts_plot_separated <- function(nowcasts_list,
#                                     D = NULL,
#                                     report_unit = "week",
#                                     # Preset methods (4 models)
#                                     methods = c("q_constant", "b_constant", "b_rw", "b_ou"),
#                                     title = NULL,
#                                     x_lab = NULL,
#                                     y_lab = NULL,
#                                     lab_comb = "Number of cases",
#                                     combine_plots = TRUE) {
#   library(ggplot2)
#   library(lubridate)
#   library(dplyr)
#   library(patchwork)
#   
#   # Basic checks
#   if (is.null(D)) {
#     stop("Parameter 'D' must be provided.")
#   }
#   if (!report_unit %in% c("week", "day")) {
#     stop("report_unit must be 'week' or 'day'.")
#   }
#   
#   factor_loc <- if (report_unit == "week") 7 else 1
#   
#   # Define colors
#   method_colors <- c(
#     "q_constant" = "#228B22",  # green
#     "b_constant" = "#b35806",  # light yellow
#     "b_rw"       = "#8446c6",  # purple
#     "b_ou"       = "#4682B4"   # steel blue
#   )
#   model_colors <- c("Actual Cases" = "red",
#                     "Reported Cases" = "black",
#                     method_colors)
#   
#   legend_breaks <- c("Actual Cases", "Reported Cases", names(method_colors))
#   legend_labels <- c("Actual Cases", 
#                      "Reported Cases", 
#                      "q-Constant", "b-Constant", "b-Random walk", "b-OU process")
#   
#   # Generate all subplots (no legend)
#   p_out <- list()
#   n_runs <- length(nowcasts_list)
#   
#   idx <- 1
#   for (i in seq_len(n_runs)) {
#     nowcasts_df <- nowcasts_list[[i]]
#     
#     now <- unique(nowcasts_df$now)
#     earliest <- unique(nowcasts_df$earliest)
#     last_date_for_delay <- unique(nowcasts_df$last_date_for_delay)
#     
#     for (model_name in methods) {
#       sub_data <- data.frame(
#         date  = nowcasts_df$date,
#         mean  = nowcasts_df[[paste0("mean_", model_name)]],
#         lower = nowcasts_df[[paste0("lower_", model_name)]],
#         upper = nowcasts_df[[paste0("upper_", model_name)]],
#         model = model_name
#       )
#       
#       p_sub <- ggplot() +
#         geom_line(data = nowcasts_df,
#                   aes(x = date, y = case_true), 
#                   color = "red", linewidth = 1.2) +
#         geom_line(data = nowcasts_df,
#                   aes(x = date, y = case_reported),
#                   color = "black", linewidth = 1.2) +
#         geom_ribbon(data = sub_data,
#                     aes(x = date, ymin = lower, ymax = upper),
#                     fill = method_colors[model_name], alpha = 0.3) +
#         geom_line(data = sub_data,
#                   aes(x = date, y = mean),
#                   color = method_colors[model_name], linewidth = 1) +
#         # annotate("text", x = now, y = -1, 
#         #          label = paste0("now: ", now), 
#         #          hjust = 1, vjust = 2, color = "red",
#         #          size = 6) +
#         geom_point(data = data.frame(x = now, y = 0),
#                    aes(x = x, y = y),
#                    shape = 17, size = 2, color = "red") +
#         {
#           if (last_date_for_delay >= earliest) {
#             list(
#               geom_vline(xintercept = last_date_for_delay, 
#                          color = "black", linetype = "dashed", size = 1),
#               annotate("text", x = last_date_for_delay, y = -1,
#                        label = as.character(last_date_for_delay),
#                        vjust = 2, color = "black")
#             )
#           } else {
#             NULL
#           }
#         } +
#         labs(
#           title = title,
#           x     = x_lab,
#           y     = y_lab
#         ) +
#         theme_minimal() +
#         theme(
#           legend.position = "none",
#           axis.text       = element_text(size = 20),
#           axis.title      = element_text(size = 22),
#           axis.text.x     = element_text(angle = 45, hjust = 1)
#         )
#       
#       p_out[[idx]] <- p_sub
#       idx <- idx + 1
#     }
#   }
#   
#   # -- OPTIONAL STEP: unify y-axis for each row of 4 plots
#   # If each row has exactly 4 subplots, we do the following:
#   # We have n_runs rows, each row = 4 subplots => total = n_runs * 4
#   # For row i in 1..n_runs, we unify p_out[ row_range ]
#   
#   for (row_i in seq_len(n_runs)) {
#     # The subplots for row row_i are indices: (row_i - 1)*4 + 1  to  row_i*4
#     idx_start <- (row_i - 1)*4 + 1
#     idx_end   <- row_i*4
#     # gather data from these 4 subplots
#     sub_list  <- p_out[idx_start:idx_end]
#     
#     # Build each ggplot to get the data
#     builds <- lapply(sub_list, ggplot_build)
#     
#     # Extract all y-values from the data frames (lines, ribbons, etc.)
#     # We combine all layers
#     all_y <- c()
#     for (b in builds) {
#       for (lyr in b$data) {
#         # gather columns that might represent y, ymin, ymax
#         ycols <- intersect(c("y", "ymin", "ymax"), names(lyr))
#         all_y <- c(all_y, unlist(lyr[ycols]))
#       }
#     }
#     # compute global range
#     y_min <- min(all_y, na.rm = TRUE)
#     y_max <- max(all_y, na.rm = TRUE)
#     
#     # Add a small margin if needed, or just use exact
#     # unify for each of the 4 subplots
#     for (k in seq_along(sub_list)) {
#       sub_list[[k]] <- sub_list[[k]] + 
#         scale_y_continuous(limits = c(y_min, y_max))
#     }
#     # put updated plots back to p_out
#     p_out[idx_start:idx_end] <- sub_list
#   }
#   
#   # Combine if requested
#   if (combine_plots) {
#     total_plots <- length(p_out)
#     # for your layout usage:
#     # ncol = 4, nrow = n_runs (since each row has 4 subplots)
#     
#     subplots <- wrap_plots(p_out, ncol = 4, nrow = n_runs)
#     
#     if (!is.null(title)) {
#       subplots <- subplots + plot_annotation(title = title)
#     }
#     
#     # Create separate legend
#     legend_items <- factor(legend_breaks, levels = legend_breaks)
#     df_legend <- data.frame(
#       group = rep(legend_items, each = 2),
#       x     = rep(c(1, 2), times = length(legend_items)),
#       y     = rep(seq(1, length(legend_items)), each = 2)
#     )
#     
#     legend_plot <- ggplot(df_legend, aes(x = x, y = y, color = group)) +
#       geom_line(aes(group = group), size = 1, alpha = 1) + 
#       scale_color_manual(
#         values = model_colors,
#         breaks = legend_breaks,
#         labels = legend_labels,
#         name   = NULL,
#       ) +
#       theme_void() +
#       theme(
#         legend.position   = "right",
#         legend.text       = element_text(size = 20),
#         legend.title      = element_text(size = 22),
#         legend.background = element_rect(fill = "white", color = NA)
#       ) +
#       guides(color = guide_legend(override.aes = list(shape = NA)))
#     
#     
#     # Here you could do your custom row/col layout
#     # For simplicity, we just do side-by-side: subplots | legend_plot
#     final_combined <- (
#       # first row
#       plot_spacer() + p_out[[1]] + p_out[[2]] + p_out[[3]] + p_out[[4]] + plot_spacer() + 
#         
#         # second row, put legend
#         plot_spacer() + p_out[[5]] + p_out[[6]] + p_out[[7]] + p_out[[8]] + legend_plot +
#         
#         #
#         plot_spacer() + p_out[[9]] + p_out[[10]] + p_out[[11]] + p_out[[12]] + plot_spacer() +
#         
#         # 
#         plot_spacer() + p_out[[13]] + p_out[[14]] + p_out[[15]] + p_out[[16]] + plot_spacer() +
#         
#         # 
#         plot_spacer() + p_out[[17]] + p_out[[18]] + p_out[[19]] + p_out[[20]] + plot_spacer()
#     ) + 
#       plot_layout(ncol = 6, widths = c(0.2,1, 1, 1, 1, 0.00001))  # use a small number to hide the last column
#     # first col for title last col for legend
#     
#     final_with_ylab <- ggdraw(final_combined) +
#       draw_label("Number of cases", 
#                  x = 0.02, y = 0.5, angle = 90,  
#                  size = 28, fontface = "bold")
#     
#     return(final_with_ylab)
#     
#   } else {
#     return(p_out)
#   }
# }

nowcasts_plot_separated <- function(results_list, scoreRange_SARI, first_date, report_unit = "week", D) {
  # Load required packages
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  library(rlang)
  
  # Define model names and titles
  models <- c("q_constant", "b_constant", "b_rw", "b_ou")
  model_titles <- c("q-Constant", "b-Constant", "b-Random walk", "b-OU process")
  
  # Define model colors
  method_colors <- c(
    "q_constant" = "#228B22",  # green
    "b_constant" = "#b35806",  # light yellow
    "b_rw"       = "#8446c6",  # purple
    "b_ou"       = "#4682B4"   # steel blue
  )
  
  n_rows <- length(results_list)  # e.g., 5
  n_cols <- length(models)         # e.g., 4
  
  # Calculate offset (in days) for the vertical dashed line:
  offset_days <- if (report_unit == "week") D * 7 else D
  
  # List to store each row's combined plot
  plot_list <- list()
  
  for (i in seq_along(results_list)) {
    df <- results_list[[i]]
    # For current row, X-axis range is from first_date to scoreRange_SARI[i]
    x_min <- first_date
    x_max <- scoreRange_SARI[i]
    
    # Calculate Y-axis upper limit for current row (based on actual and predicted upper values)
    y_max <- max(c(df$case_true, df$case_reported,
                   df$upper_q_constant, df$upper_b_constant, df$upper_b_rw, df$upper_b_ou),
                 na.rm = TRUE)
    
    # Extract current "now" date (assumed constant within df)
    current_now <- unique(df$now)
    # Calculate vertical dashed line position
    proposed_vline <- as.Date(current_now) - offset_days
    # Only add the line if the computed date is not before x_min
    vline_layer <- if (proposed_vline >= x_min) {
      geom_vline(xintercept = proposed_vline, color = "black", linetype = "dashed", size = 1)
    } else {
      NULL
    }
    
    # List to store subplots for each model (columns)
    row_plots <- list()
    for (j in seq_along(models)) {
      mod <- models[j]
      mean_col <- paste0("mean_", mod)
      lower_col <- paste0("lower_", mod)
      upper_col <- paste0("upper_", mod)
      
      # Only the first row's subplots show model titles
      current_title <- if (i == 1) model_titles[j] else NULL
      
      p <- ggplot(df, aes(x = date)) +
        # Plot actual cases (red line)
        geom_line(aes(y = case_true), color = "red", size = 1) +
        # Plot reported cases (black line)
        geom_line(aes(y = case_reported), color = "black", size = 1) +
        # Plot the prediction ribbon for the current model
        geom_ribbon(aes(ymin = !!sym(lower_col), ymax = !!sym(upper_col)),
                    fill = method_colors[mod], alpha = 0.3) +
        # Plot the model prediction (mean)
        geom_line(aes(y = !!sym(mean_col)), color = method_colors[mod], size = 1) +
        # Conditionally add the vertical dashed line
        { vline_layer } +
        scale_x_date(limits = c(x_min, x_max)) +
        scale_y_continuous(limits = c(0, y_max)) +
        labs(title = current_title, x = NULL, y = NULL) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text = element_text(size = 8)
        )
      
      row_plots[[j]] <- p
    }
    
    # Combine the 4 subplots horizontally for the current row
    row_combined <- wrap_plots(row_plots, ncol = n_cols)
    
    # Create a right-side label plot to display the current date (scoreRange_SARI[i])
    label_plot <- ggdraw() +
      draw_label(label = as.character(scoreRange_SARI[i]),
                 fontface = "bold", size = 12, angle = 270)
    # Combine row plots with the right-side label (adjust relative widths as needed)
    row_with_label <- plot_grid(row_combined, label_plot, ncol = 2, rel_widths = c(1, 0.01))
    
    plot_list[[i]] <- row_with_label
  }
  
  # Combine all rows vertically
  combined_plot <- wrap_plots(plot_list, ncol = 1)
  
  dummy_data <- data.frame(
    x = rep(1:2, 6),
    y = rep(c(1,2), 6),
    group = factor(rep(legend_breaks, each = 2), levels = legend_breaks)
  )
  
  dummy_plot <- ggplot(dummy_data, aes(x = x, y = y, color = group)) +
    geom_line(size = 1) +
    scale_color_manual(
      name = "Models",
      breaks = legend_breaks,
      labels = legend_labels,
      values = model_colors
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  g <- ggplotGrob(dummy_plot) # get legend
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  
  # Add a left-side Y axis label with additional margin.
  final_plot <- ggdraw(combined_plot) +
    draw_label("Number of cases", x = -0.01, y = 0.5, angle = 90, size = 14, fontface = "bold")
  
  # Wrap the final plot with additional margin using plot_annotation
  final_with_annotation <- final_plot + 
    plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 20)))
  
  # Combine the final plot and legend (legend at bottom)
  final_with_legend <- plot_grid(final_with_annotation, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  return(final_with_legend)
}


nowcasts_plot_single <- function(nowcasts_list,
                                 D = NULL,
                                 report_unit = "week",
                                 # Single method to be plotted
                                 method = "b_ou",
                                 title = NULL,
                                 x_lab = NULL,
                                 y_lab = NULL,
                                 lab_comb = "Number of cases",
                                 combine_plots = TRUE) {
  library(ggplot2)
  library(lubridate)
  library(dplyr)
  library(patchwork)
  library(cowplot)
  
  # Basic checks
  if (is.null(D)) {
    stop("Parameter 'D' must be provided.")
  }
  if (!report_unit %in% c("week", "day")) {
    stop("report_unit must be 'week' or 'day'.")
  }
  
  factor_loc <- if (report_unit == "week") 7 else 1
  
  # Define colors
  method_colors <- c(
    "q_constant" = "#228B22",  # green
    "b_constant" = "#b35806",  # light yellow
    "b_rw"       = "#8446c6",  # purple
    "b_ou"       = "#4682B4"   # steel blue
  )
  
  # Define legend mappings for model names
  legend_breaks <- c("Final cumulative cases", "Reported Cases", names(method_colors))
  legend_labels <- c("Final cumulative cases", 
                     "Reported Cases", 
                     "q-Constant", "b-Constant", "b-Random walk", "b-OU process")
  
  # Generate subplots for a single method
  p_out <- list()
  n_runs <- length(nowcasts_list)
  idx <- 1
  
  for (i in seq_len(n_runs)) {
    nowcasts_df <- nowcasts_list[[i]]
    
    now <- unique(nowcasts_df$now)
    earliest <- unique(nowcasts_df$earliest)
    last_date_for_delay <- unique(nowcasts_df$last_date_for_delay)
    
    sub_data <- data.frame(
      date  = nowcasts_df$date,
      mean  = nowcasts_df[[paste0("mean_", method)]],
      lower = nowcasts_df[[paste0("lower_", method)]],
      upper = nowcasts_df[[paste0("upper_", method)]],
      model = method
    )
    
    # Determine whether to show legend
    show_legend <- (i == n_runs && combine_plots)
    
    p_sub <- ggplot() +
      geom_line(data = nowcasts_df,
                aes(x = date, y = case_reported, linetype = "Reported cases"),
                color = "black", linewidth = 1) +
      geom_line(data = nowcasts_df,
                aes(x = date, y = case_true, linetype = "Final cumulative cases"), 
                color = "red", linewidth = 1) +
      geom_ribbon(data = sub_data,
                  aes(x = date, ymin = lower, ymax = upper, fill = model),
                  alpha = 0.3) +
      geom_line(data = sub_data,
                aes(x = date, y = mean, color = model),
                linewidth = 1) +
      # Removed the annotate("text", ...) that displayed "now: ..."
      # geom_point(data = data.frame(x = now, y = 0),
      #            aes(x = x, y = y),
      #            shape = 17, size = 1.5, color = "red") +
      {
        if (last_date_for_delay >= earliest) {
          list(
            geom_vline(xintercept = last_date_for_delay, 
                       color = "black", linetype = "twodash", size = 1),
            annotate("text", x = last_date_for_delay, y = -0.1,
                     label = as.character(last_date_for_delay),
                     vjust = 2, color = "black", size = 8)
          )
        } else {
          NULL
        }
      } +
      scale_color_manual(values = method_colors,
                         #name = "Model",
                         name = NULL,
                         breaks = names(method_colors),
                         labels = c("q-Constant", "b-Constant", "b-Random walk", "b-OU process")) +
      scale_fill_manual(values = method_colors, 
                        name = "Prediction interval",
                        breaks = names(method_colors),
                        labels = c("q-Constant", "b-Constant", "b-Random walk", "b-OU process")) +
      # Use proper line types that are supported in R
      scale_linetype_manual(values = c("Final cumulative cases" = "dotted", "Reported cases" = "dashed"), 
                            name = "Case type") +
      labs(
        title = paste0(now),  # Updated title to include the current date
        x     = x_lab,
        y     = y_lab
      ) +
      theme_minimal() +
      theme(
        legend.position = if(show_legend) "right" else "none",
        legend.text = element_text(size = 14, family = "serif"),
        legend.title = element_text(size = 14, family = "serif", face = "bold"),
        axis.text = element_text(size = 14, family = "serif"),
        axis.title = element_text(size = 14, family = "serif", face = "bold"),
        axis.text.x = element_text( hjust = 1, size = 14, family = "serif"),
        plot.title = element_text(size = 14, family = "serif", face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 15)
      )
    
    p_out[[idx]] <- p_sub
    idx <- idx + 1
  }
  
  # Unify y-axis across subplots
  for (row_i in seq_len(n_runs)) {
    sub_plot <- p_out[[row_i]]
    build <- ggplot_build(sub_plot)
    all_y <- c()
    for (lyr in build$data) {
      ycols <- intersect(c("y", "ymin", "ymax"), names(lyr))
      all_y <- c(all_y, unlist(lyr[ycols]))
    }
    y_min <- min(all_y, na.rm = TRUE)
    y_max <- max(all_y, na.rm = TRUE)
    
    p_out[[row_i]] <- sub_plot + scale_y_continuous(limits = c(y_min, y_max))
  }
  
  # Combine subplots if requested (one column per period)
  if (combine_plots) {
    # Set all but the last subplot to have no legend
    for (i in 1:(n_runs-1)) {
      p_out[[i]] <- p_out[[i]] + theme(legend.position = "none")
    }
    
    # Ensure the last subplot has a legend on the right
    p_out[[n_runs]] <- p_out[[n_runs]] + 
      theme(legend.position = "right",
            legend.box.margin = margin(0, 0, 0, 20))
    
    combined_plots <- wrap_plots(p_out, ncol = 3, nrow = 2)
    
    if (!is.null(title)) {
      combined_plots <- combined_plots + plot_annotation(
        title = title,
        plot.margin = margin(10, 10, 300, 50),
        theme = theme(plot.title = element_text(size = 16, family = "serif", face = "bold", hjust = 0.5))
      )
    }
    
    # Adjust y-axis label position to avoid overlap with the first column of plots
    final_with_ylab <- ggdraw(combined_plots + plot_layout(guides = "collect")) +
      draw_label(lab_comb, 
                 x = -0.01, y = 0.5, angle = 90,  
                 size = 14, fontface = "bold") +
      theme(plot.margin = margin(10, 10, 10, 30))
    
    return(final_with_ylab)
  }
  else {
    return(p_out)
  }
}

#' Plot multiple curves with different model types
#'
#' @param model A CmdStanModel fit result, used to extract draws (samples).
#' @param x_range A numeric vector for the x values you want to plot over (e.g., 1:N).
#' @param num_rows Number of rows in facet plot.
#' @param num_cols Number of columns in facet plot.
#' @param exp_grow_func A function that generates the curve (e.g., exponential growth) 
#'        given x and parameters. You can define your own logic.
#' @param raw_data A data frame or list containing “black line” parameters 
#'        (if you have them). This could store b_t, phi_t, p, or other variables 
#'        you want to show in black lines for comparison.
#' @param model_type Character string indicating which type of model we are dealing with:
#'        "fixed_q", "fixed_b", "rw_b", or "ou_b".
#' @param D (Optional) The dimension/length of p if model_type = "fixed_q" (e.g., D+1).
#' @param N_obs (Optional) The length of b_t and phi_t if model_type = "rw_b"/"ou_b".
#'
#' This function extracts parameters based on the model type, then plots both black lines 
#' (from some external source, e.g. raw_data) and red lines (from model draws) with 
#' uncertainty intervals (95% CI).
#'

plot_fitted_qd <- function(model, 
                           D = 15,
                           N_obs = 60,
                           x_range = c(0:D), 
                           num_rows = 5, 
                           num_cols = 5, 
                           raw_data = NULL, 
                           model_type = c("fixed_q", "fixed_b", "rw_b", "ou_b"),
                           alpha = 0.05
) {
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  model_type <- match.arg(model_type)
  
  exp_grow_func <- function(d, b, phi){
    1 - (1 - phi)*exp(-b * d)
  }
  
  upr <- 1 - alpha/2; lwr <- alpha/2
  #------------------------------------------------------------------
  # 1) Extract parameters from CmdStan model draws
  #    The extraction strategy depends on the model type
  #------------------------------------------------------------------
  
  if (model_type == "fixed_q") {
    # In a "fixed_q" model, we extract 'p' which is of length D+1
    # Here we assume the parameter name in the Stan model is "p".
    # You can adjust the code if your parameter is named differently.
    
    p_draws <- model$draws(variables = c("p"), format = "draws_df")
    # Suppose we only want the first 30 draws to show in plots (just as an example):
    p_samples <- p_draws[, 1:(D+1)]
    
    # Compute mean and 95% CI across the draws
    mean_values_p <- apply(p_samples, 2, mean)
    ci_lower_p <- apply(p_samples, 2, function(x) quantile(x, lwr))
    ci_upper_p <- apply(p_samples, 2, function(x) quantile(x, upr))
    
    # You can then use p to generate some curve if needed.
    # For demonstration, I'll define: y_red = exp_grow_func(x_range, p).
    # But you might have your own logic to handle p.
    # For example, if p is a vector controlling daily shape, you might have 
    # a function like daily_func(x, p_vector). 
    # Modify as needed:
    
    # Example placeholder
    # We'll treat mean_values_p[i] as a single parameter that goes into exp_grow_func,
    # but in your real code, you'd likely define a custom function that uses p as a vector.
    
  } else if (model_type == "fixed_b") {
    # In a "fixed_b" model, we extract fixed b and phi (both scalars, presumably)
    b_draws   <- model$draws(variables = "b",   format = "draws_df")
    phi_draws <- model$draws(variables = "phi", format = "draws_df")
    
    b_samples   <- b_draws$b
    phi_samples <- phi_draws$phi
    
    mean_values_b   <- mean(b_samples)
    ci_lower_b      <- quantile(b_samples, lwr)  
    ci_upper_b      <- quantile(b_samples, upr)  
    
    mean_values_phi   <- mean(phi_samples)
    ci_lower_phi      <- quantile(phi_samples, lwr)  
    ci_upper_phi      <- quantile(phi_samples, upr)  
    
  } else if (model_type %in% c("rw_b", "ou_b")) {
    # In a "rw_b" (random walk) model, we assume b_t and phi_t are length = N_obs
    b_t_draws   <- model$draws(variables = "b_t",   format = "draws_df")
    phi_t_draws <- model$draws(variables = "phi_t", format = "draws_df")
    
    # Suppose each column in b_t_draws or phi_t_draws is a sample dimension 
    # (b_t[1], b_t[2], ..., b_t[N_obs]) repeated across MCMC draws.
    b_t_samples   <- b_t_draws[, 1:N_obs]
    phi_t_samples <- phi_t_draws[, 1:N_obs]
    
    mean_values_b_t   <- apply(b_t_samples, 2, mean)
    ci_lower_b_t      <- apply(b_t_samples, 2, function(x) quantile(x, lwr))
    ci_upper_b_t      <- apply(b_t_samples, 2, function(x) quantile(x, upr))
    
    mean_values_phi_t <- apply(phi_t_samples, 2, mean)
    ci_lower_phi_t    <- apply(phi_t_samples, 2, function(x) quantile(x, lwr))
    ci_upper_phi_t    <- apply(phi_t_samples, 2, function(x) quantile(x, upr))
    
  }
  #------------------------------------------------------------------
  # 2) Construct data for plotting
  #    We will demonstrate a generic approach:
  #    - "Black line" data might come from raw_data (your external data).
  #    - "Red line" data and ribbons (95% CI) come from the means/CI of the model draws.
  #------------------------------------------------------------------
  
  # Let's assume raw_data has the same structure as in your original code example: 
  # b_t[c(1:30)], phi[c(1:30)], etc.
  # In practice, you might store it differently depending on the model_type.
  #
  # We'll create a list of data frames, one per facet (num_rows * num_cols).
  # Each data frame has x, y, group, lower, upper, plot_id. Then we'll bind_rows them.
  
  # Here is a generic placeholder approach that tries to mimic your original code 
  # with black lines from 'raw_data' and red lines from the mean/CI. 
  # You will likely adjust to your actual needs (particularly for "fixed_q").
  
  # Number of subplots we want to create:
  n_plots <- num_rows * num_cols
  
  plot_data_list <- lapply(seq_len(n_plots), function(i) {
    #---------------------------------------
    # BLACK CURVE:
    #---------------------------------------
    # For demonstration, let's pretend that raw_data$b_t and raw_data$phi 
    # each store at least 30 values that we can index with [i].
    # (In real usage, you might do something else or might not have black lines for all models.)
    
    if (model_type != "fixed_q") {
      # Example:
      b_for_black   <- raw_data$b[i]
      phi_for_black <- raw_data$phi[i]
      
      # Evaluate black curve
      y_black <- exp_grow_func(x_range, b_for_black, phi_for_black)
    } else {
      y_black <- raw_data$qd
    }
    
    #---------------------------------------
    # RED CURVE & UNCERTAINTY:
    #---------------------------------------
    # Here we do a simple demonstration for "fixed_b", "rw_b", or "ou_b" 
    # using b/phi or b_t[i], phi_t[i] to produce a red curve. 
    # For "fixed_q", you'd define your own logic to incorporate p.
    #---------------------------------------
    
    if (model_type == "fixed_q") {
      # If p is length D+1, you'd likely have a function that uses all of those p 
      # values at once. Below is a simple placeholder:
      # mean_values_p[i] is a single number, but p is actually a vector of length D+1.
      # So you'd use the entire vector for your real function. Here, we just demo:
      p_mean  <- mean_values_p
      p_lower <- ci_lower_p
      p_upper <- ci_upper_p
      
      y_red        <- cumsum(p_mean)  # placeholder usage
      y_red_lower  <- cumsum(p_lower) # placeholder usage
      y_red_upper  <- pmin(cumsum(p_upper), 1)# placeholder usage
      
    } else if (model_type == "fixed_b") {
      
      y_red       <- exp_grow_func(x_range, mean_values_b,   mean_values_phi)
      y_red_lower <- exp_grow_func(x_range, ci_lower_b,  ci_lower_phi)
      y_red_upper <- pmin(exp_grow_func(x_range, ci_upper_b,  ci_upper_phi), 1)
      
    } else if (model_type %in% c("rw_b", "ou_b")) {
      # b_t and phi_t are length N_obs each. We are just indexing i for demonstration.
      # In practice, you might want to visualize all b_t or phi_t across time, 
      # or multiple subplots, etc.
      if (i <= length(mean_values_b_t)) {
        b_mean_t   <- mean_values_b_t[i]
        b_lower_t  <- ci_lower_b_t[i]
        b_upper_t  <- ci_upper_b_t[i]
        
        phi_mean_t <- mean_values_phi_t[i]
        phi_lower_t<- ci_lower_phi_t[i]
        phi_upper_t<- ci_upper_phi_t[i]
        
        y_red       <- exp_grow_func(x_range, b_mean_t,   phi_mean_t)
        y_red_lower <- exp_grow_func(x_range, b_lower_t,  phi_lower_t)
        y_red_upper <- exp_grow_func(x_range, b_upper_t,  phi_upper_t)
      } else {
        y_red       <- rep(NA, length(x_range))
        y_red_lower <- rep(NA, length(x_range))
        y_red_upper <- rep(NA, length(x_range))
      }
      
    } else {
      # Fallback
      y_red       <- rep(NA, length(x_range))
      y_red_lower <- rep(NA, length(x_range))
      y_red_upper <- rep(NA, length(x_range))
    }
    
    #---------------------------------------
    # Assemble data frame for subplot i
    #---------------------------------------
    data.frame(
      x = rep(x_range, 2),  # repeat x for black + red
      y = c(y_black, y_red),
      group = c(rep("Black", length(x_range)), rep("Red", length(x_range))),
      lower = c(rep(NA, length(x_range)), y_red_lower),  # only valid for red
      upper = c(rep(NA, length(x_range)), y_red_upper),  # only valid for red
      plot_id = i  # used to distinguish subplots
    )
  })
  
  # Combine list of data frames
  plot_data <- bind_rows(plot_data_list)
  
  #------------------------------------------------------------------
  # 3) Plot using ggplot2
  #------------------------------------------------------------------
  ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_line(aes(size = group), show.legend = FALSE) +  
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), 
                alpha = 0.2, color = NA) +
    scale_color_manual(values = c("Black" = "black", "Red" = "red")) +
    scale_fill_manual(values = c("Black" = NA, "Red" = "red")) +
    scale_size_manual(values = c("Black" = 1.5, "Red" = 1)) +  # width for lines
    facet_wrap(~plot_id, nrow = num_rows, ncol = num_cols) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = paste("Comparison of Curves with Uncertainty Intervals -", model_type),
         x = "X", y = "Y")
}



