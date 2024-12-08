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

nowcasts_plot <- function(results_list, models_to_run = c("fixed_q", "fixed_b", "b_poly", "b_spline"),
                          title = NULL, x_lab = NULL, y_lab = "Cases / Nowcast") {
  p_out <- list(); nowcasts_out <- list();
  for (i in 1:length(results_list$case_true)) {
    # Extract data for the specific time window
    case_true <- results_list[["case_true"]][[i]]
    case_reported <- results_list[["case_reported"]][[i]]
    dates <- results_list[["dates"]][[i]]
    
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
      "fixed_b" = "#FF8C00",  # 
      "b_poly" = "#DAA520",  # 
      "b_spline" = "#4682B4" # 
    )
    
    # Create the base ggplot object
    p <- ggplot(nowcasts, aes(x = date)) +
      geom_line(aes(y = case_true, color = "Real Cases")) +
      geom_line(aes(y = case_reported, color = "Reported Cases"))
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
    
    p <- p +
      scale_color_manual(values = model_colors) +
      labs(title = title,
           x = x_lab,
           y = y_lab,
           color = NULL) +
      theme_minimal() +
      theme(
        legend.position = c(0.1, 0.9),  # Legend position
        legend.justification = c(0, 1),  # Top-left alignment
        legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"), # Legend border
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)
      )
    p_out[[i]] <- p
  }
  return(list(plots = p_out,
              nowcasts = nowcasts_out))
}

