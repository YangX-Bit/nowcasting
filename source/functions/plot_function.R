nowcasts_plot <- function(sample_time_varying, sample_fixed_q,
                          N_obs,
                          title = NULL, x_lab = NULL, y_lab = "Cases / Nowcast"){
  
  nowcasts <- data.frame(mean = apply(sample_time_varying, 2, mean),
                         lower = apply(sample_time_varying, 2, quantile, probs = 0.025),
                         upper = apply(sample_time_varying, 2, quantile, probs = 0.975),
                         #
                         mean_fixped_q = apply(sample_fixed_q, 2, mean),
                         lower_fixped_q = apply(sample_fixed_q, 2, quantile, probs = 0.025),
                         upper_fixped_q = apply(sample_fixed_q, 2, quantile, probs = 0.975),
                         #
                         date =  data$date[1:N_obs],
                         case_true = data$case_true[1:N_obs],
                         case_reported = apply(data_trunc[1:N_obs,], 1, max))
  
  p <- ggplot(nowcasts, aes(x = date)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.5) +
    geom_line(aes(y = mean, , color = "Nowcasts(time varying q)")) +
    geom_ribbon(aes(ymin = lower_fixped_q, ymax = upper_fixped_q), fill = "green", alpha = 0.5) +
    geom_line(aes(y = mean_fixped_q, , color = "Nowcasts(fixed q)")) +
    geom_vline(xintercept = as.Date(now - D), color = "orange") +
    annotate("text", x = as.Date(now - D), y = 0, label = format(as.Date(now - D), "%b %d"),
             vjust = 1, hjust = -0.1, color = "orange") +
    geom_line(aes(y = case_true, color = "Real Cases")) +
    geom_line(aes(y = case_reported, color = "Reported Cases")) +
    scale_color_manual(values = c("Real Cases" = "red", "Reported Cases" = "black",
                                  "Nowcasts(time varying q)" = "blue",
                                  "Nowcasts(fixed q)" = "green")) +
    labs(title = title,
         x = x_lab,
         y = y_lab,
         color = NULL) +
    theme_minimal() +
    theme(
      legend.position = c(0.1, 0.9),  # 
      legend.justification = c(0, 1),  # on left-up
      legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"), # border
      legend.key = element_rect(fill = "white", color = NA),
      
      legend.text = element_text(size = 16),       
      legend.title = element_text(size = 16),       
      axis.text = element_text(size = 16),         
      axis.title = element_text(size = 16)      
    )
  return(p)
}
