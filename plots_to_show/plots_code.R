path_proj = here::here()
path_source = file.path(path_proj, "source")


######### exp curve ######### 

library(ggplot2)
library(latex2exp)

# Define function for q_{t,d}
q_func <- function(d, phi, b) {
  return(1 - phi * exp(-b * d))
}

# Define parameter values
phi_values <- c(0.9, 0.5)    # Two phi values
b_values <- c(0.2, 0.5, 1)   # Three b values
d_values <- seq(0, 10, length.out = 100)  # Define d range

# Create data frame for plotting
plot_data <- expand.grid(d = d_values, phi = phi_values, b = b_values)
plot_data$q_td <- mapply(q_func, plot_data$d, plot_data$phi, plot_data$b)

# Convert phi and b to factor for better legend formatting
plot_data$phi <- factor(plot_data$phi, labels = c(expression(phi == 0.9), expression(phi == 0.5)))
plot_data$b <- factor(plot_data$b, labels = c("b = 0.2", "b = 0.5", "b = 1"))

# Generate plot
exp_curve <- ggplot(plot_data, aes(x = d, y = q_td, color = b, linetype = phi)) +
  geom_line(size = 1) + 
  labs(
    title = TeX("$q_{t,d} = (1 - \\phi \\cdot e^{-b \\cdot d})$"),
    x = "Delay (d)", 
    y = TeX("$q_{t,d}$"),
    color = "b values",     # Legend for colors (b)
    linetype = "phi values" # Legend for line types (phi)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(filename = file.path(path_proj, "plots_to_show", "exp_curve.png"),
       plot = exp_curve,
       width = 14, height = 5, dpi = 300)

########################### 


############### qtd for each scenarios ##################
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load your simulation functions
source(file.path(path_source, "simulation", "simulations_functions_final.R"))

# ---- SETUP PARAMETERS ----
alpha_increase_seq_1 <- seq(10, 750, by = 30)
alpha_decrease_seq_1 <- seq(750, 10, by = -30)
alpha_lamb <- c(rep(10,5), alpha_increase_seq_1 + rnorm(length(alpha_increase_seq_1), 10, 10),
                alpha_decrease_seq_1 + rnorm(length(alpha_decrease_seq_1), 10, 10),
                rep(10,5))
beta_lamb <- 0.5
T <- 60
D <- 15

# ---- SCENARIOS ----
scenarios <- list(
  "FR Constant" = list(method = "q_constant", params = list(b = 0.7, phi = 0.9)),
  "FR RW" = list(method = "b_rw", params = list(mu_logb = log(0.7), sigma_logb = 0.1, mu_logitphi = 1, sigma_logitphi = 0.1)),
  "FR OU" = list(method = "b_ou", params = list(theta_logb = 0.3, mu_logb = log(0.7), init_logb = log(0.7), sigma_logb = 0.2,
                                                theta_logitphi = 0.3, mu_logitphi = 1, init_logitphi = 1, sigma_logitphi = 0.2)),
  "NFR Constant" = list(method = "q_constant", params = list(b = 0.7, phi = 0.9)),
  "NFR RW" = list(method = "b_rw", params = list(mu_logb = log(0.2), sigma_logb = 0.1, mu_logitphi = 1.5, sigma_logitphi = 0.1)),
  "NFR OU" = list(method = "b_ou", params = list(theta_logb = 0.2, mu_logb = log(0.2), init_logb = log(0.2), sigma_logb = 0.15,
                                                 theta_logitphi = 0.2, mu_logitphi = 1.5, init_logitphi = 1.5, sigma_logitphi = 0.15))
)

# ---- RUN SIMULATIONS ----
plot_data <- data.frame()

for (scenario_name in names(scenarios)) {
  params <- scenarios[[scenario_name]]
  
  sim_params <- list(
    data = list(
      alpha_lamb = alpha_lamb,
      beta_lamb  = beta_lamb,
      T          = T,
      date_start = as.Date("2024-01-01"),
      D          = D
    ),
    q_model = list(
      method        = params$method,
      method_params = params$params
    )
  )
  
  sim_result <- simulateData(sim_params)
  
  if (params$method == "q_constant") {
    # Constant Model: Single Line
    temp_q <- data.frame(d = 1:D, q_td = sim_result$q[1:D])
    temp_q$scenario <- scenario_name
    temp_q$replicate <- 1  # Only one realization
    
  } else {
    # RW & OU Models: Multiple Lines (matplot style)
    q_matrix <- as.data.frame(sim_result$q)  # Convert matrix to data frame
    q_matrix$t <- 1:T  # Add time column
    q_long <- pivot_longer(q_matrix, cols = -t, names_to = "d", values_to = "q_td")
    q_long$scenario <- scenario_name
    
    # Convert d column to numeric
    q_long$d <- as.numeric(sub("V", "", q_long$d))
    
    temp_q <- q_long
  }
  
  plot_data <- bind_rows(plot_data, temp_q)
}

# Define correct row-wise scenario order
scenario_order <- c(
  "FR Constant", "FR RW", "FR OU",  
  "NFR Constant", "NFR RW", "NFR OU"
)

# Apply correct factor order to both facets and legend
plot_data$scenario <- factor(plot_data$scenario, levels = scenario_order)

# ---- PLOT q_td ----
qtd <- ggplot(plot_data, aes(x = d, y = q_td, group = t, color = scenario)) +
  geom_line(alpha = 0.7, size = 0.8) +
  facet_wrap(~scenario, nrow = 2) +  # Arrange facets row-wise
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "brown")) +  
  labs(
    title = "Evolution of q over Delay (D)",
    x = "Delay (d)", 
    y = "q_{t,d}",
    color = "Scenario"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +  # Force legend to match row-wise order
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

ggsave(filename = file.path(path_proj, "plots_to_show", "sims_scenarios.png"),
       plot = qtd,
       width = 12, height = 8, dpi = 300)

###################################################### 
# b shape#
########### RW ##########

params_b_rw_FR <- list(
  data = list(
    alpha_lamb = alpha_lamb,
    beta_lamb  = beta_lamb,
    T       = T,
    date_start = as.Date("2024-01-01"),
    D = D
  ),
  q_model = list(
    method        = "b_rw",
    method_params = list(mu_logb = log(0.7), sigma_logb = 0.1, mu_logitphi = 1, sigma_logitphi = 0.1)
  )
)

b_rw_FR <- simulateData(params_b_rw_FR)


page_num <- ceiling(nrow(b_rw_FR$case_reported_cumulated)/16)
exp_plot_b_rw <- fit_exp_plot(b_rw_FR$case_reported_cumulated,ncol = 4, nrow = 4, page = c(1:page_num), if_fit = T)
exp_b_data_b_rw<- data.frame( date = as.Date(rownames(b_rw_FR$case_reported_cumulated)),
                              b = exp_plot_b_rw$coefficients$b)
exp_b_plot_rw <- ggplot(exp_b_data_b_rw, aes(x = date, y = b)) +
  geom_point(color = "black", size = 1.5) +
  geom_smooth(method = "loess", se = TRUE,
              color = "blue", fill = "grey", alpha = 0.5) +
  theme_minimal() +
  labs(x = NULL, y = "Y", title = NULL)

ggsave(filename = file.path(path_proj, "plots_to_show", "rw_b_shape_FR.png"),
       plot = exp_b_plot_rw,
       width = 16, height = 6, dpi = 300)

########### OU ##########
params_b_ou_FR <- list(
  data = list(
    alpha_lamb = alpha_lamb,
    beta_lamb  = beta_lamb,
    T       = T,
    date_start = as.Date("2024-01-01"),
    D = D
  ),
  q_model = list(
    method        = "b_ou",
    method_params = list(theta_logb = 0.3, mu_logb = log(0.7), init_logb = log(0.7), sigma_logb = 0.2,
                         theta_logitphi = 0.3, mu_logitphi = 1, init_logitphi = 1, sigma_logitphi = 0.2)
  )
)
b_ou_FR <- simulateData(params_b_ou_FR)

page_num <- ceiling(nrow(b_ou_FR$case_reported_cumulated)/16)
exp_plot_ou <- fit_exp_plot(b_ou_FR$case_reported_cumulated,ncol = 4, nrow = 4, page = c(1:page_num), if_fit = T)
exp_b_data_ou<- data.frame( date = as.Date(rownames(b_ou_FR$case_reported_cumulated)),
                            b = exp_plot_ou$coefficients$b)
exp_b_plot_ou <- ggplot(exp_b_data_ou, aes(x = date, y = b)) +
  geom_point(color = "black", size = 1.5) +
  geom_smooth(method = "loess", se = TRUE,
              color = "blue", fill = "grey", alpha = 0.5) +
  theme_minimal() +
  labs(x = NULL, y = "Y", title = NULL)

ggsave(filename = file.path(path_proj, "plots_to_show", "ou_b_shape_FR.png"),
       plot = exp_b_plot_ou,
       width = 16, height = 6, dpi = 300)
