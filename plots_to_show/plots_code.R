path_proj = here::here()
path_source = file.path(path_proj, "source")


######### exp curve ######### 
library(ggplot2)
library(latex2exp)
library(tidyr)

# Define function for q_{t,d}
q_exp <- function(d, phi, b) { return(1 - (1 - phi) * exp(-b * d)) }
q_gom <- function(d, phi, b) { return(exp(log(phi) * exp(-b * d))) }

# Define parameter values
phi_values <- c(0.01, 0.25)    # Two phi values: 0.5 and 0.9
b_values <- c(0.3, 0.6, 1.2)     # Three b values
d_values <- seq(0, 10, length.out = 100)  # Define d range

# Create data frame for plotting
plot_data <- expand.grid(d = d_values, phi = phi_values, b = b_values)
plot_data$Exponential <- mapply(q_exp, plot_data$d, plot_data$phi, plot_data$b)
plot_data$Gompertz <- mapply(q_gom, plot_data$d, plot_data$phi, plot_data$b)
plot_data <- pivot_longer(plot_data, Exponential:Gompertz, names_to = "model", values_to = "q_td")

# Convert phi and b to factor with appropriate labels.
# For phi, use strings that can be parsed into plotmath expressions.
plot_data$phi <- factor(plot_data$phi, levels = phi_values, 
                        labels = c("phi == 0.01", "phi == 0.25"))
plot_data$b <- factor(plot_data$b, levels = b_values, 
                      labels = c("b = 0.3", "b = 0.6", "b = 1.2"))
plot_data$group <- paste0(plot_data$model, ":~~", plot_data$phi)


# Generate plot with facets for phi values
exp_curve <- ggplot(plot_data, aes(x = d, y = q_td, linetype = factor(b))) +
  geom_line(size = 1, linewidth = rel(0.4), color = "red") +
  labs(title = NULL, x = TeX("$d$"), y = TeX("$q(d)$"), color = NULL, linetype = NULL) +
  facet_wrap(~ group, labeller = label_parsed, scales = "free", nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("darkblue", "red", "black")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme_classic(9) +
  theme(
    legend.position = c(0.95, 0.3),
    legend.key.height = unit(0.3, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
exp_curve

ggsave(filename = file.path(path_proj, "plots_to_show", "qd_models.png"),
       plot = exp_curve,
       width = 8, height = 2.2, dpi = 300)

########################### 


############### qtd for each scenarios ##################
set.seed(1)
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
  "S1: FR - Constant" = list(method = "q_constant", params = list(b = 0.7, phi = 0.9)),
  "S2: FR - RW" = list(method = "b_rw", params = list(mu_logb = log(0.7), sigma_logb = 0.1, mu_logitphi = 1, sigma_logitphi = 0.1)),
  "S3: FR - OU" = list(method = "b_ou", params = list(theta_logb = 0.3, mu_logb = log(0.7), init_logb = log(0.7), sigma_logb = 0.2,
                                                theta_logitphi = 0.3, mu_logitphi = 1, init_logitphi = 1, sigma_logitphi = 0.2)),
  "S4: NFR - Constant" = list(method = "q_constant", params = list(b = 0.2, phi = 0.9)),
  "S5: NFR - RW" = list(method = "b_rw", params = list(mu_logb = log(0.2), sigma_logb = 0.1, mu_logitphi = 1.5, sigma_logitphi = 0.1)),
  "S6: NFR - OU" = list(method = "b_ou", params = list(theta_logb = 0.2, mu_logb = log(0.2), init_logb = log(0.2), sigma_logb = 0.15,
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
  "S1: FR - Constant", "S2: FR - RW", "S3: FR - OU",  
  "S4: NFR - Constant", "S5: NFR - RW", "S6: NFR - OU"
)

# Apply correct factor order to both facets and legend
plot_data$scenario <- factor(plot_data$scenario, levels = scenario_order)

plot_data <- separate_wider_regex(plot_data, scenario, cols_remove = FALSE,
    c(label = ".*", ": ", fr = ".*", " - ", model = ".*"))
plot_data$fr <- factor(plot_data$fr, levels = c("FR", "NFR"),
    labels = c("Fully reported (FR)", "Non-fully reported (NFR)"))
plot_data$model <- factor(plot_data$model, levels = c("Constant", "RW", "OU"),
    labels = c("Constant", "Random walks (RW)", "Ornstein-Uhlenbeck (OU) processes"))

df_labs <- data.frame(label = paste0("(S", 1:6, ")"), x = rep(0.1, 6), y = rep(0.99, 6),
    fr = rep(c("Fully reported (FR)", "Non-fully reported (NFR)"), each = 3),
    model = rep(c("Constant", "Random walks (RW)", "Ornstein-Uhlenbeck (OU) processes"), 2)
)

# ---- PLOT q_td ----
qtd <- ggplot(plot_data, aes(x = d, y = q_td)) +
  geom_line(aes(color = scenario, group = t), alpha = 0.3, linewidth = rel(0.6)) +
  geom_text(aes(x, y, label = label), df_labs, vjust = 1, hjust = 0, size = 3, fontface = "bold") +
  facet_grid(fr ~ model, scale = "free") + 
  # facet_wrap(~scenario, nrow = 2, scales = "free") + 
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = c("black", "darkblue", "red3", "seagreen4", "mediumpurple4", "darkgoldenrod3")) +  
  labs(title = NULL, x = TeX("$d$"), y = TeX("$q_t(d)$"), color = "Scenario") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +  
  # theme_classic(9) +  
  theme_bw(9) +
  theme(
    legend.position = "none",  
    # axis.title.x = element_text(size = 20, face = "bold"), 
    # axis.title.y = element_text(size = 20, face = "bold"),
    # strip.text = element_text(size = 20, face = "bold"), 
    # axis.text.x = element_text(size = 20),  
    # axis.text.y = element_text(size = 20) 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
qtd

ggsave(filename = file.path(path_proj, "plots_to_show", "sims_scenarios.png"),
       plot = qtd,
       width = 8, height = 4.5, dpi = 300)

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
