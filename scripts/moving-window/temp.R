#nowcast 

# alpha
N_up <- 30    # 
N_down <- 30  # 

# 
alpha_up <- seq(10, 300, length.out = N_up)^2 / 100  # control the speed

# 
alpha_down <- rev(alpha_up)

# 
alpha <- c(alpha_up, alpha_down)
alpha


seed <- 123
# length for nowcast
N_obs <- 60
D <- 15
k_back <- 10  

data <- simsDataGenQ( #alpha = alpha,
                      alpha =  c(1:10, seq(1, 197, by = 4)),
                      days = N_obs, method = "random_walk",
                      b = 0.1, sigma_rw = 0.1, seed = seed)

now <- as.Date("2024-02-01")
when <- seq(now-k_back+1, now, by="1 day")

case_data <- data$case_reported



########### loop to fit ##########

seeds = 123
scoreRange <- seq(as.Date("2024-01-20"),as.Date("2024-02-28"),by="1 day")

plot_list <- list()
list_i <- c(40)
i = 40
#c(1,5,10,15,20,25,30,35,40)
for (i in list_i) {
  #What's "today"
  now <- scoreRange[i]
  # show the status
  cat(paste("====================\nnow=",now,
            " (",i,"/",length(scoreRange),")\n====================\n",sep=""))
  
  #when <- seq(now-k-safePredictLag+1, now-safePredictLag, by="1 day")
  
  # Nowcast #
  # get data
  data_use <- head(data$case_reported, N_obs - length(scoreRange) + i)
  
  # truncated version of data
  data_trunc <- create_triangular_data(data_use, if_zero = F)
  N_obs_local <- nrow(data_trunc)
  
  # coordinates for non-NAs
  indices_data_trunc <- find_non_na_coords(data_trunc)
  
  data_trunc[is.na(data_trunc)] <- 0 # to avoid NAs in data
  
  X_spline <- create_basis(N_obs_local, n_knots = 5)

  # input list
  stan_data_trunc <- list(N_obs = N_obs_local, D = D + 1, Y = data_trunc,
                          K = nrow(indices_data_trunc), obs_index = indices_data_trunc,
                          J = ncol(X_spline),
                          X_spline = X_spline,
                          sigma_b = 0.1)
  
  fit_trunc <- stan(
    file = file.path(path_proj, "source", "models",
                     "trunc", "stan_model_dep-spline-trunc-modi.stan"),  
    data = stan_data_trunc, 
    iter = 2000, chains = 3, seed = seeds,
    #control = list(adapt_delta = 0.96, max_treedepth = 15),
    refresh = 500
  )
  
  
  fit_trunc_fixped_q <- stan(
    file = file.path(path_proj, "source", "models", "trunc",
                     "stan_model_time-ind-p-trunc.stan"),  
    data = stan_data_trunc, 
    iter = 2000, chains = 3, seed = seeds,
    refresh = 500
  )
  
  # extract parameters
  samples_nt <- rstan::extract(fit_trunc, pars = "N_t")$N_t
  samples_nt_fixped_q <- rstan::extract(fit_trunc_fixped_q, pars = "N_t")$N_t
  
  nowcasts <- data.frame(mean = apply(samples_nt, 2, mean),
                         lower = apply(samples_nt, 2, quantile, probs = 0.025),
                         upper = apply(samples_nt, 2, quantile, probs = 0.975),
                         #
                         mean_fixped_q = apply(samples_nt_fixped_q, 2, mean),
                         lower_fixped_q = apply(samples_nt_fixped_q, 2, quantile, probs = 0.025),
                         upper_fixped_q = apply(samples_nt_fixped_q, 2, quantile, probs = 0.975),
                         #
                         date =  data$date[1:N_obs_local],
                         case_true = data$case_true[1:N_obs_local],
                         case_reported = apply(data_trunc[1:N_obs_local,], 1, max))
  
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
    labs(title = NULL,
         x = NULL,
         y = "Cases / Nowcast",
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
  print(p)
  
  plot_list[[i]] <- p
}

################
#diagnose
summary(fit_trunc)$summary[, "Rhat"]
summary(fit_trunc)$summary[,"n_eff"]
stan_trace(fit_trunc, pars = c("b_t"))
stan_rdump(fit_trunc)



fit_trunc_fixped_q


#################

filenames <- c("nowcast1.png", "nowcast2.png", "nowcast3.png", "nowcast4.png")
count_temp <- 1
for (i in list_i) {
  ggsave(
    filename = filenames[count_temp],           # file name
    plot = plot_list[[i]],                 # object to save
    path = file.path(path_proj, "poster"),  # save dir
    width = 10,                        # width
    height = 8,                        # height
    units = "in",                      # size, "in", "cm", "mm"
    dpi = 300                          # DPI
  )
  count_temp <- count_temp + 1 
}

### note the environment
data$case_reported
data$date
data_paper_form <- dataTransform(data$case_reported, start_date = as.Date("2024-02-29"))

now_date <- as.Date("2024-01-20")
nc <- nowcast(now=now_date,when=now_date,
              dEventCol="dHosp",dReportCol="dReport",data=data_paper_form,D=D,method="lawless")
plotReportingTriangle(nc)

# plot for abstarct
ggplot(nowcasts, aes(x = date)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.5) +
  # geom_line(aes(y = mean), color = "blue") +
  geom_vline(xintercept = as.Date(now - D), color = "orange") +
  
  geom_line(aes(y = case_true, color = "Real Cases")) +
  geom_line(aes(y = case_reported, color = "Reported Cases")) +
  
  scale_color_manual(values = c("Real Cases" = "red", "Reported Cases" = "black")) +
  annotate("text", x = as.Date(now - D), y = 0, label = format(as.Date(now - D), "%b %d"),
           vjust = 1, hjust = -0.1, color = "orange") +
  labs(
    # title = "Nowcast with True and Reported Cases",
    x = NULL,
    y = "Cases",
    color = NULL  # legend title
  ) +
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

