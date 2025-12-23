library(ggplot2)
library(dplyr)
library(scales) 

# --- Data Preparation ---
# Assuming matrix_cases is correctly loaded and the first row is:
#matrix_cases_row1 <- c(22603, 73716, 104034, 118006, 127744, 134735, 139322, 142852, 144742, 146727, 148479)
data_raw <- c(117, 125, 39, 35, 44, 16, 12, 19, 43, 6, 6, 13, 7, 7, 3, 11, 0, 3, 1, 1, 0)
matrix_cases_row1 <- cumsum(data_raw)
final_reported_cases <- tail(matrix_cases_row1, 1) # 148479

plot_data <- data.frame(
  Delay_Days = 0:20,
  Cumulative_Cases = matrix_cases_row1
)

percentage_milestones <- c(0,0.25, 0.50, 0.75, 0.90, 0.95,1.00)

# Define the scaling factor for the secondary axis
# The formula is: Secondary_Value = (Primary_Value / Final_Value) * 100
# Since the scale is 0 to 148479 on the left, and 0 to 1 on the right (for percentage as ratio)
# The ratio needed is Final_Value (148479)
ratio <- final_reported_cases / 1 # 1 is the max percentage (100%)

# --- ggplot Visualization ---
p1 <- ggplot(plot_data, aes(x = Delay_Days, y = Cumulative_Cases)) +
  
  geom_hline(yintercept = percentage_milestones * final_reported_cases, 
             linetype = "dashed", # Use dashed line
             color = "grey60", 
             linewidth = 0.5,
             alpha = 0.7) +
  
  # 2. Lines for trends
  geom_line(color = "#0072B2", linewidth = 1.2) +
  
  # 3. Points for actual data (keep the points for clarity)
  geom_point(color = "#D55E00", size = 3, shape = 21, fill = "white", stroke = 1) +
  
  # 4. Data labels for cases (removed percentage to reduce clutter)
  geom_text(aes(label = formatC(Cumulative_Cases, format = "d", big.mark = ",")),
            vjust = -0.5, # Move text slightly above the point
            size = 3.5,
            color = "#000000") + 
  
  # Set titles and axis labels
  labs(
    title = "", 
    x = "Reporting Delay (In weeks)",
    y = "Cumulative Reported Cases (Count)" # Clarify left axis
  ) +
  
  # 5. Customize Axes
  scale_x_continuous(breaks = seq(min(plot_data$Delay_Days), max(plot_data$Delay_Days), by = 1)) +
  
  # Primary Y-axis (Cases)
  scale_y_continuous(
    labels = scales::comma,
    
    # Secondary Y-axis (Percentage)
    sec.axis = sec_axis(
      trans = ~ . / ratio, # Transformation: Cases / Final_Cases
      name = "Percentage of Final Reported Cases",
      labels = scales::percent, # Format labels as percentage
      breaks = percentage_milestones # *** Use the new percentage values as breaks ***
    )
  ) +
  
  # Use classic theme
  theme_classic() +
  
  # Fine-tune theme elements
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold", color = "#0072B2"), # Left axis color
    axis.title.y.right = element_text(size = 16, face = "bold", color = "#D55E00"), # Right axis color
    axis.text = element_text(size = 10),
    legend.position = "none",
    plot.margin = unit(c(1, 1.5, 1, 1), "cm") # Increased right margin for secondary axis title
  )

print(p1)


#get_date(week=36,year=2009)
