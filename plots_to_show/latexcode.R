library(dplyr)

library(dplyr)

library(dplyr)

library(dplyr)

highlight_metrics <- function(tables, method_names = NULL, date_labels = NULL) {
  # 确保输入是列表且长度为4
  # if (!is.list(tables) || length(tables) != 4) {
  #   stop("Input must be a list of 4 data frames")
  # }
  
  # 如果没有提供日期向量，使用默认数字
  if (is.null(date_labels)) {
    date_labels <- paste0("Scenario ", 1:4)
  }
  
  # 如果没有提供method_names，使用默认名称
  if (is.null(method_names)) {
    method_names <- c("Method 1", "Method 2", "Method 3", "Method 4", "Method 5")
  }
  
  # 合并表格并添加场景列和自定义方法名
  combined_df <- do.call(rbind, lapply(seq_along(tables), function(i) {
    df <- tables[[i]]
    df$Scenario <- i
    # 替换Method列为自定义方法名
    df$Method <- method_names[1:nrow(df)]
    return(df)
  }))
  
  # 定义需要标记的列
  cols_to_process <- c("RMSE", "RMSPE", "MAE", "MAPE", "Interval_Width")
  
  # 处理每个列
  highlighted_df <- combined_df %>%
    group_by(Scenario) %>%
    mutate(across(all_of(cols_to_process), 
                  ~ ifelse(. == min(.), 
                           paste0("\\textcolor{red}{", formatC(., format = "f", digits = 5), "}"),
                           formatC(., format = "f", digits = 5)))) %>%
    # Coverage Rate特殊处理（最大值标记）
    mutate(`Coverage_Rate` = 
             ifelse(`Coverage_Rate` == max(`Coverage_Rate`), 
                    paste0("\\textcolor{red}{", formatC(`Coverage_Rate`, format = "f", digits = 5), "}"),
                    formatC(`Coverage_Rate`, format = "f", digits = 5)))
  
  # 生成LaTeX表格
  latex_table <- "\\begin{table}[htbp]\n\\centering\n"
  latex_table <- paste0(latex_table, "\\caption{Metrics Comparison}\n")
  latex_table <- paste0(latex_table, "\\begin{tabular}{c|c|c|c|c|c|c|c}\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  latex_table <- paste0(latex_table, "Scenario & RMSE & RMSPE & MAE & MAPE & Interval Width & Coverage Rate & Method \\\\\n")
  latex_table <- paste0(latex_table, "\\hline\n")
  
  # 添加每个场景的数据
  for(i in seq_along(unique(highlighted_df$Scenario))) {
    scenario <- unique(highlighted_df$Scenario)[i]
    scenario_data <- highlighted_df[highlighted_df$Scenario == scenario,]
    
    # 添加日期行
    latex_table <- paste0(latex_table, "\\multicolumn{8}{l}{\\textbf{Now is ", date_labels[i], "}} \\\\\n")
    latex_table <- paste0(latex_table, "\\hline\n")
    
    # 添加该场景的数据行
    for(j in 1:nrow(scenario_data)) {
      row <- scenario_data[j,]
      row_str <- paste(
        "", # 第一列留空
        row$RMSE,
        row$RMSPE,
        row$MAE,
        row$MAPE,
        row$`Interval_Width`,
        row$`Coverage_Rate`,
        row$Method,
        sep = " & "
      )
      latex_table <- paste0(latex_table, row_str, " \\\\\n")
    }
    
    # 在场景之间添加分隔线
    latex_table <- paste0(latex_table, "\\hline\n")
  }
  
  latex_table <- paste0(latex_table, "\\end{tabular}\n")
  latex_table <- paste0(latex_table, "\\end{table}")
  
  return(latex_table)
}

scoreRange[c(2,3,4,6)]

cat(highlight_metrics(list_table, method_names = c("Fixed q", "Fixed b", "Linear b", "Spline b"),
                      date_labels = scoreRange[c(2,3,4,6)]))

scoreRange[c(3,5)]

cat(highlight_metrics(hus_list, method_names = c("Fixed q", "Fixed b", "Linear b", "Spline b"),
                      date_labels = scoreRange[c(3,5)]))
