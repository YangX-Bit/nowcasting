library(aweek)
library(dplyr)
library(readr)

get_infodengue_data <- function(root_dir, start, end, states, D = 0, 
                                week_start = "Sunday", fill_missing = TRUE) {
  
  # 参数验证函数
  validate_inputs <- function() {
    if (!dir.exists(root_dir)) {
      stop("Root directory does not exist: ", root_dir)
    }
    if (!is.character(states) || length(states) == 0) {
      stop("states must be a non-empty character vector")
    }
    if (!is.numeric(D) || D < 0 || D != as.integer(D)) {
      stop("D must be a non-negative integer")
    }
    if (!is.logical(fill_missing) || length(fill_missing) != 1) {
      stop("fill_missing must be a single logical value")
    }
  }
  
  is_ew_format <- function(x) {
    grepl("^\\d{6,7}$", as.character(x))
  }
  
  is_date_format <- function(x) {
    d <- try(as.Date(x), silent=TRUE)
    !inherits(d, "try-error") && !is.na(d)
  }
  
  date_to_ew <- function(date_obj) {
    aw_obj <- date2week(date_obj, week_start = "Sunday")
    year <- substr(aw_obj, 1, 4)
    weeknum <- substr(aw_obj, 7, 8)
    paste0(year, weeknum)
  }
  
  normalize_ew <- function(x) {
    if (is_ew_format(x)) {
      return(as.character(x))
    } else if (is_date_format(x)) {
      return(date_to_ew(as.Date(x)))
    } else {
      stop("Invalid start/end format. Use YYYYWW or YYYY-MM-DD")
    }
  }
  
  # 将"YYYYWW"转换为对应的周起始日期
  get_date_from_ew <- function(ew_str) {
    year <- substr(ew_str,1,4)
    weeknum <- substr(ew_str,5,6)
    weeknum <- sprintf("%02d", as.integer(weeknum))
    aw_code_day <- paste0(year, "-W", weeknum, "-1")
    aw_obj <- as.aweek(aw_code_day, week_start = week_start)
    get_date(aw_obj)
  }
  
  # 验证输入
  validate_inputs()
  
  set_week_start(week_start)
  
  start_ew <- normalize_ew(start)
  end_ew <- normalize_ew(end)
  
  #start_date <- get_date_from_ew(start_ew)
  #end_date <- get_date_from_ew(end_ew)
  start_date <- get_date(year = substr(start_ew, 1, 4), week = substr(start_ew, 5, 6))
  end_date <- get_date(year = substr(end_ew, 1, 4), week = substr(end_ew, 5, 6))
  
  if (end_date < start_date) {
    stop("End date cannot be earlier than start date")
  }
  
  # 生成从start到end的周序列
  all_dates <- seq(from = start_date, to = end_date, by = "week")
  all_ews <- sapply(all_dates, date_to_ew, USE.NAMES = FALSE)
  
  # 如果D>0则扩展week序列
  if (D > 0) {
    extra_dates <- seq(from = end_date + 7, by = "week", length.out = D)
    extra_ews <- sapply(extra_dates, date_to_ew, USE.NAMES = FALSE)
    all_ews_extended <- c(all_ews, extra_ews)
    all_dates_extended <- c(all_dates, extra_dates)
  } else {
    all_ews_extended <- all_ews
    all_dates_extended <- all_dates
  }
  
  N <- length(all_ews)  # 矩阵行数
  total_weeks <- length(all_ews_extended) # 扩展后周数
  data_log <- list()
  
  result_list <- vector("list", length(states))
  names(result_list) <- states
  
  for (st_idx in seq_along(states)) {
    st <- states[st_idx]
    
    # 初始化矩阵
    M <- matrix(NA, nrow = N, ncol = D+1)
    
    # 遍历每一周的数据文件
    for (i in seq_len(total_weeks)) {
      ew_code <- all_ews_extended[i]
      current_ew_start_date <- format(all_dates_extended[i], "%Y-%m-%d")
      # 延迟一周发布
      delayed_date <- as.Date(current_ew_start_date) + 7
      file_name <- paste0(st, "_", format(delayed_date, "%Y-%m-%d"), "_infodengue.csv")
      folder_path <- file.path(root_dir, ew_code)
      file_path <- file.path(folder_path, file_name)
      
      if (file.exists(file_path)) {
        df <- read_csv(file_path, show_col_types = FALSE)
        # df中应有ew和sum_of_cases列
        if ("ew" %in% names(df) && "sum_of_cases" %in% names(df)) {
          for (row_i in seq_len(nrow(df))) {
            W_code <- df$ew[row_i]
            cases_val <- df$sum_of_cases[row_i]
            W_idx <- match(W_code, all_ews)
            if (!is.na(W_idx)) {
              col <- i - W_idx + 1
              if (col >= 1 && col <= (D+1)) {
                M[W_idx, col] <- cases_val
              }
            }
          }
        }
        
        data_log[[length(data_log) + 1]] <- list(
          state = st,
          ew = ew_code,
          file = file_name,
          status = "Found"
        )
      } else {
        data_log[[length(data_log) + 1]] <- list(
          state = st,
          ew = ew_code,
          file = file_name,
          status = "Missing"
        )
      }
    }
    
    # 缺失值填补：从左到右填补
    if (fill_missing) {
      for (r in 1:N) {
        for (c in 2:(D+1)) {
          if (is.na(M[r,c])) {
            M[r,c] <- M[r,c-1]
          }
        }
      }
    }
    
    # 根据要求：从倒数第一列(delayD)开始，到倒数第D列(delay1)，
    # 第delayX列末尾X行置为NA
    # delay0列(第一列)不需要填NA，因为X=0
    for (X in 1:D) {
      # 列为X+1列（因为X从0开始映射列索引= X+1）
      # 需要将末尾X行设为NA
      start_row <- N - X + 1
      end_row <- N
      M[start_row:end_row, X+1] <- NA
    }
    
    # all_dates 对应每一行的实际日期
    rownames(M) <- format(all_dates, "%Y-%m-%d")
    # 设置列名为 delay0 ... delayD
    colnames(M) <- paste0("delay", 0:D)
    
    result_list[[st]] <- M
  }
  
  structure(result_list, data_log = data_log)
}

brazil_ufs <- c(
  "AC", "AL", "AP",
  "AM",
  "BA", "CE", "DF", "ES",
  "GO",
  "MA", "MT", "MS", "MG", "PA",
  "PB", "PR", "PE",
  "PI",
  "RJ", "RN", "RS", "RO", "RR",
  "SC",
  "SP", "SE", "TO"
)


