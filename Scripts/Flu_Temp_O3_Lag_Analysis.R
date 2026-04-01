# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(patchwork)))
# 
# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Absolute Humidity
calculate_ah <- function(temp, rh) {
  e_s <- 6.112 * exp((17.67 * temp) / (temp + 243.5))
  e <- e_s * (rh / 100)
  ah <- (216.7 * e) / (temp + 273.15)
  return(ah)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Load data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
path_hk <- glue('{dat.dir}/HK/FLU-CL-AQ.xlsx')
path_macau <- glue('{dat.dir}/Macau/FLU-CL-AQ.xlsx')

process_data <- function(file_path, city_name) {
  
  # Read the Excel file
  df <- readxl::read_excel(file_path)
  
  # Check for the date column
  if ('date' %in% names(df)) {
    df$date <- as.Date(df$date)
  } else {
    stop("Error: 'date' column not found in the Excel file.")
  }
  
  df$month <- as.numeric(format(df$date, '%m'))
  
  # Ensure temperature and humidity are numeric
  df$ave.temp <- suppressWarnings(as.numeric(as.character(df$ave.temp)))
  df$humidity <- suppressWarnings(as.numeric(as.character(df$humidity)))
  
  # Calculate Absolute Humidity (AH)
  if ('ave.temp' %in% names(df) && 'humidity' %in% names(df)) {
    df$abs.humidity <- calculate_ah(df$ave.temp, df$humidity)
  } else {
    warning("Missing 'ave.temp' or 'humidity' columns. Absolute Humidity cannot be calculated!")
    df$abs.humidity <- NA
  }
  
  # Handle binary variables (impute NA with 0)
  if ('holiday' %in% names(df)) df$holiday <- ifelse(is.na(df$holiday), 0, df$holiday)
  if ('pandemic' %in% names(df)) df$pandemic <- ifelse(is.na(df$pandemic), 0, df$pandemic)
  
  # Define continuous variables requiring interpolation
  continuous_cols <- c('FLUA', 'FLUB', 'FLUAB', 'ave.temp', 'humidity', 'abs.humidity', 'PM10', 'PM2.5', 'NO2', 'SO2', 'O3')
  
  # Intersect to prevent errors if specific columns are missing in the dataset
  continuous_cols <- intersect(continuous_cols, names(df))
  
  # Linear Interpolation for missing values
  for (col in continuous_cols) {
    if (sum(!is.na(df[[col]])) > 2) {
      df[[col]] <- zoo::na.approx(df[[col]], rule = 2) 
    }
  }
  
  # Save raw data
  df_raw <- df
  df_raw$city <- city_name
  
  df_final <- df
  df_final$city <- city_name
  
  # Create a numeric time index for trend control
  df_final$time_index <- 1:nrow(df_final)
  
  print(paste('Successfully processed:', city_name, '| Sample size:', nrow(df_final)))
  return(list(processed = df_final, raw = df_raw))
}

# 加载数据
data_hk <- process_data(path_hk, 'Hong Kong')
data_macau <- process_data(path_macau, 'Macau')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Calculate Lag-Response (Core Upgrade: Multivariable Adjustment)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
calculate_lag_response_adjusted <- function(data, target_var, driver_var, max_lag = 14) {
  
  if (!target_var %in% names(data) || !driver_var %in% names(data)) {
    print(paste('Error: Variable', target_var, 'or', driver_var, 'not found. Skipping calculation.'))
    return(data.frame())
  }
  
  results <- data.frame()
  n <- nrow(data)
  
  # Identify available control variables
  has_ah <- 'abs.humidity' %in% names(data)
  has_hol <- 'holiday' %in% names(data)
  has_pan <- 'pandemic' %in% names(data)
  has_time <- 'time_index' %in% names(data)
  
  for (lag in 0:max_lag) {
    idx_start <- lag + 2 
    if (idx_start > n) break
    
    # 1. Construct core vectors
    y_vec <- data[[target_var]][idx_start:n]                   
    x_vec <- data[[driver_var]][(idx_start - lag):(n - lag)]
    
    # 2. Construct control variable matrix (Z)
    z_auto <- data[[target_var]][(idx_start - 1):(n - 1)]
    z_df <- data.frame(Auto = z_auto)
    
    # Add Absolute Humidity (same lag as driver)
    if (has_ah) {
      z_df$AH <- data$abs.humidity[(idx_start - lag):(n - lag)]
    }
    
    # Add Holiday (current week)
    if (has_hol) {
      z_df$Holiday <- data$holiday[idx_start:n]
    }
    
    # Add Pandemic (current week)
    if (has_pan) {
      z_df$Pandemic <- data$pandemic[idx_start:n]
    }
    
    # Add Time Trend Control
    if (has_time) {
      z_df$Time_Trend <- data$time_index[idx_start:n]
    }
    
    # 3. Data alignment and cleaning
    temp_df <- cbind(Y = y_vec, X = x_vec, z_df)
    temp_df <- na.omit(temp_df) 
    
    len <- nrow(temp_df)
    
    if (len > 30) { 
      tryCatch({
        clean_x <- temp_df$X
        clean_y <- temp_df$Y
        clean_z <- temp_df[, 3:ncol(temp_df)] 
        
        test <- ppcor::pcor.test(x = clean_x, y = clean_y, z = clean_z, method = 'pearson')
        
        k <- ncol(clean_z)
        se <- 1 / sqrt(len - 2 - k) 
        
        results <- rbind(results, data.frame(
          Lag = lag,
          Correlation = test$estimate,
          P_value = test$p.value,
          LowerCI = test$estimate - 1.96 * se,
          UpperCI = test$estimate + 1.96 * se
        ))
      }, error = function(e) { message(paste('Error at lag', lag, ':', e$message)) })
    }
  }
  return(results)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Analysis Loop
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
target_var_name <- 'FLUAB'
driver_list <- c('ave.temp', 'O3')
max_lag_weeks <- 14

# Use lapply to iterate over drivers
plot_list <- lapply(driver_list, function(driver_var_name) {
  
  print(paste('Processing driver:', driver_var_name))
  
  # A. Main Analysis
  res_hk <- calculate_lag_response_adjusted(data_hk$processed, target_var_name, driver_var_name, max_lag = max_lag_weeks)
  res_macau <- calculate_lag_response_adjusted(data_macau$processed, target_var_name, driver_var_name, max_lag = max_lag_weeks)
  
  if (nrow(res_hk) > 0) res_hk$City <- 'Hong Kong'
  if (nrow(res_macau) > 0) res_macau$City <- 'Macau'
  df_plot_full <- rbind(res_hk, res_macau)
  
  # B. Winter/Spring Only
  target_months <- c(12, 1, 2, 3, 4)
  hk_seasonal <- data_hk$processed %>% filter(month %in% target_months)
  macau_seasonal <- data_macau$processed %>% filter(month %in% target_months)
  
  res_hk_sea <- calculate_lag_response_adjusted(hk_seasonal, target_var_name, driver_var_name, max_lag = max_lag_weeks)
  if (nrow(res_hk_sea) > 0) {
    res_hk_sea$City <- 'Hong Kong'
    res_hk_sea$Type <- 'Winter/Spring Only'
  }
  res_macau_sea <- calculate_lag_response_adjusted(macau_seasonal, target_var_name, driver_var_name, max_lag = max_lag_weeks)
  if (nrow(res_macau_sea) > 0) {
    res_macau_sea$City <- 'Macau'
    res_macau_sea$Type <- 'Winter/Spring Only'
  }
  
  # Prepare sensitivity data
  if (nrow(res_hk) > 0) {
    res_hk_full <- res_hk 
    res_hk_full$City <- 'Hong Kong'
    res_hk_full$Type <- 'Full Year'
  }
  if (nrow(res_macau) > 0) {
    res_macau_full <- res_macau
    res_macau_full$City <- 'Macau'
    res_macau_full$Type <- 'Full Year'
  }
  
  df_sensitivity <- rbind(res_hk_full, res_macau_full, res_hk_sea, res_macau_sea)
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # 4. Plot one
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_plot_full <- df_plot_full %>%
    mutate(
      sig_label = case_when(
        P_value < 0.001 ~ '***',
        P_value < 0.01  ~ '**',
        P_value < 0.05  ~ '*',
        TRUE            ~ ''
      )
    )
  
  df_plot_full$City <- factor(df_plot_full$City, levels = c('Macau', 'Hong Kong'))
  
  # Dynamic Y-axis label
  y_label <- glue('Partial Correlation ({driver_var_name})')
  
  gga <- ggplot(df_plot_full, aes(x = Lag, y = Correlation, color = City, fill = City)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.15, color = NA) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_text(
      data = subset(df_plot_full, sig_label != ''), 
      aes(
        label = sig_label, 
        y = Correlation + ifelse(grepl('Hong Kong', City), 0.01, -0.04),
      ), 
      size = 6, 
      color = 'black', 
      show.legend = FALSE,
      vjust = 0
    ) +
    scale_color_manual(values = c('Hong Kong' = '#00A087', 'Macau' = '#E50914')) +
    scale_fill_manual(values = c('Hong Kong' = '#00A087', 'Macau' = '#E50914')) +
    scale_x_continuous(breaks = seq(0, max_lag_weeks, 1), expand = c(0, 0), limits = c(-0.1, 14.1)) +
    labs(x = '', y = y_label) +
    theme_bw(base_size = base.size, base_family = base.family) +
    theme(
      legend.position = 'inside',
      legend.position.inside = c(0.85, 0.15),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.box = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0.5, 0, unit = 'cm'),
      axis.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      axis.title.x = element_text(color = base.col, size = base.size * 0.8, face = 'bold', family = base.family, margin = margin(t = 0, unit = 'cm')),
      axis.title.y = element_text(color = base.col, size = base.size * 0.8, face = 'bold', family = base.family),
      plot.title = element_text(hjust = 0.5, face = 'bold')
    )
  
  Table_S_LagCorr <- df_plot_full %>%
    rename(
      `Time Lag (Weeks)` = Lag,
      `Partial Correlation Coefficient` = Correlation,
      `P-value` = P_value,
      `Lower 95% CI` = LowerCI,
      `Upper 95% CI` = UpperCI,
      `City` = City,
      `Significance Level` = sig_label
    )
  writexl::write_xlsx(Table_S_LagCorr, glue('{tbl.dir}/Table_S_Lagged_Partial_Correlations_{driver_var_name}.xlsx'))
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # 5. Plot two
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_sensitivity$City <- factor(df_sensitivity$City, levels = c('Macau', 'Hong Kong'))
  ggb <- ggplot(df_sensitivity, aes(x = Lag, y = Correlation, color = Type)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    geom_point(data = subset(df_sensitivity, P_value < 0.05), shape = 8, size = 2, color = 'black') +
    facet_wrap(~City,) +
    scale_color_manual(values = c('Full Year' = 'gray60', 'Winter/Spring Only' = '#DC0000')) +
    scale_x_continuous(breaks = seq(0, max_lag_weeks, 1), expand = c(0, 0), limits = c(-0.1, 14.1)) +
    labs(
      x = 'Lag (Weeks)',
      y = 'Partial Correlation (Adjusted)'
    ) +
    theme_bw(base_size = base.size, base_family = base.family) +
    theme(
      legend.position = 'inside',
      legend.position.inside = c(0.85, 0.15),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.box = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      axis.title.x = element_text(color = base.col, size = base.size * 0.8, face = 'bold', family = base.family, margin = margin(t = 0.3, unit = 'cm')),
      axis.title.y = element_text(color = base.col, size = base.size * 0.8, face = 'bold', family = base.family)
    )
  
  # Modify strip background colors
  hh <- ggplot_gtable(ggplot_build(ggb))
  stripr <- which(grepl('strip-t', hh$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', hh$grobs[[i]]$grobs[[1]]$childrenOrder))
    hh$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(c('Macau' = '#E50914', 'Hong Kong' = '#00A087')[k], 1, 7), '7F')
    k <- k + 1
  }
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # 6. Combine Plot
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  design <- '
    AAA
    BBB
    '
  gg_combined <- gga + free(wrap_elements(hh)) + 
    plot_layout(design = design, heights = c(0.8, 1))
  
  # Return the combined plot object
  return(gg_combined)
})

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 7. Save Outputs
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save each plot corresponding to the driver list
for (i in seq_along(driver_list)) {
  driver_name <- driver_list[i]
  current_plot <- plot_list[[i]]
  output_filename <- glue('{fig.dir}/{ofig}_{driver_name}.pdf')
  ggsave(current_plot, filename = output_filename, width = 10, height = 7, units = 'in', bg = '#FFFFFF')
  print(paste('Saving plot for:', driver_name, '->', output_filename))
}
