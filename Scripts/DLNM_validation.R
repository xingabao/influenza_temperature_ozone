# Load R packages
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(patchwork)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'
gc()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data Loading and Preprocessing
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Macau
dat.raw <- data.table::fread(glue('{dat.dir}/Macau/INFLU-U2.csv.gz'), select = c('Year', 'Week', 'FLUA', 'FLUB'))
dat.raw <- dat.raw[FLUA != -1 & FLUB != -1]
dat.week <- dat.raw[, .(
  Total_Tests = .N,
  Positive_Cases = sum(FLUA == 1 | FLUB == 1, na.rm = TRUE)
), by = .(Year, Week)]
dat.week[, Pos_Rate := Positive_Cases / Total_Tests]
rm(dat.raw); gc()

dat.cli.mo <- readxl::read_excel(glue('{dat.dir}/Macau/Climate.xlsx')) %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()
dat.air.mo <- readxl::read_excel(glue('{dat.dir}/Macau/Air.Quality.xlsx')) %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.m.mo <- merge(dat.cli.mo, dat.air.mo[, .(Year, Week, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO)], by = c('Year', 'Week'), all.x = TRUE)
dat.m.mo <- merge(dat.m.mo, dat.week, by = c('Year', 'Week'), all.x = TRUE)

dat.m.mo <- na.omit(dat.m.mo, cols = c('ave.temp', 'Positive_Cases', 'Total_Tests', 'humidity', 'pressure', 'PM2.5', 'PM10', 'O3', 'SO2', 'NO2', 'CO'))
dat.m.mo <- dat.m.mo[Year >= 2010 & Year <= 2025]
dat.m.mo[, Time := as.numeric(date - min(date)) / 7 + 1]
dat.m.mo[is.na(PM2.5), PM2.5 := mean(PM2.5, na.rm = TRUE)]
dat.m.mo[is.na(PM10), PM10 := mean(PM10, na.rm = TRUE)]
dat.m.mo[is.na(O3), O3 := mean(O3, na.rm = TRUE)]
dat.m.mo[is.na(SO2), SO2 := mean(SO2, na.rm = TRUE)]
dat.m.mo[is.na(NO2), NO2 := mean(NO2, na.rm = TRUE)]
dat.m.mo[is.na(CO), CO := mean(CO, na.rm = TRUE)]
dat.m.mo[, Season := ifelse(month(date) %in% c(12, 1, 2, 3, 4), 'Winter/Spring', 'Summer/Autumn')]
dat.m.mo$Pandemic = 0
dat.m.mo[dat.m.mo$date > '2020-01-01' & dat.m.mo$date < '2023-01-01', ]$Pandemic = 1
dat.m.mo <- dat.m.mo %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.mo <- dat.m.mo %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))

# Hong Kong
dat.raw.hk <- readxl::read_excel(glue('{dat.dir}/HK/flux_data.xlsx'), skip = 2)
dat.flux.hk <- data.table(
  Year = as.numeric(dat.raw.hk[[1]]),
  Week = as.numeric(dat.raw.hk[[2]]),
  Positive_Cases = as.numeric(dat.raw.hk[[10]]),
  Pos_Rate = as.numeric(dat.raw.hk[[14]])
) %>% mutate(
  Total_Tests = ceiling(Positive_Cases / Pos_Rate)
)
dat.flux.hk <- na.omit(dat.flux.hk, cols = c('Year', 'Week', 'Positive_Cases', 'Total_Tests'))
dat.flux.hk[, Date := as.Date(paste(Year, Week, 1, sep = '-'), format = '%Y-%U-%u')]

dat.cli.hk <- readxl::read_excel(glue('{dat.dir}/HK/Climate.xlsx')) %>% 
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()

dat.air.hk <- readxl::read_excel(glue('{dat.dir}/HK/Air.Quality.xlsx')) %>% 
  dplyr::select(date, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()

dat.cli.hk <- dat.cli.hk[, .(ave.temp = mean(ave.temp, na.rm = TRUE), humidity = mean(humidity, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]
dat.air.hk <- dat.air.hk[, .(O3 = mean(O3, na.rm = TRUE), PM2.5 = mean(PM2.5, na.rm = TRUE), PM10 = mean(PM10, na.rm = TRUE), CO = mean(CO, na.rm = TRUE), SO2 = mean(SO2, na.rm = TRUE), NO2 = mean(NO2, na.rm = TRUE)), by = .(Year, Week)]

dat.m.hk <- merge(dat.flux.hk, dat.cli.hk, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk <- merge(dat.m.hk, dat.air.hk, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk <- na.omit(dat.m.hk, cols = c('ave.temp', 'Positive_Cases', 'Total_Tests', 'PM2.5', 'O3'))
dat.m.hk <- dat.m.hk[Year >= 2014 & Year <= 2025] 
dat.m.hk[, Time := as.numeric(Date - min(Date)) / 7 + 1]
dat.m.hk[, Season := ifelse(month(Date) %in% c(12, 1, 2, 3, 4), 'Winter/Spring', 'Summer/Autumn')]
dat.m.hk$Pandemic = 0
dat.m.hk[dat.m.hk$date > '2020-01-01' & dat.m.hk$date < '2023-01-01', ]$Pandemic = 1
dat.m.hk <- dat.m.hk %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.hk <- dat.m.hk %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))

# Zhuhai
dat.raw.zh <- readxl::read_excel(glue('{dat.dir}/Zhuhai/FLU-CL-AQ.xlsx'))
dat.raw.zh <- data.table(
  Date = as.Date(dat.raw.zh$date),
  Positive_Cases = dat.raw.zh$nFLUAB,
  Total_Tests = dat.raw.zh$nSample
) %>% mutate(
  Pos_Rate = Positive_Cases / Total_Tests,
  Year = year(Date),
  Week = week(Date)
)
dat.raw.zh <- na.omit(dat.raw.zh, cols = c('Year', 'Week', 'Positive_Cases', 'Total_Tests'))

dat.cli.zh <- readxl::read_excel(glue('{dat.dir}/Zhuhai/Climate.xlsx')) %>% 
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()

dat.air.zh <- readxl::read_excel(glue('{dat.dir}/Zhuhai/Air.Quality.xlsx')) %>% 
  dplyr::select(date, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()

dat.cli.zh <- dat.cli.zh[, .(ave.temp = mean(ave.temp, na.rm = TRUE), humidity = mean(humidity, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]
dat.air.zh <- dat.air.zh[, .(O3 = mean(O3, na.rm = TRUE), PM2.5 = mean(PM2.5, na.rm = TRUE), PM10 = mean(PM10, na.rm = TRUE), CO = mean(CO, na.rm = TRUE), SO2 = mean(SO2, na.rm = TRUE), NO2 = mean(NO2, na.rm = TRUE)), by = .(Year, Week)]

dat.m.zh <- merge(dat.raw.zh, dat.cli.zh, by = c('Year', 'Week'), all.x = TRUE)
dat.m.zh <- merge(dat.m.zh, dat.air.zh, by = c('Year', 'Week'), all.x = TRUE)
dat.m.zh <- na.omit(dat.m.zh, cols = c('ave.temp', 'Positive_Cases', 'Total_Tests', 'PM2.5', 'O3'))
dat.m.zh <- dat.m.zh[Year >= 2015 & Year <= 2025] 
dat.m.zh[, Time := as.numeric(Date - min(Date)) / 7 + 1]
dat.m.zh[, Season := ifelse(month(Date) %in% c(12, 1, 2, 3, 4), 'Winter/Spring', 'Summer/Autumn')]
dat.m.zh$Pandemic = 0
dat.m.zh[dat.m.zh$date > '2020-01-01' & dat.m.zh$date < '2023-01-01', ]$Pandemic = 1
dat.m.zh <- dat.m.zh %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.zh <- dat.m.zh %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))

# ------------------------------------------------------------------------------
# Define Core DLNM Analysis Function
# ------------------------------------------------------------------------------
run_dlnm_analysis <- function(data, var_name, outcome_type = 'counts', subset_name = 'Full Year', lag_weeks = 8) {
  
  message(sprintf('  [Run] Variable: %-10s | Subset: %s | Max Lag: %d weeks', var_name, subset_name, lag_weeks))
  var_vec <- data[[var_name]]
  
  outcome <- data$Positive_Cases
  offset_val <- log(data$Total_Tests)
  control_formula <- '+ Season * Pandemic' 
 
  valid_idx <- !is.na(var_vec) & !is.na(outcome)
  data_sub <- data[valid_idx]
  var_vec <- var_vec[valid_idx]
  
  if (length(var_vec) < 50) return(NULL)
  argvar <- list(fun = 'ns', knots = quantile(var_vec, probs = c(0.1, 0.5, 0.9), na.rm = TRUE))
  arglag <- list(fun = 'ns', df = 3)
  cb_var <- crossbasis(var_vec, lag = lag_weeks, argvar = argvar, arglag = arglag)
  form_str <- paste0('outcome ~ cb_var + ns(Time, df = 5)', control_formula)
  
  if (!is.null(offset_val)) form_str <- paste0(form_str, ' + offset(offset_val)')
  formula <- as.formula(form_str)
  fit <- tryCatch({
    glm(formula, family = quasipoisson(), data = data_sub, model = TRUE)
  }, error = function(e) { return(NULL) })
  
  if (is.null(fit)) return(NULL)
  var_pred <- seq(quantile(var_vec, 0.05, na.rm = TRUE), quantile(var_vec, 0.95, na.rm = TRUE), length.out = 50)
  cp <- crosspred(cb_var, fit, cen = median(var_vec, na.rm = TRUE), at = var_pred)
  critical_val <- NA
  crit_type <- 'Risk Threshold'
  if (var_name == 'O3') {
    idx_min <- which.min(cp$allRRfit)
    critical_val <- var_pred[idx_min]
    crit_type <- 'Protection Peak'
  } else {
    diff_RR <- diff(cp$allRRfit); diff_var <- diff(var_pred)
    slope <- diff_RR / diff_var
    var_slope <- var_pred[-length(var_pred)] + diff_var/2
    idx_low <- which(var_slope < median(var_vec, na.rm = TRUE))
    if (length(idx_low) > 0) critical_val <- var_slope[idx_low[which.min(slope[idx_low])]]
  }
  return(list(
    df = data.frame(Var = var_pred, RR = cp$allRRfit, Low = cp$allRRlow, High = cp$allRRhigh),
    var_name = var_name, subset = subset_name, critical = critical_val, type = crit_type
  ))
}

# ------------------------------------------------------------------------------
# Execute Main Analysis across all regions and core variables
# ------------------------------------------------------------------------------
results_all <- list()
vars_core <- c('ave.temp', 'pressure', 'humidity', 'Particulates', 'O3', 'Primary Pollutants')
for (v in vars_core) {
  res <- run_dlnm_analysis(dat.m.mo, v, 'counts', 'Macau_Full', lag_weeks = 8)
  if (!is.null(res)) results_all[[paste0('Macau_', v, '_Full')]] <- res
}

for (v in vars_core) {
  res <- run_dlnm_analysis(dat.m.hk, v, 'counts', 'HK_Full', lag_weeks = 8)
  if (!is.null(res)) results_all[[paste0('HK_', v, '_Full')]] <- res
}

for (v in vars_core) {
  res <- run_dlnm_analysis(dat.m.zh, v, 'counts', 'ZH_Full', lag_weeks = 8)
  if (!is.null(res)) results_all[[paste0('ZH_', v, '_Full')]] <- res
}

# ------------------------------------------------------------------------------
# Sensitivity Analysis: Test varying maximum lag structures
# ------------------------------------------------------------------------------
run_sens_maxlag <- function(data, var_name, max_lag) {
  var_vec <- data[[var_name]]
  outcome <- data$Positive_Cases
  offset_val <- log(data$Total_Tests)
  argvar <- list(fun = 'ns', df = 4)
  arglag <- list(fun = 'ns', df = 3)
  cb <- crossbasis(var_vec, lag = max_lag, argvar = argvar, arglag = arglag)
  fit <- glm(outcome ~ cb + ns(Time, df = 5), family = quasipoisson(), offset = offset_val, data = data)
  pred_var <- seq(min(var_vec, na.rm = TRUE), max(var_vec, na.rm = TRUE), length.out = 50)
  cp <- crosspred(cb, fit, cen = median(var_vec, na.rm = TRUE), at = pred_var)
  res_df <- data.frame(Var = pred_var, RR = cp$allRRfit, Low = cp$allRRlow, High = cp$allRRhigh)
  res_df$Max_Lag <- paste0('Lag ', max_lag, ' Weeks')
  res_df$Variable <- var_name
  
  return(res_df)
}

# Validation of PLS-PM Pathways
val.plot <- rbindlist(lapply(results_all, function(x) {
  d <- x$df
  d$Variable <- x$var_name
  d$Region <- ifelse(grepl('HK', x$subset), 'Hong Kong', ifelse(grepl('Macau', x$subset), 'Macau', 'Zhuhai'))
  d$Role <- case_when(
    x$var_name == 'ave.temp' ~ 'Temperature',
    x$var_name == 'pressure' ~ 'Atmosphere',
    x$var_name == 'humidity' ~ 'Moisture',
    x$var_name == 'Particulates' ~ 'Particulates',
    x$var_name == 'Primary Pollutants' ~ 'Primary Pollutants',
    x$var_name == 'O3' ~ 'Ozone'
  )
  d
}))

# Common ggplot2 theme
theme.com <- theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(family = base.family, color = base.col, size = base.size),
    axis.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
    axis.title = element_text(color = base.col, size = base.size * 0.9, family = base.family, face = 'bold')
  )

pals <- c(
  'Temperature' = '#E50611', 
  'Moisture' = '#404040',
  'Atmosphere' = '#9D57B8',
  'Particulates' = '#E67E24', 
  'Primary Pollutants' = '#19A955',
  'Ozone' = '#3D99DA'
)
for (place in c('Macau', 'HK', 'Zhuhai')) {
  
  region <- ifelse(place == 'HK', 'Hong Kong', place)
  
  dat.m. <- val.plot[Region == region]
  dat.m.$Role <- factor(dat.m.$Role, levels = names(pals))
  gg.val. <- ggplot(dat.m., aes(x = Var, y = RR)) +
    geom_ribbon(aes(ymin = Low, ymax = High), fill = 'grey70', alpha = 0.4) +
    geom_line(aes(color = Role), linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, ifelse(place == 'Macau', 5, ifelse(place == 'Zhuhai', 4, 10)))) +
    facet_wrap(~Role, scales = 'free_x') +
    scale_color_manual(values = pals) +
    labs(title = NULL, x = 'Exposure Level', y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'none',
      strip.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      strip.background = element_rect(linewidth = 0)
    ); gg.val.
  
  gg.val <- ggplot_gtable(ggplot_build(gg.val.))
  stripr <- which(grepl('strip-t', gg.val$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg.val$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.val$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(pals[k], 1, 7), '7F')
    k <- k + 1
  }
  plot(gg.val)
  ggsave(gg.val, filename = glue('{fig.dir}/{ofig}_{place}.pdf'), width = 9, height = 7.2, bg = NULL)
  
  # Sensitivity Analysis
  sens_maxlag_list <- list()
  lags_to_test <- 3:9
  if (place == 'Macau') {
    dat.m. <- dat.m.mo
  } else if (place == 'Zhuhai') {
    dat.m. <- dat.m.zh
  } else {
    dat.m. <- dat.m.hk
  }
  
  for (v in c('ave.temp', 'pressure', 'humidity', 'Particulates', 'O3', 'Primary Pollutants')) {
    for (lg in lags_to_test) {
      sens_maxlag_list[[paste0(v, '_lag', lg)]] <- run_sens_maxlag(dat.m., v, lg)
    }
  }
  df_sens_maxlag <- rbindlist(sens_maxlag_list)
  df_sens_maxlag[, Max_Lag := factor(Max_Lag, levels = paste0('Lag ', lags_to_test, ' Weeks'))]
  
  # Plot
  gg.sens.temp <- ggplot(df_sens_maxlag[Variable == 'ave.temp'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_y_continuous(expand = c(0.1, 0.1), labels = c('0.0', '0.5', '1.0', '1.5', '2.0', '2.5'), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5), limits = c(0, ceiling(max(df_sens_maxlag[Variable == 'ave.temp']$RR)))) +
    labs(title = NULL, x = expression(bold(paste('Temperature (', degree, 'C)'))),  y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'none'
    ); gg.sens.temp
  
  gg.sens.hum <- ggplot(df_sens_maxlag[Variable == 'humidity'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_y_continuous(expand = c(0.1, 0.1), labels = c('0.0', '0.5', '1.0', '1.5', '2.0', '2.5'), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5), limits = c(0, 5)) +
    labs(title = NULL, x = 'Humidity (%)', y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'none'
    ); gg.sens.hum
  
  gg.sens.pressure <- ggplot(df_sens_maxlag[Variable == 'pressure'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_y_continuous(expand = c(0.1, 0.1), limits = c(0, ifelse(place == 'HK', 8, 3))) +
    labs(title = NULL, x = 'Pressure (hPa)', y = 'Relative Risk (RR)') +
    theme.com +
    theme(
      legend.position = 'none'
    ); gg.sens.pressure
  
  gg.sens.pm <- ggplot(df_sens_maxlag[Variable == 'Particulates'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 20)) +
    labs(title = NULL, x = 'Particulates', y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'none'
    ); gg.sens.pm
  
  gg.sens.pp <- ggplot(df_sens_maxlag[Variable == 'Primary Pollutants'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_x_continuous(limits = c(0, 110)) +
    labs(title = NULL, x = 'Primary Pollutants', y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'none'
    ); gg.sens.pp
  
  gg.sens.o3 <- ggplot(df_sens_maxlag[Variable == 'O3'], aes(x = Var, y = RR, color = Max_Lag)) +
    geom_line(linewidth = 1) + geom_hline(yintercept = 1, linetype = 'dashed') +
    scale_color_manual(values = c('#E50914', '#00A087', '#FDBF6F', '#FF7F00', '#33A02C', '#A6CEE3', '#08306B')) + 
    scale_y_continuous(expand = c(0., 0.), labels = c('0.0', '1.0', '2.0', '3.0', '4.0', '5.0'), limits = c(0, 5.5)) +
    labs(title = NULL, x = 'Ozone', y = 'Relative Risk') +
    theme.com +
    theme(
      legend.position = 'inside',
      legend.position.inside = c(0.65, 0.65),
      legend.title = element_blank(),
      legend.background = element_blank()
    ); gg.sens.o3
  
  gg.combine <- (gg.sens.temp + gg.sens.hum + gg.sens.pressure) /
    (gg.sens.pp + gg.sens.pm + gg.sens.o3) +
    plot_annotation(
      tag_levels = 'a',
      tag_prefix = '(',
      tag_suffix = ')'
    ) & theme(
      plot.tag = element_text(color = base.col, size = base.size * 1.0, family = base.family, face = 'bold'),
    ); gg.combine
  
  ggsave(gg.combine, filename = glue('{fig.dir}/{ofig}_sensitivity_MaxLag_{place}.pdf'), width = 12, height = 7, bg = NULL)
}

# Regional Validation: Consistency of Driver and Protector Effects
dat.compare <- val.plot[Variable %in% c('ave.temp', 'pressure', 'humidity', 'Particulates', 'O3', 'Primary Pollutants')]
dat.compare$Variable <- factor(dat.compare$Variable, levels = c('ave.temp', 'humidity', 'pressure', 'Particulates', 'Primary Pollutants', 'O3'))
dat.compare$Region <- factor(dat.compare$Region, levels = c('Macau', 'Hong Kong', 'Zhuhai'))
gg.compare <- ggplot(dat.compare, aes(x = Var, y = RR, color = Region, fill = Region)) +
  geom_ribbon(aes(ymin = Low, ymax = High), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  facet_wrap(
    ~ Variable,
    scales = 'free_x', 
    labeller = labeller(Variable = c('ave.temp' = 'Temperature', 'O3' = 'Ozone', 'humidity' = 'Moisture', 'pressure' = 'Atmosphere'))
  ) +
  scale_color_manual(values = c('Macau' = '#e41a1c', 'Hong Kong' = '#1f78b4', 'Zhuhai' = '#00A087')) +
  scale_fill_manual(values = c('Macau' = '#e41a1c', 'Hong Kong' = '#1f78b4', 'Zhuhai' = '#00A087')) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = NULL, x = 'Exposure Level', y = 'Relative Risk') +
  theme.com +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.28),
    legend.title = element_blank(),
    legend.background = element_blank(),
    strip.text = element_text(color = base.col, size = base.size * 0.9, family = base.family),
    strip.background = element_blank()
  ); gg.compare

if (FALSE) ggsave(gg.compare, filename = glue('{fig.dir}/{ofig}_comparison.pdf'), width = 10, height = 7, bg = NULL)

# Print summary
message('\n关键临界点摘要：')
for (v in vars_core) {
  key <- paste0('Macau_', v, '_Full')
  if (key %in% names(results_all)) {
    res <- results_all[[key]]
    val <- res$critical
    type <- res$type
    if (!is.na(val)) message(sprintf('  - %-10s: %.2f (%s)', v, val, type))
  }
}
