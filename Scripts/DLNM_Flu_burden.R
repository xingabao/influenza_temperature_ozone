# Load R packages
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))

# Plotting Theme
base.col <- '#000000'
base.size <- 18
base.family <- 'serif'

theme.com <- theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(family = base.family, color = 'black', size = base.size),
    axis.text = element_text(color = 'black', size = base.size * 0.8, family = base.family),
    axis.title.x = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold', margin = margin(t = 0.15, unit = 'in')),
    axis.title.y = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.text = element_text(size = base.size * 0.9, face = 'bold'),
    legend.position = 'right'
  )

theme.bw <- theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(family = base.family, color = 'black', size = base.size),
    axis.text = element_text(color = 'black', size = base.size * 0.8, family = base.family),
    axis.title.x = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold', margin = margin(t = 0.15, unit = 'in')),
    axis.title.y = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.text = element_text(size = base.size * 0.9, face = 'bold'),
    legend.position = 'right'
  )

# Calculate QAIC (Quasi-AIC)
get_qaic <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, lambda = fitted(model), log = TRUE))
  qaic <- -2 * loglik / phi + 2 * length(coef(model))
  return(qaic)
}

# Sensitivity Analysis (keep as 95th vs median for temp and O3, as in original)
run_sensitivity_check <- function(data, var_name, age_group_name, test_q = 0.95) {
  dfs_var <- 2:5
  dfs_lag <- 3:5
  max_lags <- c(4, 6, 8, 10, 12)
  
  data <- data[order(Year, Week)]
  if (!'Time' %in% colnames(data)) data[, Time := 1:.N]
  
  var_vec <- data[[var_name]]
  outcome <- data$Positive_Cases
  valid_idx <- !is.na(var_vec) & !is.na(outcome)
  data_sub <- data[valid_idx]
  var_vec <- var_vec[valid_idx]
  
  if (nrow(data_sub) < 50) return(NULL)
  
  if ('Total_Tests' %in% colnames(data_sub)) {
    offset_val <- log(data_sub$Total_Tests)
  } else {
    offset_val <- log(data_sub$Population_Base)
  }
  
  n_years <- as.numeric(difftime(max(data_sub$date), min(data_sub$date), units = 'days')) / 365.25
  df_time <- ceiling(7 * n_years)
  
  test_val <- quantile(var_vec, test_q, na.rm = TRUE)
  
  sens_results <- list()
  counter <- 1
  
  for (ml in max_lags) {
    for (dv in dfs_var) {
      for (dl in dfs_lag) {
        cb_sens <- crossbasis(
          var_vec, lag = ml,
          argvar = list(fun = 'ns', df = dv),
          arglag = list(fun = 'ns', df = dl)
        )
        
        form_str <- 'outcome ~ cb_sens + ns(Time, df = df_time) + Pandemic + offset(offset_val)'
        fit_sens <- tryCatch({
          glm(as.formula(form_str), family = quasipoisson(), data = data_sub)
        }, error = function(e) return(NULL))
        
        if (!is.null(fit_sens)) {
          qaic <- get_qaic(fit_sens)
          ref_val <- median(var_vec, na.rm = TRUE)
          cp_est <- crosspred(cb_sens, fit_sens, cen = ref_val, at = test_val)
          
          sens_results[[counter]] <- data.frame(
            AgeGroup = age_group_name,
            Variable = var_name,
            Max_Lag = ml,
            DF_Var = dv,
            DF_Lag = dl,
            QAIC = qaic,
            RR_Est = cp_est$allRRfit,
            RR_Low = cp_est$allRRlow,
            RR_High = cp_est$allRRhigh
          )
          counter <- counter + 1
        }
      }
    }
  }
  return(data.table::rbindlist(sens_results))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Data Loading and Pre-processing
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Macau Data
dat.raw.mo <- data.table::fread(file.path(dat.dir, 'Macau/INFLU-U2.csv.gz'), select = c('Age', 'Year', 'Week', 'FLUA', 'FLUB'))
dat.raw.mo <- dat.raw.mo[FLUA != -1 & FLUB != -1]
dat.raw.mo[, Positive := as.integer(FLUA == 1 | FLUB == 1)]

age.breaks.mo <- c(0, 6, 12, 18, 50, 65, Inf)
age.labels.mo <- c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 50)', '[50, 65)', '65+')
dat.raw.mo[, AgeGroup := cut(Age, breaks = age.breaks.mo, labels = age.labels.mo, right = FALSE)]

dat.week.mo <- dat.raw.mo[, .(Positive_Cases = sum(Positive, na.rm = TRUE), Total_Tests = .N), by = .(Year, Week, AgeGroup)]

dat.cli.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Climate.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.air.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Air.Quality.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.m.mo <- merge(dat.cli.mo, dat.air.mo[, .(Year, Week, O3 = IAQI_O3)], by = c('Year', 'Week'), all.x = TRUE)
dat.m.mo <- merge(dat.m.mo, dat.week.mo, by = c('Year', 'Week'), all.x = TRUE)

dat.m.mo[, Pandemic := 0]
dat.m.mo[date >= as.Date('2020-03-01') & date <= as.Date('2023-03-31'), Pandemic := 1]

dat.m.mo <- dat.m.mo[Year >= 2010 & Year <= 2025]
dat.m.mo <- dat.m.mo[!is.na(AgeGroup)] %>% na.omit()

# 1.2 Hong Kong Data
dat.raw.hk <- suppressMessages(readxl::read_excel(file.path(dat.dir, 'HK/flux_data.xlsx'), skip = 2))
dat.sub.hk <- dat.raw.hk[, c(1, 2)]
hk.admission.cols <- dat.raw.hk[, 17:22]
age.labels.hk <- c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 50)', '[50, 65)', '65+')
colnames(hk.admission.cols) <- age.labels.hk

dat.admission.hk <- cbind(dat.sub.hk, hk.admission.cols) %>% na.omit()
dat.admission.hk.long <- melt(
  as.data.table(dat.admission.hk),
  id.vars = c(names(dat.sub.hk)),
  variable.name = 'AgeGroup',
  value.name = 'Rate'
)

dat.admission.hk.long$Population_Base <- 10000
dat.admission.hk.long$Positive_Cases <- ceiling(dat.admission.hk.long$Rate * dat.admission.hk.long$Population_Base / 100)
colnames(dat.admission.hk.long)[1:2] <- c('Year', 'Week')
dat.admission.hk.long$Year <- as.numeric(dat.admission.hk.long$Year)

dat.cli.hk <- readxl::read_excel(file.path(dat.dir, 'HK/Climate.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table() %>%
  .[, .(ave.temp = mean(ave.temp, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]

dat.air.hk <- readxl::read_excel(file.path(dat.dir, 'HK/Air.Quality.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table() %>%
  .[, .(O3 = mean(IAQI_O3, na.rm = TRUE)), by = .(Year, Week)]

dat.m.hk <- merge(dat.cli.hk, dat.air.hk, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk <- merge(dat.m.hk, dat.admission.hk.long, by = c('Year', 'Week'), all.x = TRUE)

dat.m.hk[, Pandemic := 0]
dat.m.hk[date >= as.Date('2020-03-01') & date <= as.Date('2023-03-31'), Pandemic := 1]

dat.m.hk <- dat.m.hk[Year >= 2014 & Year <= 2025]
dat.m.hk <- dat.m.hk[!is.na(AgeGroup)]

# 1.3 England Data
dat.week.en <- readxl::read_excel(file.path(dat.dir, 'England/influenza-age.xlsx'))
dat.week.en <- dat.week.en %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.week.en <- dat.week.en %>%
  dplyr::select(Year, Week, Age, Positive_Rate = FLUAB)

dat.week.en <- dat.week.en %>%
  mutate(
    AgeGroup = case_when(
      Age == "00-04" ~ '[0, 6)',
      Age == "05-14" ~ '[6, 12)',
      Age == "15-44" ~ '[12, 18)',
      Age == "45-64" ~ '[18, 50)',
      Age == "65-79" ~ '[50, 65)',
      Age == "80+" ~ '65+',
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(Year, Week, AgeGroup, Positive_Rate)

dat.cli.en <- readxl::read_excel(file.path(dat.dir, 'England/Climate.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.air.en <- readxl::read_excel(file.path(dat.dir, 'England/Air.Quality.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.m.en <- merge(dat.cli.en, dat.air.en[, .(Year, Week, O3 = IAQI_O3)], by = c('Year', 'Week'), all.x = TRUE)
dat.m.en <- merge(dat.m.en, dat.week.en, by = c('Year', 'Week'), all.x = TRUE)

dat.m.en[, Pandemic := 0]
dat.m.en[date >= as.Date('2020-03-23') & date <= as.Date('2022-02-24'), Pandemic := 1]
dat.m.en <- dat.m.en[!is.na(AgeGroup)]
dat.m.en$Population_Base <- 10000
dat.m.en$Positive_Cases <- dat.m.en$Positive_Rate * 10000
dat.m.en$ave.temp <- (dat.m.en$ave.temp - 32) * 5/9

# 1.4 Zhuhai Data
dat.week.zh <- readxl::read_excel(file.path(dat.dir, 'Zhuhai/influenza-age.xlsx'))
dat.cli.zh <- readxl::read_excel(file.path(dat.dir, 'Zhuhai/Climate.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.air.zh <- readxl::read_excel(file.path(dat.dir, 'Zhuhai/Air.Quality.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()

dat.m.zh <- merge(dat.cli.zh, dat.air.zh[, .(Year, Week, O3 = IAQI_O3)], by = c('Year', 'Week'), all.x = TRUE)
dat.m.zh <- merge(dat.m.zh, dat.week.zh, by = c('Year', 'Week'), all.x = TRUE)

dat.m.zh[, Pandemic := 0]
dat.m.zh[date >= as.Date('2020-03-01') & date <= as.Date('2023-03-31'), Pandemic := 1]
dat.m.zh <- dat.m.zh[!is.na(AgeGroup)]
dat.m.zh <- dat.m.zh %>% select(Year, Week, AgeGroup, date, ave.temp, O3, Positive_Cases, Total_Tests, Pandemic) %>% na.omit()


if (TRUE) {
  # Sensitivity analysis
  dat.m.mo <- dat.m.mo %>% filter(Year < 2020 | Year > 2022)
  dat.m.hk <- dat.m.hk %>% filter(Year < 2020 | Year > 2022)
  dat.m.en <- dat.m.en %>% filter(Year < 2020 | Year > 2022)
  dat.m.zh <- dat.m.zh %>% filter(Year < 2020 | Year > 2022)
  ofig <- glue('{ofig}_remove_covid')
}
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. DLNM Analysis Function (returns both tails summaries)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
run_age_dlnm <- function(data, var_name, age_group_name, lag_weeks, city = 'xxx', n_sim = 1000) {
  
  data <- data[order(Year, Week)]
  data[, Time := 1:.N]
  
  var_vec <- data[[var_name]]
  outcome <- data$Positive_Cases
  valid_idx <- !is.na(var_vec) & !is.na(outcome) & outcome >= 0
  data_sub <- data[valid_idx]
  var_vec <- var_vec[valid_idx]
  
  if (length(var_vec) < 50) return(NULL)
  
  # Offset
  if ('Total_Tests' %in% colnames(data_sub)) {
    offset_val <- log(data_sub$Total_Tests)
  } else {
    offset_val <- log(data_sub$Population_Base)
  }
  
  # Season
  if (city == 'England') {
    if (!'Season' %in% colnames(data_sub)) {
      data_sub[, Season := ifelse(month(date) %in% c(10, 11, 12, 1, 2, 3, 4), 'Winter', 'Summer')]
    }
  } else {
    if (!'Season' %in% colnames(data_sub)) {
      data_sub[, Season := ifelse(month(date) %in% c(12, 1, 2, 3, 4), 'Winter', 'Summer')]
    }
  }
  
  argvar <- list(fun = 'ns', df = 3)
  arglag <- list(fun = 'ns', df = 3)
  
  cb_var <- crossbasis(var_vec, lag = lag_weeks, argvar = argvar, arglag = arglag)
  
  # Formula
  form_str <- 'outcome ~ cb_var + ns(Time, df = 5) + Season + Pandemic + offset(offset_val)'
  formula <- as.formula(form_str)
  
  fit <- tryCatch({
    glm(formula, family = quasipoisson(), data = data_sub)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  # ---------------------------------------------------------
  # A. Attributable Burden (split tails: <5th and >95th), ref = median
  # Uses attrdl() to properly account for lagged exposure history
  # ---------------------------------------------------------
  cen_val <- median(var_vec, na.rm = TRUE)
  x_cf_all <- rep(cen_val, length(var_vec))
  cases_obs <- outcome[valid_idx]
  total_cases <- sum(cases_obs, na.rm = TRUE)
  
  compute_an_allmedian <- function(x, x_cf, cb, fit) {
    # Rebuild crossbasis under observed and counterfactual series
    cb_obs <- crossbasis(
      x, lag = attr(cb, "lag"),
      argvar = attr(cb, "argvar"),
      arglag = attr(cb, "arglag")
    )
    cb_cf <- crossbasis(
      x_cf, lag = attr(cb, "lag"),
      argvar = attr(cb, "argvar"),
      arglag = attr(cb, "arglag")
    )
    
    # Pull the cb_var block of coefficients
    bname <- names(coef(fit))
    ind_cb <- grep("^cb_var", bname)
    if (length(ind_cb) == 0) ind_cb <- grep("cb_var", bname)
    if (length(ind_cb) == 0) stop("Could not find cb_var coefficients in model.")
    
    beta_cb <- coef(fit)[ind_cb]
    V_cb <- vcov(fit)[ind_cb, ind_cb, drop = FALSE]
    
    X_obs <- as.matrix(cb_obs)
    X_cf <- as.matrix(cb_cf)
    
    # Align columns just in case
    if (!identical(colnames(X_obs), colnames(X_cf))) {
      X_cf <- X_cf[, colnames(X_obs), drop = FALSE]
    }
    
    # Observed linear predictor already includes everything
    eta_obs <- as.numeric(predict(fit, type = "link"))
    
    # Swap only the cb contribution to get counterfactual eta
    cb_contrib_obs <- as.numeric(X_obs %*% beta_cb)
    cb_contrib_cf <- as.numeric(X_cf %*% beta_cb)
    
    eta_cf <- eta_obs - cb_contrib_obs + cb_contrib_cf
    
    mu_obs <- exp(eta_obs)
    mu_cf <- exp(eta_cf)
    
    an_est <- sum(mu_obs - mu_cf, na.rm = TRUE)
    
    # Simulated CI using MVN for cb coefficients
    set.seed(42)
    beta_sim <- mvtnorm::rmvnorm(n_sim, mean = beta_cb, sigma = V_cb)
    X_diff <- X_obs - X_cf
    delta_cb_sim <- X_diff %*% t(beta_sim) # n_time x n_sim
    
    an_sim <- apply(delta_cb_sim, 2, function(d) {
      sum(mu_obs - exp(eta_obs - d), na.rm = TRUE)
    })
    
    ci <- as.numeric(stats::quantile(an_sim, c(0.025, 0.975), na.rm = TRUE))
    
    list(an_est = an_est, an_low = ci[1], an_high = ci[2])
  }
  
  an_res <- compute_an_allmedian(var_vec, x_cf_all, cb_var, fit)
  
  burden_res <- data.frame(
    AgeGroup = age_group_name,
    Variable = var_name,
    Total_Cases = total_cases,
    Cen = cen_val,
    Attributable_Cases = an_res$an_est,
    AN_Low = an_res$an_low,
    AN_High = an_res$an_high,
    PAF = an_res$an_est / total_cases,
    PAF_Low = an_res$an_low / total_cases,
    PAF_High = an_res$an_high / total_cases
  )
  
  # ---------------------------------------------------------
  # B. Curve Prediction (keep 5-95% curve)
  # ---------------------------------------------------------
  pred_min <- quantile(var_vec, 0.05, na.rm = TRUE)
  pred_max <- quantile(var_vec, 0.95, na.rm = TRUE)
  var_pred <- seq(pred_min, pred_max, length.out = 50)
  
  cp <- crosspred(cb_var, fit, cen = cen_val, at = var_pred)
  
  res_df <- data.frame(
    Var = var_pred,
    RR = cp$allRRfit,
    Low = cp$allRRlow,
    High = cp$allRRhigh,
    AgeGroup = age_group_name,
    Variable = var_name,
    Reference = cen_val
  )
  
  # ---------------------------------------------------------
  # Summary at both tails: 5th vs median AND 95th vs median
  # ---------------------------------------------------------
  low_idx <- which.min(abs(var_pred - pred_min))
  high_idx <- which.min(abs(var_pred - pred_max))
  
  low_log_rr <- log(cp$allRRfit[low_idx])
  low_se_rr <- (log(cp$allRRhigh[low_idx]) - log(cp$allRRlow[low_idx])) / (2 * 1.96)
  
  high_log_rr <- log(cp$allRRfit[high_idx])
  high_se_rr <- (log(cp$allRRhigh[high_idx]) - log(cp$allRRlow[high_idx])) / (2 * 1.96)
  
  summary_res <- data.frame(
    AgeGroup = age_group_name,
    Variable = var_name,
    Low_RR = cp$allRRfit[low_idx],
    Low_Low = cp$allRRlow[low_idx],
    Low_High = cp$allRRhigh[low_idx],
    Low_logRR = low_log_rr,
    Low_SE = low_se_rr,
    High_RR = cp$allRRfit[high_idx],
    High_Low = cp$allRRlow[high_idx],
    High_High = cp$allRRhigh[high_idx],
    High_logRR = high_log_rr,
    High_SE = high_se_rr
  )
  
  return(list(Curve = res_df, Summary = summary_res, Burden = burden_res))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Execution
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
results.curve.list <- list()
results.summary.list <- list()
results.burden.list <- list()
results.sens.list <- list()

age.order.unified <- c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 50)', '[50, 65)', '65+')
region.order <- c('Macau', 'Hong Kong', 'Zhuhai', 'England')

# 3.1 Macau: sensitivity (for lag selection)
unique_ages_mo <- intersect(age.order.unified, unique(dat.m.mo$AgeGroup))

for (ag in unique_ages_mo) {
  dt_sub <- dat.m.mo[AgeGroup == ag]
  
  sens_res <- run_sensitivity_check(dt_sub, 'ave.temp', ag, test_q = 0.95)
  if (!is.null(sens_res)) {
    sens_res$Region <- 'Macau'
    results.sens.list[[paste0('Sens_Temp_', ag)]] <- sens_res
  }
  
  sens_res <- run_sensitivity_check(dt_sub, 'O3', ag, test_q = 0.95)
  if (!is.null(sens_res)) {
    sens_res$Region <- 'Macau'
    results.sens.list[[paste0('MO_O3_', ag)]] <- sens_res
  }
}

dat.sensitivity <- rbindlist(results.sens.list)

# ----------------------------------------------------------------------------- #
# Determine optimal lags SEPARATELY for Temperature and Ozone
# ----------------------------------------------------------------------------- #
# Select the best model (min QAIC) for each Variable-AgeGroup pair
best.models <- dat.sensitivity[order(QAIC), .SD[1], by = .(Region, AgeGroup, Variable)]

# Calculate the representative optimal lag for each variable
# Using median across age groups ensures robustness against outliers
optimal_lags <- best.models[, .(Optimal_Lag = as.integer(median(Max_Lag))), by = .(Variable)]

max_lag_temp <- optimal_lags[Variable == "ave.temp", Optimal_Lag]
max_lag_o3   <- optimal_lags[Variable == "O3", Optimal_Lag]

# Fallback in case one variable is missing (unlikely but safe)
if (is.na(max_lag_temp)) max_lag_temp <- 8
if (is.na(max_lag_o3))   max_lag_o3 <- 8

print(glue("Optimal Lag for Temperature: {max_lag_temp} weeks"))
print(glue("Optimal Lag for Ozone: {max_lag_o3} weeks"))

# Run region-wide
run_region <- function(region_name, dat_region, city_flag = "xxx", lag_temp, lag_o3) {
  unique_ages <- intersect(age.order.unified, unique(dat_region$AgeGroup))
  
  for (ag in unique_ages) {
    dt_sub <- dat_region[AgeGroup == ag]
    
    # Temperature Analysis with specific lag
    res <- run_age_dlnm(dt_sub, "ave.temp", ag, lag_temp, city = city_flag)
    if (!is.null(res)) {
      res$Curve$Region <- region_name
      res$Summary$Region <- region_name
      res$Burden$Region <- region_name
      
      results.curve.list[[paste0(region_name, "_Temp_", ag)]] <<- res$Curve
      results.summary.list[[paste0(region_name, "_Temp_Sum_", ag)]] <<- res$Summary
      results.burden.list[[paste0(region_name, "_Temp_Bur_", ag)]] <<- res$Burden
    }
    
    # Ozone Analysis with specific lag
    res <- run_age_dlnm(dt_sub, "O3", ag, lag_o3, city = city_flag)
    if (!is.null(res)) {
      res$Curve$Region <- region_name
      res$Summary$Region <- region_name
      res$Burden$Region <- region_name
      
      results.curve.list[[paste0(region_name, "_O3_", ag)]] <<- res$Curve
      results.summary.list[[paste0(region_name, "_O3_Sum_", ag)]] <<- res$Summary
      results.burden.list[[paste0(region_name, "_O3_Bur_", ag)]] <<- res$Burden
    }
  }
}

run_region('Macau', dat.m.mo, city_flag = 'xxx', lag_temp = max_lag_temp, lag_o3 = max_lag_o3)
run_region('Hong Kong', dat.m.hk, city_flag = 'xxx', lag_temp = max_lag_temp, lag_o3 = max_lag_o3)
run_region('Zhuhai', dat.m.zh, city_flag = 'xxx', lag_temp = max_lag_temp, lag_o3 = max_lag_o3)
run_region('England', dat.m.en, city_flag = 'England', lag_temp = max_lag_temp, lag_o3 = max_lag_o3)

# Combine Results
dat.curve <- rbindlist(results.curve.list)
dat.summary <- rbindlist(results.summary.list)
dat.burden <- rbindlist(results.burden.list)

# Set Factor Levels
dat.curve[, AgeGroup := factor(AgeGroup, levels = age.order.unified)]
dat.curve[, Region := factor(Region, levels = region.order)]

dat.summary[, AgeGroup := factor(AgeGroup, levels = age.order.unified)]
dat.summary[, Region := factor(Region, levels = region.order)]

dat.burden[, AgeGroup := factor(AgeGroup, levels = age.order.unified)]
dat.burden[, Region := factor(Region, levels = region.order)]

# England age group relabel for plotting consistency (as your original did)
dat.summary_plot <- copy(dat.summary)
dat.summary_plot[Region == 'England', AgeGroup := c('[0, 5)', '[0, 5)', '[5, 15)', '[5, 15)', '[15, 45)', '[15, 45)', '[45, 65)', '[45, 65)', '[65, 80)', '[65, 80)', '80+', '80+')]
dat.summary_plot[, AgeGroup := factor(
  AgeGroup,
  levels = c('[0, 5)', '[0, 6)', '[5, 15)', '[6, 12)', '[12, 18)', '[15, 45)', '[18, 50)', '[45, 65)', '[50, 65)', '65+', '[65, 80)', '80+')
)]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4. Build combined forest dataset (4 facets: low/high temp & low/high O3)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
make_forest_long <- function(dt_sum) {
  # dt_sum is data.table with columns:
  # Low_RR Low_Low Low_High High_RR High_Low High_High and Variable
  low_dt <- dt_sum[, .(
    Region, AgeGroup, Variable, Tail = 'Low',
    RR = Low_RR, CI_Low = Low_Low, CI_High = Low_High
  )]
  
  high_dt <- dt_sum[, .(
    Region, AgeGroup, Variable, Tail = 'High',
    RR = High_RR, CI_Low = High_Low, CI_High = High_High
  )]
  
  out <- rbind(low_dt, high_dt)
  
  out[, ExposureFacet := fifelse(Variable == 'ave.temp' & Tail == 'Low', 'Low temperature (5th vs median)',
                                 fifelse(Variable == 'ave.temp' & Tail == 'High', 'High temperature (95th vs median)',
                                         fifelse(Variable == 'O3' & Tail == 'Low', 'Low ozone (5th vs median)',
                                                 'High ozone (95th vs median)')))]
  
  out[, ExposureFacet := factor(
    ExposureFacet,
    levels = c('Low temperature (5th vs median)', 'High temperature (95th vs median)',
               'Low ozone (5th vs median)', 'High ozone (95th vs median)')
  )]
  
  out[, Sig_Label := fifelse(CI_Low > 1 | CI_High < 1, '*', '')]
  
  return(out[])
}

dat.forest <- make_forest_long(dat.summary_plot)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5. Visualization
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pals <- c(
  '[0, 6)' = '#E50611', '[6, 12)' = '#404040', '[12, 18)' = '#9D57B8',
  '[18, 50)' = '#E67E24', '[50, 65)' = '#19A955', '65+' = '#3D99DA',
  '[0, 5)' = '#E50611', '[5, 15)' = '#404040', '[15, 45)' = '#9D57B8',
  '[45, 65)' = '#E67E24', '[65, 80)' = '#19A955', '80+' = '#3D99DA'
)

facet_fill_cols <- c(
  'Low temperature (5th vs median)' = '#E50611',
  'High temperature (95th vs median)' = '#00A087',
  'Low ozone (5th vs median)' = '#E50611',
  'High ozone (95th vs median)' = '#00A087'
)

ppcols2 <- c('Temperature' = '#FF7F00', 'Ozone' = '#0530AD')
plcols <- c('Macau' = '#E50914', 'Hong Kong' = '#00A087', 'Zhuhai' = '#FFE25B', 'England' = '#50B6F6')

# Forest plot function: single region per file; facet by ExposureFacet only
forest_one_region <- function(region_name) {
  dt <- dat.forest %>% filter(Region == region_name)
  
  dt <- dt %>% mutate(
    FacetRow = case_when(
      ExposureFacet %in% c("Low temperature (5th vs median)", "Low ozone (5th vs median)") ~ "Low (5th vs median)",
      ExposureFacet %in% c("High temperature (95th vs median)", "High ozone (95th vs median)") ~ "High (95th vs median)",
      TRUE ~ NA_character_
    ),
    FacetCol = case_when(
      grepl("temperature", ExposureFacet, ignore.case = TRUE) ~ "Temperature",
      grepl("ozone|O3", ExposureFacet, ignore.case = TRUE) ~ "Ozone",
      TRUE ~ NA_character_
    ),
    FacetRow = factor(FacetRow, levels = c("Low (5th vs median)", "High (95th vs median)")),
    FacetCol = factor(FacetCol, levels = c("Temperature", "Ozone"))
  )
  
  gg <- ggplot(dt, aes(x = AgeGroup, y = RR, color = AgeGroup)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_text(aes(label = Sig_Label, y = CI_High + 0.10), size = 5, color = 'black', fontface = 'bold', nudge_x = 0.15, family = base.family) +
    coord_flip() +
    facet_grid(FacetRow ~ FacetCol, scales = "free_y") +
    scale_color_manual(values = pals) +
    labs(x = 'Age Group', y = 'Relative Risk') +
    theme.bw +
    theme(
      legend.position = 'none',
      strip.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      strip.background = element_rect(linewidth = 0),
      panel.grid = element_blank(),
      plot.title = element_blank()
    )
  
  # Manually tint facet strips (ggplot_gtable hack)
  ggt <- ggplot_gtable(ggplot_build(gg))
  
  stript <- which(grepl('strip-t', ggt$layout$name))
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', ggt$grobs[[i]]$grobs[[1]]$childrenOrder))
    ggt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(ppcols2[k], 1, 7), '7F')
    k <- k + 1
  }
  
  stripr <- which(grepl('strip-r', ggt$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', ggt$grobs[[i]]$grobs[[1]]$childrenOrder))
    ggt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(facet_fill_cols[k], 1, 7), '4F')
    k <- k + 1
  }
  
  return(ggt)
}

# 4.2 Attributable Burden Plot (PAF)
# Population Attributable Fraction (PAF)
# Positive: Risk Factor (Excess Cases); Negative: Protective Factor (Prevented Cases)
places <- c('Macau', 'Zhuhai')
places <- c('Hong Kong', 'England')
if ('Hong Kong' %in% places) {
  ooo = '_attributable_burden_plot_HK_EN.pdf'
} else {
  ooo = '_attributable_burden_plot.pdf'
}
gg.burden. <- ggplot(dat.burden %>% filter(Region %in% places), aes(x = AgeGroup, y = PAF, fill = AgeGroup)) +
  geom_bar(stat = 'identity', width = 0.88) +
  geom_text(
    aes(
      label = sprintf('%+.1f%%', PAF*100),
      y = ifelse(
        (Region == 'Macau' & AgeGroup == '65+' & Variable == 'ave.temp') |
          (Region == 'Macau' & AgeGroup == '[6, 12)'),
        PAF / 1.2,
        ifelse(
          (Region == 'Macau' & AgeGroup == '[18, 50)' & Variable == 'O3') |
            (Region == 'Macau' & AgeGroup == '[50, 65)' & Variable == 'O3') |
            (Region == 'Zhuhai' & AgeGroup == '65+' & Variable == 'O3') |
            (Region == 'Hong Kong' & AgeGroup == '[12, 18)' & Variable == 'O3'),
          -0.15,
          ifelse((Region == 'Hong Kong' & AgeGroup == '[12, 18)' & Variable == 'ave.temp'), 0.12, PAF)
        )
      ),
      vjust = ifelse(
        (Region == 'Macau' & AgeGroup == '65+' & Variable == 'ave.temp'),
        0.5,
        ifelse(PAF >= 0, -0.5, 1.5)
      ),
      color = ifelse(
        (Region == 'Macau' & AgeGroup == '[6, 12)') |
          (Region == 'Hong Kong' & AgeGroup == '[12, 18)' & Variable == 'ave.temp') |
          (Region == 'Hong Kong' & AgeGroup == '[12, 18)' & Variable == 'O3'), '#FFFFFF', base.col),
    ),
    size = 3.5
  ) +
  scale_color_identity() +
  facet_grid(Variable ~ Region, scales = 'free_y', labeller = labeller(Variable = c(ave.temp = 'Temperature', O3 = 'Ozone'))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = pals) +
  labs(x = 'Age Group', y = 'Attributable Fraction (%)') +
  theme.bw +
  theme(
    legend.position = 'none',
    strip.text = element_text(color = base.col, size = base.size * 0.9, family = base.family),
    strip.background = element_rect(linewidth = 0),
    axis.text.x = element_text(color = 'black', size = base.size * 0.8, family = base.family, angle = 45, hjust = 1),
    panel.grid = element_blank()
  ); gg.burden.

gg.burden <- ggplot_gtable(ggplot_build(gg.burden.))
stript <- which(grepl('strip-t', gg.burden$layout$name))
k <- 1
for (i in stript) {
  j <- which(grepl('rect', gg.burden$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.burden$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(plcols[places][k], 1, 7), '7F')
  k <- k + 1
}
stripr <- which(grepl('strip-r', gg.burden$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg.burden$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.burden$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(ppcols2[k], 1, 7), '4F')
  k <- k + 1
}

# 
ggsave(gg.burden, filename = glue('{fig.dir}/{ofig}{ooo}'), width = 9, height = 7.8, bg = NULL)

# 4.3 Faceted Exposure-Response Curves
# Temperature
# Temperature Effect on Influenza Risk (Harmonized Age Groups)
# Reference: Median Temperature (Vertical Dotted Line)
dat.curve <- dat.curve %>% mutate(AgeGroup = case_when(
  Region == 'England' & AgeGroup == '[0, 6)' ~ '[0, 5)',
  Region == 'England' & AgeGroup == '[6, 12)' ~ '[5, 15)',
  Region == 'England' & AgeGroup == '[12, 18)' ~ '[15, 45)',
  Region == 'England' & AgeGroup == '[18, 50)' ~ '[45, 65)',
  Region == 'England' & AgeGroup == '[50, 65)' ~ '[65, 80)',
  Region == 'England' & AgeGroup == '65+' ~ '80+',
  TRUE ~ AgeGroup
))
dat.curve$AgeGroup <- factor(dat.curve$AgeGroup, levels = c('[0, 5)', '[0, 6)', '[5, 15)', '[6, 12)', '[12, 18)', '[15, 45)', '[18, 50)', '[45, 65)', '[50, 65)', '65+', '[65, 80)', '80+'))

curve.temp <- function(cities = c(names(plcols)), ym = 3) {
  gg.curve.temp. <- ggplot(dat.curve %>% filter(Region %in% cities, Variable == 'ave.temp'), aes(x = Var, y = RR)) +
    geom_vline(aes(xintercept = Reference), linetype = 'dotted', color = 'grey50') +
    geom_ribbon(aes(ymin = Low, ymax = High, fill = AgeGroup), alpha = 0.2) +
    geom_line(aes(color = AgeGroup), linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    coord_cartesian(ylim = c(0.5, ym)) +
    scale_color_manual(values = pals) +
    scale_fill_manual(values = pals) +
    labs(x = expression(bold(paste('Temperature (', degree, 'C)'))), y = 'Relative Risk') +
    facet_grid(Region ~ AgeGroup) +
    theme.com +
    theme(
      legend.position = 'none',
      strip.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      strip.background = element_rect(linewidth = 0)
    )
  
  gg.curve.temp <- ggplot_gtable(ggplot_build(gg.curve.temp.))
  stript <- which(grepl('strip-t', gg.curve.temp$layout$name))
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', gg.curve.temp$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.curve.temp$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(pals[k], 1, 7), '7F')
    k <- k + 1
  }
  stripr <- which(grepl('strip-r', gg.curve.temp$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg.curve.temp$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.curve.temp$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(plcols[cities][k], 1, 7), '4F')
    k <- k + 1
  }
  return(gg.curve.temp)
}

# Ozone
# Ozone Effect on Influenza Risk (Harmonized Age Groups)
# Reference: Median Ozone Level
curve.o3 <- function(cities = c(names(plcols)), ym = 3) {
  gg.curve.o3. <- ggplot(dat.curve %>% filter(Region %in% cities, Variable == 'O3'), aes(x = Var, y = RR)) +
    geom_vline(aes(xintercept = Reference), linetype = 'dotted', color = 'grey50') +
    geom_ribbon(aes(ymin = Low, ymax = High, fill = AgeGroup), alpha = 0.2) +
    geom_line(aes(color = AgeGroup), linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    coord_cartesian(ylim = c(0.5, ym)) +
    scale_color_manual(values = pals) +
    scale_x_continuous(limits = c(5, 40)) +
    scale_fill_manual(values = pals) +
    labs(x = 'Ozone (IAQI)', y = 'Relative Risk') +
    facet_grid(Region ~ AgeGroup) +
    theme.com +
    theme(
      legend.position = 'none',
      strip.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
      strip.background = element_rect(linewidth = 0)
    )
  
  gg.curve.o3 <- ggplot_gtable(ggplot_build(gg.curve.o3.))
  stript <- which(grepl('strip-t', gg.curve.o3$layout$name))
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', gg.curve.o3$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.curve.o3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(pals[k], 1, 7), '7F')
    k <- k + 1
  }
  stripr <- which(grepl('strip-r', gg.curve.o3$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg.curve.o3$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.curve.o3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(plcols[cities][k], 1, 7), '4F')
    k <- k + 1
  }
  return(gg.curve.o3)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6. Save Output (one forest file per region)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for (rg in region.order) {
  g <- forest_one_region(rg)
  ggsave(g, filename = glue('{fig.dir}/{ofig}_forest_{gsub(" ", "_", rg)}.pdf'), width = 8, height = 6.5, bg = NULL)
}

ggsave(curve.temp(cities = c('Macau', 'Hong Kong'), ym = 3), filename = glue('{fig.dir}/{ofig}_temp_curves_MO_HK.pdf'), width = 10, height = 7, bg = NULL)
ggsave(curve.temp(cities = c('Zhuhai'), ym = 2), filename = glue('{fig.dir}/{ofig}_temp_curves_ZH.pdf'), width = 10, height = 3.5, bg = NULL)
ggsave(curve.temp(cities = c('England'), ym = 5), filename = glue('{fig.dir}/{ofig}_temp_curves_EN.pdf'), width = 10, height = 3.5, bg = NULL)

ggsave(curve.o3(cities = c('Macau', 'Hong Kong'), ym = 3), filename = glue('{fig.dir}/{ofig}_o3_curves_MO_HK.pdf'), width = 10, height = 7, bg = NULL)
ggsave(curve.o3(cities = c('Zhuhai'), ym = 2), filename = glue('{fig.dir}/{ofig}_o3_curves_ZH.pdf'), width = 10, height = 3.5, bg = NULL)
ggsave(curve.o3(cities = c('England'), ym = 4), filename = glue('{fig.dir}/{ofig}_o3_curves_EN.pdf'), width = 10, height = 3.5, bg = NULL)

# 4.4 Sensitivity Analysis, Boxplot of RR stability
# Stability of Effect Estimates
# RR (at 95th percentile) across different Lags (x-axis) and DF (colors)
dat.sensitivity$AgeGroup <- factor(dat.sensitivity$AgeGroup, levels = age.order.unified)

gg.sensitivity. <- ggplot(dat.sensitivity, aes(x = as.factor(Max_Lag), y = RR_Est)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(DF_Var)), width = 0.2, alpha = 0.6) +
  facet_grid(Variable ~ AgeGroup, scales = 'free_y', labeller = labeller(Variable = c(ave.temp = 'Temperature', O3 = 'Ozone'))) +
  scale_color_manual(values = c('#E50611', '#404040', '#3D99DA', '#E67E24', '#19A955')) +
  labs(
    x = 'Maximum Lag (Weeks)',
    y = 'Estimated Relative Risk',
    color = 'DF of Variable'
  ) +
  theme.com +
  theme(
    legend.position = 'none',
    strip.text = element_text(color = base.col, size = base.size * 0.9, family = base.family),
    strip.background = element_rect(linewidth = 0)
  )

gg.sensitivity <- ggplot_gtable(ggplot_build(gg.sensitivity.))
stript <- which(grepl('strip-t', gg.sensitivity$layout$name))
k <- 1
for (i in stript) {
  j <- which(grepl('rect', gg.sensitivity$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.sensitivity$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(pals[k], 1, 7), '7F')
  k <- k + 1
}
stripr <- which(grepl('strip-r', gg.sensitivity$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg.sensitivity$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.sensitivity$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(ppcols2[k], 1, 7), '4F')
  k <- k + 1
}
ggsave(gg.sensitivity, filename = glue('{fig.dir}/{ofig}_sensitivity_plot.pdf'), width = 10, height = 7, bg = NULL)

# ============================================================================= #
# 7. Tables (Formatted for Publication)
# ============================================================================= #
# Define a helper function for formatting: Estimate (Low-High)
# Handles NA values gracefully
fmt_ci <- function(est, low, high) {
  ifelse(is.na(est), NA_character_, 
         sprintf("%.3f (%.3f - %.3f)", est, low, high))
}

# -----------------------------------------------------------------------------
# 7.1 Summary Table (Effect Estimates)
# -----------------------------------------------------------------------------
dat.summary_export <- as.data.frame(dat.summary) %>%
  mutate(
    # Format Low Exposure (5th vs Median)
    `Low Exposure (RR, 95% CI)` = fmt_ci(Low_RR, Low_Low, Low_High),
    # Format High Exposure (95th vs Median)
    `High Exposure (RR, 95% CI)` = fmt_ci(High_RR, High_Low, High_High)
  ) %>%
  # Select and reorder columns for publication
  select(
    Region,
    `Age Group` = AgeGroup,
    Variable,
    `Low Exposure (RR, 95% CI)`,
    `High Exposure (RR, 95% CI)`
  ) %>%
  # Optional: Recode variable names for better readability
  mutate(Variable = case_when(
    Variable == "ave.temp" ~ "Temperature",
    Variable == "O3" ~ "O3",
    TRUE ~ Variable
  ))

writexl::write_xlsx(dat.summary_export, glue('{tbl.dir}/Table_S_DLNM_Summary_LowHigh.xlsx'))

# -----------------------------------------------------------------------------
# 7.2 Burden Table (Attributable Burden)
# -----------------------------------------------------------------------------
dat.burden_export <- as.data.frame(dat.burden) %>%
  mutate(
    # Format PAF (Population Attributable Fraction)
    `PAF (95% CI)` = fmt_ci(PAF, PAF_Low, PAF_High)
  ) %>%
  # Select and reorder columns for publication
  select(
    Region,
    `Age Group` = AgeGroup,
    Variable,
    `PAF (95% CI)`
  ) %>%
  # Optional: Recode variable names
  mutate(Variable = case_when(
    Variable == "ave.temp" ~ "Temperature",
    Variable == "O3" ~ "O3",
    TRUE ~ Variable
  ))

writexl::write_xlsx(dat.burden_export, glue('{tbl.dir}/Table_S_DLNM_Burden.xlsx'))

# -----------------------------------------------------------------------------
# 7.3 QAIC Table (from Sensitivity Analysis)
# -----------------------------------------------------------------------------
# QAIC table from sensitivity
Table_S_QAIC <- dat.sensitivity %>%
  group_by(Region, AgeGroup, Variable) %>%
  slice_min(order_by = QAIC, n = 1) %>%
  select(Region, `Age Group` = AgeGroup, Variable, Max_Lag, DF_Var, DF_Lag, QAIC) %>%
  mutate(Variable = case_when(
    Variable == "ave.temp" ~ "Average Temperature",
    TRUE ~ Variable
  ))
writexl::write_xlsx(Table_S_QAIC, glue('{tbl.dir}/Table_S_QAIC_Comparison.xlsx'))

gc()
