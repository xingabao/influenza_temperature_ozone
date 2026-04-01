# Load R packages
suppressMessages(suppressWarnings(library(mice)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(purrr)))
suppressMessages(suppressWarnings(library(qgraph)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(psychonetrics)))

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
    axis.title = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.text = element_text(size = base.size * 0.9, face = 'bold'),
    legend.position = 'right'
  )

add_holiday_flag <- function(dat_week, dat_holiday) {
  holiday_vec <- as.Date(unique(dat_holiday$date))
  week_flags <- dat_week %>%
    dplyr::select(Year, Week, date) %>%
    distinct() %>%
    mutate(temp_daily_date = map(date, ~seq.Date(.x - 6, .x, by = "day"))) %>%
    unnest(temp_daily_date) %>%
    group_by(Year, Week) %>%
    summarise(
      has_holiday = any(temp_daily_date %in% holiday_vec),
      .groups = "drop"
    ) %>%
    mutate(Holiday = if_else(has_holiday, 1, -1))

  result <- dat_week %>%
    left_join(week_flags, by = c("Year", "Week"))
  
  return(result)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Data Loading and Pre-processing
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Macau Data
dat.raw.mo <- data.table::fread(file.path(dat.dir, 'Macau/INFLU-U2.csv.gz'), select = c('Age', 'Year', 'Week', 'FLUA', 'FLUB'))
dat.raw.mo <- dat.raw.mo[FLUA != -1 & FLUB != -1]
dat.raw.mo[, Positive := as.integer(FLUA == 1 | FLUB == 1)]

# Harmonized Age Grouping
age.breaks.mo <- c(0, 6, 12, 18, 50, 65, Inf)
age.labels.mo <- c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 50)', '[50, 65)', '65+')
dat.raw.mo[, AgeGroup := cut(Age, breaks = age.breaks.mo, labels = age.labels.mo, right = FALSE)]
dat.week.mo <- dat.raw.mo[, .(Positive_Cases = sum(Positive, na.rm = TRUE), Total_Tests = .N), by = .(Year, Week, AgeGroup)]

dat.cli.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Climate.xlsx')) %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()
dat.air.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Air.Quality.xlsx')) %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()
dat.hol.mo <- readr::read_tsv(file.path(dat.dir, 'Macau/holidays.txt'), show_col_types = FALSE, col_names = FALSE) %>% dplyr::select(date = X1, holiday = X2) %>% as.data.table()

dat.m.mo <- merge(dat.cli.mo, dat.air.mo[, .(Year, Week, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO)], by = c('Year', 'Week'), all.x = TRUE)
dat.m.mo <- merge(dat.m.mo, dat.week.mo, by = c('Year', 'Week'), all.x = TRUE)
dat.m.mo <- dat.m.mo[!is.na(AgeGroup)] %>% na.omit()

dat.m.mo[is.na(PM2.5), PM2.5 := mean(PM2.5, na.rm = TRUE)]
dat.m.mo[is.na(PM10), PM10 := mean(PM10, na.rm = TRUE)]
dat.m.mo[is.na(O3), O3 := mean(O3, na.rm = TRUE)]
dat.m.mo[is.na(SO2), SO2 := mean(SO2, na.rm = TRUE)]
dat.m.mo[is.na(NO2), NO2 := mean(NO2, na.rm = TRUE)]
dat.m.mo[is.na(CO), CO := mean(CO, na.rm = TRUE)]
dat.m.mo[, Season := ifelse(month(date) %in% c(12, 1, 2, 3, 4), -1, 1)]
dat.m.mo$Pandemic = -1
dat.m.mo[dat.m.mo$date > '2020-01-01' & dat.m.mo$date < '2023-01-01', ]$Pandemic = 1
dat.m.mo <- dat.m.mo[Year >= 2010 & Year <= 2025]
dat.m.mo <- dat.m.mo %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.mo <- dat.m.mo %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))
dat.m.mo$Pos_Rate <- dat.m.mo$Positive_Cases / dat.m.mo$Total_Tests
dat.m.mo <- add_holiday_flag(dat.m.mo, dat.hol.mo)

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
colnames(dat.admission.hk.long)[1:2] <- c('Year', 'Week')
dat.admission.hk.long$Year <- as.numeric(dat.admission.hk.long$Year)

dat.cli.hk <- readxl::read_excel(glue('{dat.dir}/HK/Climate.xlsx')) %>% 
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()
dat.air.hk <- readxl::read_excel(glue('{dat.dir}/HK/Air.Quality.xlsx')) %>% 
  dplyr::select(date, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()
dat.cli.hk <- dat.cli.hk[, .(ave.temp = mean(ave.temp, na.rm = TRUE), humidity = mean(humidity, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]
dat.air.hk <- dat.air.hk[, .(O3 = mean(O3, na.rm = TRUE), PM2.5 = mean(PM2.5, na.rm = TRUE), PM10 = mean(PM10, na.rm = TRUE), CO = mean(CO, na.rm = TRUE), SO2 = mean(SO2, na.rm = TRUE), NO2 = mean(NO2, na.rm = TRUE)), by = .(Year, Week)]
dat.hol.hk <- readr::read_tsv(file.path(dat.dir, 'HK/holidays.txt'), show_col_types = FALSE, col_names = FALSE) %>% dplyr::select(date = X1, holiday = X2) %>% as.data.table()

dat.m.hk <- merge(dat.cli.hk, dat.air.hk, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk <- merge(dat.m.hk, dat.admission.hk.long, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk[, Season := ifelse(month(date) %in% c(12, 1, 2, 3, 4), -1, 1)]
dat.m.hk$Pandemic = -1
dat.m.hk[dat.m.hk$date > '2020-01-01' & dat.m.hk$date < '2023-01-01', ]$Pandemic = 1
dat.m.hk <- dat.m.hk[Year >= 2014 & Year <= 2025]
dat.m.hk <- dat.m.hk[!is.na(AgeGroup)]
dat.m.hk <- dat.m.hk %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.hk <- dat.m.hk %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))
dat.m.hk <- add_holiday_flag(dat.m.hk, dat.hol.hk)

# 1.2 England Data
dat.week.en <- readxl::read_excel(file.path(dat.dir, 'England/influenza-age.xlsx'))
dat.week.en <- dat.week.en %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()
dat.week.en <- dat.week.en %>% dplyr::select(Year, Week, Age, Positive_Rate = FLUAB)
dat.week.en <- dat.week.en %>% mutate(
  AgeGroup = case_when(
    Age == "00-04" ~ '[0, 5)',
    Age == "05-14" ~ '[5, 15)',
    Age == "15-44" ~ '[15, 45)',
    Age == "45-64" ~ '[45, 65)',
    Age == "65-79" ~ '[65, 80)',
    Age == "80+" ~ '80+',
    TRUE ~ NA_character_
  )
) %>% dplyr::select(Year, Week, AgeGroup, Rate = Positive_Rate)

dat.cli.en <- readxl::read_excel(glue('{dat.dir}/England/Climate.xlsx')) %>% 
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()
dat.air.en <- readxl::read_excel(glue('{dat.dir}/England/Air.Quality.xlsx')) %>% 
  dplyr::select(date, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()
dat.cli.en <- dat.cli.en[, .(ave.temp = mean(ave.temp, na.rm = TRUE), humidity = mean(humidity, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]
dat.air.en <- dat.air.en[, .(O3 = mean(O3, na.rm = TRUE), PM2.5 = mean(PM2.5, na.rm = TRUE), PM10 = mean(PM10, na.rm = TRUE), CO = mean(CO, na.rm = TRUE), SO2 = mean(SO2, na.rm = TRUE), NO2 = mean(NO2, na.rm = TRUE)), by = .(Year, Week)]

dat.m.en <- merge(dat.cli.en, dat.air.en, by = c('Year', 'Week'), all.x = TRUE)
dat.m.en <- merge(dat.m.en, dat.week.en, by = c('Year', 'Week'), all.x = TRUE)
dat.m.en[, Season := ifelse(month(date) %in% c(10, 11, 12, 1, 2, 3, 4), -1, 1)]
dat.m.en$Pandemic = -1
dat.m.en[dat.m.en$date > '2020-03-23' & dat.m.en$date < '2022-02-24', ]$Pandemic = 1
dat.m.en <- dat.m.en[Year >= 2014 & Year <= 2026]
dat.m.en <- dat.m.en[!is.na(AgeGroup)]
dat.m.en <- dat.m.en %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.en <- dat.m.en %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))

# 1.4 Zhuhai Data
dat.week.zh <- readxl::read_excel(file.path(dat.dir, 'Zhuhai/influenza-age.xlsx'))
dat.week.zh$Rate <- dat.week.zh$Positive_Cases / dat.week.zh$Total_Tests
dat.week.zh$Positive_Cases <- NULL
dat.week.zh$Total_Tests <- NULL
dat.cli.zh <- readxl::read_excel(file.path(dat.dir, 'Zhuhai/Climate.xlsx')) %>% mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% as.data.table()
dat.air.zh <- readxl::read_excel(glue('{dat.dir}/Zhuhai/Air.Quality.xlsx')) %>% 
  dplyr::select(date, PM2.5 = IAQI_PM2.5, PM10 = IAQI_PM10, O3 = IAQI_O3, SO2 = IAQI_SO2, NO2 = IAQI_NO2, CO = IAQI_CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>% 
  as.data.table()

dat.cli.zh <- dat.cli.zh[, .(ave.temp = mean(ave.temp, na.rm = TRUE), humidity = mean(humidity, na.rm = TRUE), pressure = mean(pressure, na.rm = TRUE), date = min(date)), by = .(Year, Week)]
dat.air.zh <- dat.air.zh[, .(O3 = mean(O3, na.rm = TRUE), PM2.5 = mean(PM2.5, na.rm = TRUE), PM10 = mean(PM10, na.rm = TRUE), CO = mean(CO, na.rm = TRUE), SO2 = mean(SO2, na.rm = TRUE), NO2 = mean(NO2, na.rm = TRUE)), by = .(Year, Week)]

dat.m.zh <- merge(dat.cli.zh, dat.air.zh, by = c('Year', 'Week'), all.x = TRUE)
dat.m.zh <- merge(dat.week.zh, dat.m.zh, by = c('Year', 'Week'), all.x = TRUE)
dat.m.zh <- dat.m.zh[!is.na(dat.m.zh$date), ]
dat.m.zh <- dat.m.zh %>% as.data.table()
dat.m.zh[, Season := ifelse(month(date) %in% c(12, 1, 2, 3, 4), -1, 1)]
dat.m.zh$Pandemic = -1
dat.m.zh[dat.m.zh$date > '2020-01-01' & dat.m.zh$date < '2023-01-01', ]$Pandemic = 1
dat.m.zh <- dat.m.zh[!is.na(AgeGroup)]
dat.m.zh <- dat.m.zh %>% mutate(`Primary Pollutants` = pmax(SO2, NO2, CO, na.rm = TRUE))
dat.m.zh <- dat.m.zh %>% mutate(`Particulates` = pmax(PM10, PM2.5, na.rm = TRUE))
dat.m.zh <- dat.m.zh %>% na.omit()

#
dat.m.mo <- dat.m.mo %>% dplyr::select(Age = AgeGroup, Influenza = Pos_Rate, O3, Temperature = ave.temp, Moisture = humidity, Atmosphere = pressure, Particulates, `Primary Pollutants`, Pandemic, Holiday, Season)
dat.m.hk <- dat.m.hk %>% dplyr::select(Age = AgeGroup, Influenza = Rate, O3, Temperature = ave.temp, Moisture = humidity, Atmosphere = pressure, Particulates, `Primary Pollutants`, Pandemic, Holiday, Season)
dat.m.en <- dat.m.en %>% dplyr::select(Age = AgeGroup, Influenza = Rate, O3, Temperature = ave.temp, Moisture = humidity, Atmosphere = pressure, Particulates, `Primary Pollutants`, Pandemic, Season)
dat.m.zh <- dat.m.zh %>% dplyr::select(Age = AgeGroup, Influenza = Rate, O3, Temperature = ave.temp, Moisture = humidity, Atmosphere = pressure, Particulates, `Primary Pollutants`, Pandemic, Season)

# Multigroup Ising Network Modeling
for (place in c('England', 'Zhuhai', 'Macau', 'HK')) {
  if (place == 'HK') {
    df_full <- dat.m.hk
  } else if (place == 'Zhuhai') {
    df_full <- dat.m.zh
  } else if (place == 'England') {
    df_full <- dat.m.en
  } else {
    df_full <- dat.m.mo
  }
  
  # Define node variables to include in the network
  if (place %in% c('Macau', 'HK')) {
    node_vars <- c('Influenza', 'O3', 'Temperature', 'Moisture', 'Atmosphere', 'Particulates', 'Primary Pollutants', 'Pandemic', 'Holiday', 'Season')
    vvvv_vars <- c('Influenza', 'O3', 'Temperature', 'Moisture', 'Atmosphere', 'Particulates', 'Primary Pollutants')
  } else {
    node_vars <- c('Influenza', 'O3', 'Temperature', 'Moisture', 'Atmosphere', 'Particulates', 'Primary Pollutants', 'Pandemic', 'Season')
    vvvv_vars <- c('Influenza', 'O3', 'Temperature', 'Moisture', 'Atmosphere', 'Particulates', 'Primary Pollutants')
  }
  
  # Ising models require complete binary data; filter and prepare dataset
  df_prepared <- df_full %>%
    dplyr::select(Age, all_of(node_vars)) %>%
    na.omit()
  
  # Dichotomize continuous variables using a global median split (High=1, Low=-1)
  df_binary <- df_prepared %>%
    mutate(
      across(
        all_of(vvvv_vars),
        ~ case_when(
          . > median(., na.rm = TRUE) ~ 1,
          . <= median(., na.rm = TRUE) ~ -1,
          TRUE ~ NA_real_
        )
      ))
  
  # Fit a series of nested multigroup Ising models
  force = FALSE
  print(node_vars)
  
  # Model 1: Saturated model (all parameters free to vary by age group)
  start <- Sys.time()
  model1 <- Ising(df_binary, vars = node_vars, groups = "Age") %>% runmodel(approximate_SEs = TRUE)
  print('model1')
  print(Sys.time() - start)
  
  model1b <- model1 %>% prune(alpha = 0.05) %>% stepup(alpha = 0.05)
  print('model1b')
  print(Sys.time() - start)
  
  # Model 2: Constrain network structure (omega) to be equal across age groups
  model2 <- model1 %>% groupequal("omega") %>% runmodel(approximate_SEs = TRUE)
  print('model2')
  print(Sys.time() - start)
  
  model2b <- model2 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
  print('model2b')
  print(Sys.time() - start)
  
  # Model 3: Constrain both network structure (omega) and thresholds (tau) to be equal
  model3 <- model2 %>% groupequal("tau") %>% runmodel(approximate_SEs = TRUE)
  print('model3')
  print(Sys.time() - start)
  
  model3b <- model3 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
  print('model3b')
  print(Sys.time() - start)
  
  # Model 4: Constrain all parameters (omega, tau, beta) to be equal across age groups
  model4 <- model3 %>% groupequal("beta") %>% runmodel(approximate_SEs = TRUE)
  print('model4')
  print(Sys.time() - start)
  
  model4b <- model4 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
  print('model4b')
  print(Sys.time() - start)
  
  print(Sys.time())
  
  # Compare models using Information Criteria to find the optimal balance of fit and parsimony
  model_comparison <- psychonetrics::compare(
    `1. all parameters free (dense)` = model1,
    `2. all parameters free (sparse)` = model1b,
    `3. equal networks (dense)` = model2,
    `4. equal networks (sparse)` = model2b,
    `5. equal networks and thresholds (dense)` = model3,
    `6. equal networks and thresholds (sparse)` = model3b,
    `7. all parameters equal (dense)` = model4,
    `8. all parameters equal (sparse)` = model4b
  ) %>% arrange(BIC, AIC)
  
  print("模型比较结果 (按 BIC, AIC 升序排列):")
  print(model_comparison)
  
  # ===================================================================
  # Result Extraction and Visualization
  # ===================================================================
  # Select the best-fitting model
  best.model <- model3 
  
  # Extract the shared network adjacency matrix and plot using qgraph
  network_matrix_all <- getmatrix(best.model, "omega")
  network_matrix <- getmatrix(best.model, "omega")[[1]]
  
  # extract edge weights and plot the network
  pdf(glue('{fig.dir}/Network_{place}.pdf'), width = 6, height = 6, family = 'Times')
  hh <- qgraph(
    network_matrix, 
    borders = FALSE,
    layout = 'spring', 
    labels = node_vars,
    theme = 'colorblind', 
    label.prop = 0.9, 
    node.width = 1.2, 
    vsize = 18,
    vsize2 = 9,
    label.cex = 1.2,
    label.norm = '000000',
    shape = 'ellipse',
    label.scale = FALSE,
    posCol = '#E50914',
    negCol = '#00A087'
  )
  dev.off()
  
  # Plot Network Temperature across age groups
  # Network temperature (1/beta) reflects system stochasticity: higher values indicate weaker constraints (more randomness)
  beta_estimates <- getmatrix(best.model, "beta")
  if (is.list(beta_estimates)) {
    betas_est <- sapply(beta_estimates, mean)
  } else {
    betas_est <- beta_estimates
  }
  print(betas_est)
  
  # Extract standard errors to compute 95% Confidence Intervals
  beta_params <- best.model@parameters %>% filter(matrix == "beta")
  betas_se <- beta_params$se
  
  # Define age group order based on region
  if (place == 'England') {
    age_groups_ordered <- c("[0, 5)", "[5, 15)", "[15, 45)", "[45, 65)", "[65, 80)", "80+" )
  } else {
    age_groups_ordered <- c("[0, 6)", "[6, 12)", "[12, 18)", "[18, 50)", "[50, 65)", "65+" )
  }
  
  z <- qnorm(0.975)
  upperCI <- betas_est + (z * betas_se)
  lowerCI <- betas_est - (z * betas_se)
  
  # Calculate inverse temperature (Note: CI bounds are inversely mapped)
  temp_data <- data.frame(
    Age = age_groups_ordered,
    Temperature = 1/betas_est,
    UpperCI = 1/lowerCI,
    LowerCI = 1/upperCI
  )
  print(temp_data)
  
  age_order <- age_groups_ordered
  temp_data$Age <- factor(temp_data$Age, levels = age_order)
  
  # Plot
  base_size = 16
  base_color = '#000000'
  base_family = 'serif'
  gg <- ggplot(temp_data, aes(x = Age, y = Temperature, group = 1)) +
    geom_point(shape = 16, size = 3) +
    labs(x = 'Age Group', y = 'Network Temperature') +
    scale_y_continuous(limits = c(0.8, 1.3), expand = c(0, 0)) +
    geom_line() +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
    theme.com 
  
  ggsave(gg, filename = glue('{fig.dir}/{ofig}_{place}.pdf'), width = 6, height = 6, dpi = 300, bg = '#FFFFFF')
}


