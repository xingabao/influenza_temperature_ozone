# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(data.table)))

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
ofig <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting Theme
base.col <- '#000000'
base.size <- 18
base.family <- 'serif'
theme.com <- theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(family = base.family, color = 'black', size = base.size),
    axis.text = element_text(color = 'black', size = base.size * 0.8, family = base.family),
    axis.title = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.text = element_text(size = base.size * 0.8, face = 'bold'),
    legend.position = 'top'
  )

iso_week_date <- function(year, week) {
  suppressWarnings(lubridate::ymd(paste0(year, '-W', sprintf('%02d', week), '-1')))
}

zscore <- function(x) {
  x <- as.numeric(x)
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}

# Macau: clinical + env -> weekly dataset
dat.raw.mo <- data.table::fread(
  file.path(dat.dir, 'Macau/INFLU-U2.csv.gz'),
  select = c('Age', 'Year', 'Week', 'FLUA', 'FLUB')
)
dat.raw.mo <- dat.raw.mo[FLUA != -1 & FLUB != -1]
dat.raw.mo[, Positive := as.integer(FLUA == 1 | FLUB == 1)]

age.breaks.mo <- c(0, 6, 12, 18, 50, 65, Inf)
age.labels.mo <- c('[0, 6)', '[6, 12)', '[12, 18)', '[18, 50)', '[50, 65)', '65+')
dat.raw.mo[, AgeGroup := cut(Age, breaks = age.breaks.mo, labels = age.labels.mo, right = FALSE)]

dat.week.mo <- dat.raw.mo[
  , .(Positive_Cases = sum(Positive, na.rm = TRUE), Total_Tests = .N),
  by = .(Year, Week, AgeGroup)
]

dat.cli.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Climate.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>%
  as.data.table()

dat.air.mo <- readxl::read_excel(file.path(dat.dir, 'Macau/Air.Quality.xlsx')) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>%
  as.data.table()

dat.m.mo <- merge(
  dat.cli.mo,
  dat.air.mo[, .(Year, Week, PM2.5, PM10, O3, SO2, NO2, CO)],
  by = c('Year', 'Week'), all.x = TRUE
)
dat.m.mo <- merge(dat.m.mo, dat.week.mo, by = c('Year', 'Week'), all.x = TRUE)

dat.m.mo <- dat.m.mo[Year >= 2010 & Year <= 2025]
dat.m.mo <- dat.m.mo[!is.na(AgeGroup)]
dat.m.mo <- dat.m.mo %>%
  dplyr::select(
    Age = AgeGroup, nCase = Positive_Cases, nSample = Total_Tests,
    Year, Week, date,
    pressure, max.temp, ave.temp, min.temp, dew.temp,
    humidity, sunshine, wind.speed, precipitation,
    `PM2.5`, PM10, O3, SO2, NO2, CO
  ) %>% as.data.table()

# Hong Kong: admissions + env -> weekly dataset
dat.raw.hk <- readxl::read_excel(file.path(dat.dir, 'HK/flux_data.xlsx'), skip = 2)
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
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>%
  as.data.table()

dat.air.hk <- readxl::read_excel(file.path(dat.dir, 'HK/Air.Quality.xlsx')) %>%
  dplyr::select(date, `PM2.5`, PM10, O3, SO2, NO2, CO) %>%
  mutate(date = as.Date(date), Year = year(date), Week = week(date)) %>%
  as.data.table()

dat.air.hk <- dat.air.hk[
  , .(O3 = mean(O3, na.rm = TRUE),
      `PM2.5` = mean(`PM2.5`, na.rm = TRUE),
      PM10 = mean(PM10, na.rm = TRUE),
      CO = mean(CO, na.rm = TRUE),
      SO2 = mean(SO2, na.rm = TRUE),
      NO2 = mean(NO2, na.rm = TRUE)),
  by = .(Year, Week)
]

dat.m.hk <- merge(dat.cli.hk, dat.air.hk, by = c('Year', 'Week'), all.x = TRUE)
dat.m.hk <- merge(dat.m.hk, dat.admission.hk.long, by = c('Year', 'Week'), all.x = TRUE)

dat.m.hk <- dat.m.hk[Year >= 2014 & Year <= 2025]
dat.m.hk <- dat.m.hk[!is.na(AgeGroup)]
dat.m.hk <- dat.m.hk %>%
  dplyr::select(
    Age = AgeGroup, nCase = Positive_Cases, nSample = Population_Base,
    Year, Week, date,
    pressure, max.temp, ave.temp, min.temp, dew.temp,
    humidity, sunshine, wind.speed, precipitation,
    `PM2.5`, PM10, O3, SO2, NO2, CO
  ) %>% as.data.table()

# Weekly ISO date + paired weeks
dat.m.mo[, wk_date := ISOweek::ISOweek2date(paste(Year, sprintf("W%02d", Week), 1, sep = "-"))]
dat.m.hk[, wk_date := ISOweek::ISOweek2date(paste(Year, sprintf("W%02d", Week), 1, sep = "-"))]
dat.m.mo <- dat.m.mo[!is.na(wk_date)]
dat.m.hk <- dat.m.hk[!is.na(wk_date)]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Macau vs Hong Kong overlay (15 vars, z-score by region)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mo_env <- unique(dat.m.mo[, .(
  wk_date, pressure, max.temp, ave.temp, min.temp, dew.temp, humidity, sunshine, wind.speed,
  precipitation, `PM2.5`, PM10, O3, SO2, NO2, CO
)])
hk_env <- unique(dat.m.hk[, .(
  wk_date, pressure, max.temp, ave.temp, min.temp, dew.temp, humidity, sunshine, wind.speed,
  precipitation, `PM2.5`, PM10, O3, SO2, NO2, CO
)])

dat_pair <- merge(mo_env, hk_env, by = 'wk_date', all = FALSE, suffixes = c('_mo', '_hk'))

vars_map15 <- data.frame(
  facet = c(
    'Atmospheric Pressure', 'Maximum Temperature', 'Average Temperature', 'Minimum Temperature', 'Dew Point Temperature', 'Humidity', 'Sunshine Duration', 'Wind Speed', 'Precipitation',
    'PM2.5', 'PM10', 'O3', 'SO2', 'NO2', 'CO'
  ),
  mo_col = c(
    'pressure_mo', 'max.temp_mo', 'ave.temp_mo', 'min.temp_mo', 'dew.temp_mo', 'humidity_mo',
    'sunshine_mo', 'wind.speed_mo', 'precipitation_mo',
    'PM2.5_mo', 'PM10_mo', 'O3_mo', 'SO2_mo', 'NO2_mo', 'CO_mo'
  ),
  hk_col = c(
    'pressure_hk', 'max.temp_hk', 'ave.temp_hk', 'min.temp_hk', 'dew.temp_hk', 'humidity_hk',
    'sunshine_hk', 'wind.speed_hk', 'precipitation_hk',
    'PM2.5_hk', 'PM10_hk', 'O3_hk', 'SO2_hk', 'NO2_hk', 'CO_hk'
  ),
  stringsAsFactors = FALSE
)

plot_df <- lapply(seq_len(nrow(vars_map15)), function(i) {
  m <- vars_map15[i, ]
  tmp <- dat_pair[, .(wk_date, mo = get(m$mo_col), hk = get(m$hk_col))]
  tmp <- tmp[is.finite(mo) & is.finite(hk)]
  if (nrow(tmp) == 0) return(NULL)
  
  tmp[, mo_z := zscore(mo)]
  tmp[, hk_z := zscore(hk)]
  
  rbind(
    tmp[, .(wk_date, region = 'Macau', value_z = mo_z, facet = m$facet)],
    tmp[, .(wk_date, region = 'Hong Kong', value_z = hk_z, facet = m$facet)]
  )
}) %>% data.table::rbindlist(fill = TRUE)

plot_df$facet <- factor(plot_df$facet, levels = vars_map15$facet)

# Macau vs Hong Kong: Weekly synchrony in meteorology and air pollution
gg. <- ggplot(plot_df, aes(x = wk_date, y = value_z, color = region)) +
  geom_line(linewidth = 0.6, alpha = 0.95) +
  facet_wrap(~ facet, ncol = 3, scales = 'free_y') +
  scale_color_manual(values = c('Macau' = '#D95F02', 'Hong Kong' = '#00A087')) +
  labs(
    x = 'Date (weekly)',
    y = 'Z-score',
    color = NULL,
  ) +
  theme.com +
  theme(
    legend.position = 'none',
    strip.text = element_text(color = base.col, size = base.size * 0.9, family = base.family),
    strip.background = element_rect(linewidth = 0),
    panel.grid = element_blank()
  ); gg.

pals <- c(rep('#00A087', 6), rep('#E50914', 9))
gg <- ggplot_gtable(ggplot_build(gg.))
stript <- which(grepl('strip-t', gg$layout$name))
k <- 1
for (i in stript) {
  j <- which(grepl('rect', gg$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- paste0(substr(pals[k], 1, 7), '7F')
  k <- k + 1
}

ggsave(glue('{fig.dir}/{ofig}_meteo_air.pdf'), plot = gg, width = 11, height = 11)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Cross-Regional Concordance Analysis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Aggregate weekly influenza activity for comparison (use wk_date as weekly index)
dat.ts.mo <- dat.week.mo[, .(nCase    = sum(Positive_Cases, na.rm = TRUE), nSample  = sum(Total_Tests, na.rm = TRUE)), by = .(Year, Week)]
dat.ts.mo[, PositiveRate := nCase / nSample]
dat.ts.mo[, Region := "Macau"]

dat.ts.hk <- dat.admission.hk.long[, .(nCase    = sum(Positive_Cases, na.rm = TRUE), nSample  = sum(Population_Base, na.rm = TRUE)), by = .(Year, Week)]
dat.ts.hk[, PositiveRate := nCase / nSample]
dat.ts.hk[, Region := "Hong Kong"]

dat.act <- rbindlist(list(dat.ts.mo, dat.ts.hk), fill = TRUE) %>% 
  arrange(Year, Week) %>% 
  mutate(
    wk_date = as.Date(ISOweek::ISOweek2date(sprintf("%d-W%02d-1", Year, Week)))
  ) %>% dplyr::select(-Year, -Week)
dat.act <- dat.act[order(wk_date, Region)][, .SD[1], by = .(wk_date, Region)]

dat.act.wide <- dcast(dat.act, wk_date ~ Region, value.var = "PositiveRate", na.rm = TRUE)
dat.act.wide <- dat.act.wide[!is.na(Macau) & !is.na(`Hong Kong`)]

dat.act.wide[, Z_MO := scale(Macau)]
dat.act.wide[, Z_HK := scale(`Hong Kong`)]

cor.val <- cor(dat.act.wide$Z_MO, dat.act.wide$Z_HK, use = "complete.obs")
cor.label <- paste0("r = ", round(cor.val, 3))

# "Temporal Synchrony of Influenza Activity ({cor.label})"
# "Standardized Positive Rate (MO) vs. Admission Rate (HK)"
print(glue("Temporal Synchrony of Influenza Activity ({cor.label})"))
gg.flu <- ggplot(dat.act.wide, aes(x = as.POSIXct(wk_date))) +
  geom_line(aes(y = Z_HK, color = "Hong Kong"), linewidth = 0.8, alpha = 0.85) +
  geom_line(aes(y = Z_MO, color = "Macau"), linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = c("Macau" = "#D95F02", "Hong Kong" = "#00A087"), breaks = c("Macau", "Hong Kong")) +
  scale_x_datetime(
    breaks = seq(as.POSIXct("2010-01-01"), as.POSIXct("2026-01-01"), by = "1 year"),
    labels = format(seq(as.POSIXct("2010-01-01"), as.POSIXct("2026-01-01"), by = "1 year"), "%Y"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Date (weekly)",
    y = "Standardized Activity Level (Z-score)",
    color = NULL
  ) +
  theme.com +
  theme(
    legend.position = c(0.12, 0.90),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(color = 'black', size = base.size * 0.8, family = base.family)
  ); gg.flu

ggsave(glue('{fig.dir}/{ofig}_MO_vs_HK.pdf'), plot = gg.flu, width = 10, height = 4.5)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Aggregate weekly environmental variables (15 vars) and build wide table for concordance
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
vars_met <- c("pressure", "max.temp", "ave.temp", "min.temp", "dew.temp", "humidity", "sunshine", "wind.speed", "precipitation")
vars_pol <- c("PM2.5", "PM10", "O3", "SO2", "NO2", "CO")
env_vars15 <- c(vars_met, vars_pol)

facet_map <- c(
  pressure       = "Atmospheric Pressure",
  max.temp       = "Maximum Temperature",
  ave.temp       = "Average Temperature",
  min.temp       = "Minimum Temperature",
  dew.temp       = "Dew Point Temperature",
  humidity       = "Humidity",
  sunshine       = "Sunshine Duration",
  wind.speed     = "Wind Speed",
  precipitation  = "Precipitation",
  `PM2.5`        = "PM2.5",
  PM10           = "PM10",
  O3             = "O3",
  SO2            = "SO2",
  NO2            = "NO2",
  CO             = "CO"
)
facet_levels <- unname(facet_map[c(vars_met, vars_pol)])

dat.m.mo. <- dat.m.mo %>%
  arrange(Year, Week) %>%
  mutate(wk_date = as.Date(ISOweek::ISOweek2date(sprintf("%d-W%02d-1", Year, Week)))) %>%
  dplyr::select(-Year, -Week)

dat.m.hk. <- dat.m.hk %>%
  arrange(Year, Week) %>%
  mutate(wk_date = as.Date(ISOweek::ISOweek2date(sprintf("%d-W%02d-1", Year, Week)))) %>%
  dplyr::select(-Year, -Week)

dat.env.mo <- dat.m.mo.[, lapply(.SD, mean, na.rm = TRUE), by = .(wk_date), .SDcols = env_vars15]
setnames(dat.env.mo, env_vars15, paste0(env_vars15, "_Macau"))

dat.env.hk <- dat.m.hk.[, lapply(.SD, mean, na.rm = TRUE), by = .(wk_date), .SDcols = env_vars15]
setnames(dat.env.hk, env_vars15, paste0(env_vars15, "_Hong Kong"))

dat.wide <- merge(dat.env.mo, dat.env.hk, by = "wk_date", all = FALSE)

make_env_plot_one <- function(dat_wide, v, base.family = NULL) {
  xcol <- paste0(v, "_Hong Kong")
  ycol <- paste0(v, "_Macau")
  if (!all(c(xcol, ycol) %in% names(dat_wide))) return(NULL)
  
  dt <- dat_wide[, .(wk_date, x = get(xcol), y = get(ycol))]
  dt <- dt[is.finite(x) & is.finite(y)]
  if (nrow(dt) < 3) return(NULL)
  
  ct <- suppressWarnings(stats::cor.test(dt$x, dt$y, method = "pearson"))
  r_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  r_str <- sprintf("%.3f", r_val)
  v_lab <- if (!is.null(facet_map) && !is.na(facet_map[[v]])) facet_map[[v]] else v
  
  sub_lab <- if (is.na(p_val)) {
    sprintf("%s:<br>R = %s; <i>P</i> = NA", v_lab, r_str)
  } else if (p_val < 0.001) {
    sprintf("%s:<br>R = %s; <i>P</i> < 0.001", v_lab, r_str)
  } else {
    sprintf("%s:<br>R = %s; <i>P</i> = %.3f", v_lab, r_str, p_val)
  }
  
  ggplot(dt, aes(x = x, y = y)) +
    geom_point(alpha = 0.25, color = "grey30", size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, color = "#0072B2", linewidth = 0.7) +
    labs(
      subtitle = sub_lab,
      x = NULL,   # 关键：去掉每个小图的轴标题
      y = NULL
    ) +
    scale_x_continuous(labels = scales::label_number()) +
    scale_y_continuous(labels = scales::label_number()) +
    theme.com +
    theme(
      legend.position = "none",
      plot.subtitle = ggtext::element_markdown(family = base.family, size = base.size * 0.8, color = base.col, lineheight = 1.15),
      text = element_text(size = base.size * 0.8, color = base.col, family = base.family),
      axis.title = element_blank()
    )
}

plots_met <- Filter(Negate(is.null), lapply(vars_met, \(v) make_env_plot_one(dat.wide, v, base.family)))
plots_pol <- Filter(Negate(is.null), lapply(vars_pol, \(v) make_env_plot_one(dat.wide, v, base.family)))

# Cross-Regional Environmental Concordance (Meteorology, 9 Variables)
# Each panel: Macau vs Hong Kong weekly values; line = linear fit; subtitle = Pearson r and P
p.met9 <- wrap_plots(plots_met, ncol = 3) +
  plot_annotation(
    theme = theme(
      axis.title.x = element_text(family = base.family, margin = margin(t = 8)),
      axis.title.y = element_text(family = base.family, margin = margin(r = 8))
    )
  ) &
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text()
  )

p.met9 <- p.met9 + labs(x = "Hong Kong (weekly)", y = "Macau (weekly)")

# Cross-Regional Environmental Concordance (Air Pollutants, 6 Variables)
# Each panel: Macau vs Hong Kong weekly values; line = linear fit; subtitle = Pearson r and P
p.pol6 <- wrap_plots(plots_pol, ncol = 3) +
  plot_annotation(
    theme = theme(
      axis.title.x = element_text(family = base.family, margin = margin(t = 8)),
      axis.title.y = element_text(family = base.family, margin = margin(r = 8))
    )
  ) &
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text()
  )

p.pol6 <- p.pol6 + labs(x = "Hong Kong (weekly)", y = "Macau (weekly)")

ggsave(glue('{fig.dir}/{ofig}_meteorology.pdf'), plot = p.met9, width = 10.5, height = 10.5)
ggsave(glue('{fig.dir}/{ofig}_air_pollution.pdf'), plot = p.pol6, width = 10.5, height = 7)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table S: Summary of Meteorological Variables in Macau and Hong Kong
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
desc_stats <- function(data, vars, region_name) {
  res <- lapply(vars, function(v) {
    x <- data[[v]]
    data.frame(
      Region = region_name,
      Variable = v,
      Mean = mean(x, na.rm = TRUE),
      SD = sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  })
  rbindlist(res)
}

met_vars <- c("pressure", "max.temp", "ave.temp", "min.temp", "dew.temp", "humidity", "sunshine", "wind.speed", "precipitation")
pol_vars <- c("PM2.5", "PM10", "O3", "SO2", "NO2", "CO")

# Macau
stats_mo_met <- desc_stats(dat.m.mo, met_vars, "Macau")
stats_mo_pol <- desc_stats(dat.m.mo, pol_vars, "Macau")
stats_mo_met$Group <- 'Meteorology'
stats_mo_pol$Group <- 'Air Pollution'

# Hong Kong
stats_hk_met <- desc_stats(dat.m.hk, met_vars, "Hong Kong")
stats_hk_pol <- desc_stats(dat.m.hk, pol_vars, "Hong Kong")
stats_hk_met$Group <- 'Meteorology'
stats_hk_pol$Group <- 'Air Pollution'

# Merge
Table_S <- rbind(stats_mo_met, stats_hk_met, stats_mo_pol, stats_hk_pol) %>%
  mutate(Variable = case_when(
    Variable == "ave.temp"   ~ "Average Temperature",
    Variable == "dew.temp"   ~ "Dew Point Temperature",
    Variable == "max.temp"   ~ "Maximum Temperature",
    Variable == "min.temp"   ~ "Minimum Temperature",
    Variable == "wind.speed" ~ "Wind Speed",
    Variable == "pressure" ~ "Atmospheric Pressure",
    Variable == "humidity" ~ "Humidity",
    Variable == "sunshine" ~ "Sunshine Duration",
    Variable == "precipitation" ~ "Precipitation",
    TRUE ~ Variable 
  )) %>%
  mutate(Region = factor(Region, levels = c('Macau', 'Hong Kong'))) %>%
  arrange(Region, Group, Variable) %>%
  dplyr::select(Region, Group, Variable, all_of(everything()))

writexl::write_xlsx(Table_S, glue('{tbl.dir}/Table_S_Env_Descriptive.xlsx'))

# Table S: Cross-Regional Environmental Correlation Matrix
cor_list <- list()
for (v in c(met_vars, pol_vars)) {
  col_mo <- paste0(v, "_mo")
  col_hk <- paste0(v, "_hk")
  if (col_mo %in% names(dat_pair) && col_hk %in% names(dat_pair)) {
    val <- cor(dat_pair[[col_mo]], dat_pair[[col_hk]], use = "complete.obs")
    cor_list[[v]] <- data.frame(Variable = v, Pearson_R = round(val, 4))
  }
}
Table_S <- rbindlist(cor_list) %>%
  mutate(Variable = case_when(
    Variable == "ave.temp"   ~ "Average Temperature",
    Variable == "dew.temp"   ~ "Dew Point Temperature",
    Variable == "max.temp"   ~ "Maximum Temperature",
    Variable == "min.temp"   ~ "Minimum Temperature",
    Variable == "wind.speed" ~ "Wind Speed",
    Variable == "pressure" ~ "Pressure",
    Variable == "humidity" ~ "Humidity",
    Variable == "sunshine" ~ "Sunshine",
    Variable == "precipitation" ~ "Precipitation",
    TRUE ~ Variable 
  )) 
writexl::write_xlsx(Table_S, glue('{tbl.dir}/Table_S_Regional_Synchrony.xlsx'))
