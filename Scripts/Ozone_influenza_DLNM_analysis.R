# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(dlnm)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(splines)))
suppressMessages(suppressWarnings(library(ggplot2)))
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Load data
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
place <- 'Macau' # Macau or HK
dat <- readxl::read_xlsx(glue('{dat.dir}/{place}/FLU-CL-AQ.xlsx'))
if (place == 'HK') {
  dat <- dat %>% arrange(date) %>% filter(date < '2020-03-01')
}
n_years <- as.numeric(difftime(max(dat$date), min(dat$date), units = "days")) / 365.25

o3_median <- median(dat$O3, na.rm = TRUE)
o3_p5     <- quantile(dat$O3, 0.05, na.rm = TRUE)
o3_p95    <- quantile(dat$O3, 0.95, na.rm = TRUE)

# Extract lag-specific Relative Risks (RR) at a specific exposure percentile
extract_lag_rr <- function(pred, target_val, label) {
  idx <- which.min(abs(as.numeric(rownames(pred$matRRfit)) - target_val))
  ml  <- ncol(pred$matRRfit) - 1
  data.frame(
    Lag = 0:ml, 
    RR = pred$matRRfit[idx, ],
    Lower = pred$matRRlow[idx, ], 
    Upper = pred$matRRhigh[idx, ],
    Exposure = label, row.names = NULL
  )
}

# Extract cumulative RR across all lags at a specific exposure percentile
get_cumul_rr <- function(pred, target_val) {
  idx <- which.min(abs(as.numeric(names(pred$allRRfit)) - target_val))
  c(RR = unname(pred$allRRfit[idx]),
    Lo = unname(pred$allRRlow[idx]),
    Hi = unname(pred$allRRhigh[idx]))
}

# Combine lag-specific and cumulative results into structured lists
collect_results <- function(pred, model_label) {
  df_lag <- bind_rows(
    extract_lag_rr(pred, o3_p5, lab_p5),
    extract_lag_rr(pred, o3_p95, lab_p95)
  ) %>% mutate(Model = model_label)
  
  rr5  <- get_cumul_rr(pred, o3_p5)
  rr95 <- get_cumul_rr(pred, o3_p95)
  df_cum <- data.frame(
    Model = model_label,
    Exposure = c("P5", "P95"),
    O3 = c(round(o3_p5, 1), round(o3_p95, 1)),
    cumRR = c(rr5[1], rr95[1]),
    CI_lo = c(rr5[2], rr95[2]),
    CI_hi = c(rr5[3], rr95[3])
  )
  list(lag = df_lag, cumul = df_cum)
}

# Define formatted labels for the 5th and 95th percentiles
lab_p5  <- paste0("P5 (", round(o3_p5, 1), " \u00b5g/m\u00b3)")
lab_p95 <- paste0("P95 (", round(o3_p95, 1), " \u00b5g/m\u00b3)")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Main Model
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cb_o3   <- crossbasis(dat$O3, lag = 12, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 4))
cb_temp <- crossbasis(dat$ave.temp, lag = 12, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 4))

model_main <- glm(
  FLUA ~ cb_o3 + cb_temp + ns(humidity, df = 3) +
    ns(date, df = round(7 * n_years)) + as.factor(holiday),
  family = quasipoisson(), data = dat
)
cat("Overdispersion parameter:", round(summary(model_main)$dispersion, 2), "\n")

pred_main <- crosspred(cb_o3, model_main, cen = o3_median)
res_main  <- collect_results(pred_main, "Main (Lag 12, df=7/yr)")

# Generate fine-grid 3D surface data for the exposure-lag-response heatmap
pred_surf <- crosspred(
  cb_o3, model_main, cen = o3_median,
  at = seq(quantile(dat$O3, 0.01, na.rm = TRUE), quantile(dat$O3, 0.99, na.rm = TRUE), length.out = 50)
)

surf_df <- expand.grid(
  O3  = as.numeric(rownames(pred_surf$matRRfit)),
  Lag = 0:12
)
surf_df$logRR <- log(as.vector(pred_surf$matRRfit))
bnd <- quantile(abs(surf_df$logRR), 0.98, na.rm = TRUE)
surf_df$logRR_c <- pmin(pmax(surf_df$logRR, -bnd), bnd)

# Extract data for the 2D cumulative exposure-response curve
df_cumul <- data.frame(
  O3    = as.numeric(names(pred_main$allRRfit)),
  RR    = pred_main$allRRfit,
  Lower = pred_main$allRRlow,
  Upper = pred_main$allRRhigh,
  row.names = NULL
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Sensitivity analysis
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sensitivity 1: Varying maximum lag structures (4 to 16 weeks)
run_sens_lag <- function(ml) {
  cb1 <- crossbasis(dat$O3, lag = ml, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 4))
  cb2 <- crossbasis(dat$ave.temp, lag = ml, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 4))
  m <- glm(FLUA ~ cb1 + cb2 + ns(humidity, df = 3) + ns(date, df = round(7 * n_years)) + as.factor(holiday), family = quasipoisson(), data = dat)
  p <- crosspred(cb1, m, cen = o3_median)
  collect_results(p, paste0("Max Lag ", ml))
}

# Sensitivity 2: Additional adjustment for Sunshine duration
m_ss <- glm(
  FLUA ~ cb_o3 + cb_temp + ns(sunshine, df = 3) + ns(humidity, df = 3) +
     ns(date, df = round(7 * n_years)) + as.factor(holiday),
  family = quasipoisson(), data = dat
)
p_ss   <- crosspred(cb_o3, m_ss, cen = o3_median)
res_ss <- collect_results(p_ss, "Adj. Sunshine")

# Sensitivity 3: Two-pollutant model adjusting for PM2.5
m_pm <- glm(
  FLUA ~ cb_o3 + cb_temp + ns(PM2.5, df = 3) + ns(humidity, df = 3) +
    ns(sunshine, df = 3) + ns(date, df = round(7 * n_years)) + as.factor(holiday),
  family = quasipoisson(), data = dat
)
p_pm   <- crosspred(cb_o3, m_pm, cen = o3_median)
res_pm <- collect_results(p_pm, "Adj. PM2.5")

# Execute lag sensitivity analyses
res_lag4  <- run_sens_lag(4)
res_lag6  <- run_sens_lag(6)
res_lag8  <- run_sens_lag(8)
res_lag10  <- run_sens_lag(10)
res_lag12 <- run_sens_lag(12)
res_lag14 <- run_sens_lag(14)
res_lag16 <- run_sens_lag(16)

# Combine lag-specific sensitivity data
df_sens_lag <- bind_rows(
  res_lag4$lag,
  res_lag6$lag,
  res_lag8$lag,
  res_lag10$lag,
  res_lag12$lag,
  res_lag14$lag,
  res_lag16$lag,
  res_ss$lag,
  res_pm$lag
)
df_sens_lag$Model <- factor(df_sens_lag$Model, levels = c(
  "Max Lag 4",
  "Max Lag 6",
  "Max Lag 8",
  "Max Lag 10",
  "Max Lag 12",
  "Max Lag 14",
  "Max Lag 16",
  "Adj. Sunshine",
  "Adj. PM2.5"
))

# Combine and export cumulative RR summary table to Excel
tab_cumul <- bind_rows(
  res_main$cumul,
  res_lag4$cumul, 
  res_lag6$cumul, 
  res_lag8$cumul,
  res_lag10$cumul,
  res_lag12$cumul, 
  res_lag14$cumul, 
  res_lag16$cumul,
  res_ss$cumul, 
  res_pm$cumul
)
tab_cumul$Result <- sprintf("%.3f (%.3f\u2013%.3f)", tab_cumul$cumRR, tab_cumul$CI_lo, tab_cumul$CI_hi)

cat("\n========== Cumulative RR Summary (vs Median) ==========\n")
writexl::write_xlsx(tab_cumul %>% select(Model, Exposure, O3, Result), glue('{tbl.dir}/Table_S_{ofig}_{place}.xlsx'))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# FIGURE 1: Main Results (3-Panel Plot for Manuscript) # 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
theme.com <- theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(family = base.family, color = 'black', size = base.size),
    axis.text = element_text(color = 'black', size = base.size * 0.8, family = base.family),
    axis.title.x = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold', margin = margin(t = 0.15, unit = 'in')),
    axis.title.y = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.text = element_text(size = base.size * 0.9, face = 'bold'),
    legend.title = element_text(size = base.size * 0.9, face = 'bold'),
    legend.text = element_text(size = base.size * 0.9),
    plot.subtitle = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold')
  )

# Panel A: 3D Exposure-Lag-Response Surface Heatmap
p_A <- ggplot(surf_df, aes(x = Lag, y = O3, fill = logRR_c)) +
  geom_tile() +
  geom_contour(aes(z = logRR), breaks = 0, color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "log(RR)",
    guide = guide_colorbar(barwidth = 0.8, barheight = 6)
  ) +
  geom_hline(yintercept = c(o3_p5, o3_p95), linetype = "dotted", color = "grey20", linewidth = 0.4) +
  annotate("text", x = 0.3, y = o3_p5,  label = "P5", size = base.size / 3.88, vjust = -0.5, fontface = "bold", family = base.family) +
  annotate("text", x = 0.3, y = o3_p95, label = "P95", size = base.size / 3.88, vjust = -0.5, fontface = "bold", family = base.family) +
  scale_x_continuous(breaks = seq(0, 12, 3), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    subtitle = "Exposure\u2013Lag\u2013Response Surface",
    x = "Lag (Weeks)",
    y = expression(bold(O[3] ~ (µg / m^3)))
  ) +
  theme.com + theme(
    legend.position = "inside",
    legend.background = element_rect(fill = '#FFFFFF44'),
    legend.position.inside = c(ifelse(place == 'Macau', 0.8, 0.3), 0.7)
  )

# Panel B: Overall Cumulative Exposure-Response Curve
p_B <- ggplot(df_cumul, aes(x = O3, y = RR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#E50914", alpha = 0.15) +
  geom_line(color = "#E50914", linewidth = 1) +
  geom_vline(xintercept = o3_p5, linetype = "dotted", color = "#1B9E77") +
  geom_vline(xintercept = o3_median, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = o3_p95, linetype = "dotted", color = "#7570B3") +
  scale_x_continuous(limits = c(0, ifelse(place == 'Macau', 250, 120)), expand = c(0, 0)) +
  annotate("text", x = o3_p5, y = Inf, label = "P5", vjust = 1.5, color = "#1B9E77", size = base.size / 3.88, fontface = "bold", family = base.family) +
  annotate("text", x = o3_median, y = Inf, label = "Median", vjust = 1.5, color = "grey40",  size = base.size / 3.88, fontface = "bold", family = base.family) +
  annotate("text", x = o3_p95, y = Inf, label = "P95", vjust = 1.5, color = "#7570B3", size = base.size / 3.88, fontface = "bold", family = base.family) +
  labs(
    subtitle = "Overall Cumulative Exposure\u2013Response",
    x = expression(bold(O[3] ~ (µg / m^3))),
    y = "Cumulative RR (Lag 0\u201312)"
  ) +
  coord_cartesian(clip = "off") +
  theme.com

# Panel C: Lag-Specific RR at the 5th and 95th percentiles
p_C <- ggplot(res_main$lag, aes(x = Lag, y = RR, color = Exposure, fill = Exposure)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("#E50914", "#00A087")) +
  scale_fill_manual(values = c("#E50914", "#00A087")) +
  scale_x_continuous(breaks = 0:12) +
  labs(
    subtitle = "Lag-specific Relative Risk",
    x = "Lag (Weeks)", 
    y = "Relative Risk (RR)"
  ) +
  theme.com + theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.23),
    legend.key.spacing.y = unit(0.25, units = 'cm'),
    legend.background = element_blank()
  )

# Combine panels into final main figure layout
fig1 <- (p_A | p_B) / p_C +
  plot_layout(heights = c(1, 0.85)) 

print(fig1)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# FIGURE 2: Sensitivity Analyses (Supplementary Material) # 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# subtitle = paste0("Reference: Median O\u2083 (", round(o3_median, 1), " µg/m\u00b3)"),
fig2. <- ggplot(df_sens_lag, aes(x = Lag, y = RR, color = Exposure, fill = Exposure)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ Model, ncol = 3, scales = "free_x") +
  scale_color_manual(values = c("#E50914", "#00A087")) +
  scale_fill_manual(values = c("#E50914", "#00A087")) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(
    x = "Lag (Weeks)", y = "Relative Risk (RR)"
  ) +
  theme.com +
  theme(
    strip.text = element_text(color = 'black', size = base.size * 0.9, family = base.family, face = 'bold'),
    strip.background = element_rect(linewidth = 0),
    legend.position = "bottom",
    legend.margin = margin(t = -0.3, unit = 'cm')
  )

fig2 <- ggplot_gtable(ggplot_build(fig2.))
stript <- which(grepl('strip-t', fig2$layout$name))
k <- 1
ppcols <- c(rep('#FF7F00', 1), rep('#0530AD', 2), rep('#FF7F00', 6))
for (i in stript) {
  j <- which(grepl('rect', fig2$grobs[[i]]$grobs[[1]]$childrenOrder))
  fig2$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <- paste0(substr(ppcols[k], 1, 7), '4F')
  k <- k + 1
}

plot(fig2)

# Save Plot
ggsave(fig1, filename = glue('{fig.dir}/{ofig}_{place}.pdf'), width = 10, height = 8.8, bg = NULL)
ggsave(fig2, filename = glue('{fig.dir}/{ofig}_sensitivity_{place}.pdf'), width = 10, height = 10.2, bg = NULL)

