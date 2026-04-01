# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(gtable)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(geomtextpath)))

# Set Env
rt.dir <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig <- tools::file_path_sans_ext(basename(basename(this.path::this.path())))
base.size <- 18
base.family <- 'serif'
base.col <- '#000000'

# Load data
dat.raw <- data.table::fread(glue('{dat.dir}/Macau/INFLU-U2.csv.gz'))

# Prepare plot data
IDF <- dat.raw %>% dplyr::filter(FLUA != -1 & FLUB != -1) %>%
  dplyr::filter(KID >= as.Date('2010-01-01')) %>%
  dplyr::select(Year, Gender, FLUA, FLUB) %>%
  dplyr::filter(Gender %in% c('F', 'M')) %>% 
  mutate(Group = case_when(
    FLUA == 1 ~ 'Influenza A',
    FLUB == 1 ~ 'Influenza B',
    TRUE ~ 'Non-Influenza'
  )) %>%
  count(Year, Group) %>%
  rename(Count = n)

dat <- IDF %>%
  mutate(Group = factor(Group, levels = c('Influenza A', 'Influenza B', 'Non-Influenza')))

lev <- c('Influenza A', 'Influenza B', 'Non-Influenza')

# Set palette
group_colors <- c('Influenza A' = '#E50914', 'Influenza B' = '#00A087', 'Non-Influenza' = '#AAAAAA')

# Construct label text: (xx.x% / xx.x% / xx.x%)
dat.labels <- dat %>%
  group_by(Year) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(
    Pct = Count / Total * 100,
    Label_Part = glue("({sprintf('%.1f', Pct)} %)")
  ) %>%
  arrange(Year, factor(Group, levels = c('Influenza A', 'Influenza B', 'Non-Influenza'))) %>%
  group_by(Year) %>%
  summarise(
    Final_Label = paste(Label_Part, collapse = " / "),
    Max_Val = sum(Count),
    .groups = 'drop'
  )

# Construct curve paths
global_max <- max(dat.labels$Max_Val) * 1.2
buffer_size <- global_max * 0.35 

dat.text.path <- dat.labels %>%
  group_by(Year) %>%
  reframe(
    Final_Label = Final_Label,
    y_pos = c(Max_Val + 2000, Max_Val + buffer_size), 
    type = c("start", "end")
  ) %>%
  ungroup()

dat.text.path[dat.text.path$y_pos > 30000 & dat.text.path$type == 'start' & dat.text.path$Year > 2014, ]$y_pos <- 32000

# Draw Plot
gg <- dat %>%
  mutate(Group = factor(Group, levels = lev), Count = replace_na(Count, 0)) %>%
  arrange(Year, Group) %>%
  group_by(Year) %>%
  mutate(
    ymax = cumsum(Count),
    ymin = lag(ymax, default = 0)
  ) %>%
  mutate(val = cumsum(Count)) %>%
  ungroup() %>%
  ggplot() +
  annotate(
    geom = 'segment',
    x = c(2009.5, 2009.5, 2009.5, 2015, 2015, 2015, 2015, 2015, 2014.5, 2014.5, 2014.5, 2016.5, 2018.5, 2018.5, 2018.5),
    xend = c(2026.5, 2025.5, 2026, 2026, 2026, 2026, 2026, 2026, 2026, 2026, 2025.5, 2025.5, 2025.5, 2025.5, 2025.5),
    y = seq(0, 210000, 15000),
    yend = seq(0, 210000, 15000),
    size = .3, color = 'black'
  ) +
  geom_rect(aes(xmin = Year - 0.35, xmax = Year + 0.35, ymin = ymin, ymax = ymax, fill = fct_rev(Group)), color = NA) +
  geom_text(data = dat, aes(label = Year, x = Year, y = 0), size = base.size / 3.88, angle = 90, hjust = 1.2, family = base.family) +
  geom_text(data = dat, aes(label = '|', x = Year, y = 0), vjust = 1.5, size = 1.0, family = base.family) +
  scale_fill_manual(values = group_colors) +
  coord_polar(theta = 'y', clip = 'off', start = 4.71) +
  scale_x_continuous(expand = expansion(mult = c(0.10, -0.10))) +
  scale_y_continuous(
    limits = c(0, NA), 
    expand = expansion(mult = c(0, 0.05)), 
    breaks = seq(0, 210000, 15000), 
    labels = scales::comma_format()
  ) +
  geom_textpath(
    data = dat.text.path,
    aes(x = Year, y = y_pos, label = Final_Label, group = Year),
    vjust = 0.5,
    hjust = 0,
    size = base.size / 3.88,
    text_only = TRUE,
    family = base.family,
    color = base.col
  ) +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    text = element_text(),
    axis.text.x = element_text(color = base.col, size = base.size * 0.8, family = base.family),
    legend.position = c(.82, .02),
    legend.title = element_blank(),
    legend.text = element_text(color = base.col, size = base.size * 0.8, family = base.family),
    legend.key.spacing.y = unit(.25, units = 'cm'),
    plot.margin = margin(t = .05, r = .05, b = .05, l = .05, unit = 'cm')
  ) +
  guides(fill = guide_legend(reverse = TRUE))

print(gg)

# Three-Line Table Construction
dat.tbl <- dat %>%
  group_by(Year) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(
    Cell_Label = glue("{Count} ({sprintf('%.1f', Count/Total*100)}%)")
  ) %>%
  dplyr::select(Year, Group, Cell_Label) %>%
  pivot_wider(names_from = Group, values_from = Cell_Label) %>%
  dplyr::select(Year, `Influenza A`, `Influenza B`, `Non-Influenza`)

header_colors <- c('black', group_colors['Influenza A'], group_colors['Influenza B'], '#000000')

nr <- nrow(dat.tbl)
nc <- ncol(dat.tbl)

body_colors_matrix <- matrix(NA, nrow = nr, ncol = nc)
body_colors_matrix[, 1] <- "black"
body_colors_matrix[, 2] <- group_colors['Influenza A']
body_colors_matrix[, 3] <- group_colors['Influenza B']
body_colors_matrix[, 4] <- '#000000'

tbl_theme <- ttheme_minimal(
  base_size = base.size * 0.8, 
  base_family = base.family,
  core = list(
    fg_params = list(col = body_colors_matrix, hjust = 1, x = 0.95), 
    bg_params = list(fill = "white", col = NA),
    padding = unit(c(4, 4), "mm")
  ),
  colhead = list(
    fg_params = list(col = header_colors, fontface = "bold", hjust = 1, x = 0.95),
    bg_params = list(fill = "white", col = NA),
    padding = unit(c(4, 4), "mm")
  )
)

tbl_grob <- tableGrob(dat.tbl, rows = NULL, theme = tbl_theme)
tbl_grob$widths <- unit(c(1.5, 4, 4, 4), "cm") 

tbl_grob <- gtable_add_grob(
  tbl_grob,
  grobs = segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(lwd = 2, col = "black")),
  t = 1, b = 1, l = 1, r = nc
)

tbl_grob <- gtable_add_grob(
  tbl_grob,
  grobs = segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(lwd = 1, col = "black")),
  t = 1, b = 1, l = 1, r = nc
)

tbl_grob <- gtable_add_grob(
  tbl_grob,
  grobs = segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(lwd = 2, col = "black")),
  t = nr + 1, b = nr + 1, l = 1, r = nc
)

final_plot <- ggdraw() +
  draw_plot(gg, x = 0, y = 0.05) +
  draw_grob(
    tbl_grob, 
    x = -0.18,
    y = -0.25,
    hjust = 0, 
    vjust = 0, 
    scale = 1.2 
  )

print(final_plot)

# Save to file
width = 10; height = 10.5
ggsave(final_plot, filename = glue('{fig.dir}/{ofig}.pdf'), width = width, height = height, bg = '#FFFFFF')