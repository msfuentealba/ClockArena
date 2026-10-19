#!/usr/bin/env Rscript
# Run from the project root: Rscript paper/scripts/build_figure1.R
# Generates Figure 1 (PDF and PNG), Supplementary Figure 1, and Supplementary Tables 1–2.
colors <- list(
  tone = c("#f6f2ed","#dedbca","#c4c0a5","#a59f83","#878264","#5e5948"),
  grey = c("#e5e4e8","#c5c9d6","#959fb3","#6e788c","#425468","#1b2942"),
  olive = c("#f2edb2","#dbdc63","#c4c400","#95a008","#637314","#304215"),
  green = c("#d7e5c5","#9fc978","#5db342","#41912f","#1d6e29","#0e3716"),
  teal = c("#c9e4ef","#96ced3","#48bcbc","#00959f","#006479","#003547"),
  blue = c("#c5e4fb","#9bc9e8","#5495ce","#006eae","#01478c","#002259"),
  purple = c("#e9d3e7","#d1a9ce","#b678b3","#a44990","#792373","#430b4d"),
  red = c("#f5cfc9","#e9a0a4","#db6463","#c5373c","#9c241b","#730c0d"),
  orange = c("#fbdcbc","#f9bd7b","#f29741","#e96700","#b34a00","#832a00"),
  yellow = c("#ffedc1","#f7dc86","#e8c54d","#ca9a23","#9b730a","#685409"),
  skin = c("#f6e5d3","#dcbc9f","#bc9778","#906852","#734d3d","#422a17"))

# Figure 1
# Run from the Metadata project root.
# Inputs: local clock table, CpG membership matrix, and integrated metadata.
# Packages: tidyverse, cowplot, scales, readxl, writexl, maps, ggpubr, ragg.
library(tidyverse)
library(cowplot)
library(scales)
.libPaths(c(".R-library", .libPaths()))
library(ggpubr)

# Setup -----------------------------------------------------------------------
paper <- 'paper'
data_dir <- file.path(paper, 'input')
out_dir <- file.path(paper, 'output')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
palette_3 <- vapply(colors, function(x) x[3], character(1))
teal <- palette_3[['teal']]; blue <- palette_3[['blue']]; coral <- palette_3[['red']]
tissue_colors <- c('Blood' = palette_3[['grey']], 'Saliva' = coral, 'Multiple tissues' = blue,
                   'Other tissues' = palette_3[['yellow']], 'Unspecified' = teal)
theme_set(theme_pubr() + theme(legend.title = element_blank()))
# Prevent intermediate plots from creating Rplots.pdf in batch runs.
pdf(NULL)
# Diagrams and text tables retain theme_pubr(), with non-data axes hidden.
no_axes <- theme(axis.title = element_blank(), axis.text = element_blank(),
                 axis.ticks = element_blank(), axis.line = element_blank())

# A. Vector workflow: literature and methylation data converge at integration ---
# Rounded cards and editable scientific icons are drawn with base R grid.
workflow_icon <- function(type, color) {
  pen <- grid::gpar(col = color, fill = 'white', lwd = 1.5)
  line <- function(x, y) grid::linesGrob(x, y, gp = pen)
  circle <- function(x, y, r, fill = 'white')
    grid::circleGrob(x, y, r, gp = grid::gpar(col = color, fill = fill, lwd = 1.3))
  rectangle <- function(x, y, w, h, fill = 'white')
    grid::rectGrob(x, y, w, h, gp = grid::gpar(col = color, fill = fill, lwd = 1.2))
  parts <- switch(type,
    paper = list(rectangle(.5, .5, .65, .85),
                 line(c(.3, .7), c(.7, .7)), line(c(.3, .7), c(.5, .5)),
                 line(c(.3, .6), c(.3, .3))),
    curate = list(rectangle(.5, .5, .65, .85),
                  line(c(.27, .37, .52), c(.67, .57, .77)),
                  line(c(.27, .37, .52), c(.3, .2, .4)),
                  line(c(.6, .75), c(.67, .67)), line(c(.6, .75), c(.3, .3))),
    dna = list(line(c(.25, .35, .65, .75), c(.1, .35, .65, .9)),
               line(c(.75, .65, .35, .25), c(.1, .35, .65, .9)),
               line(c(.3, .7), c(.22, .22)), line(c(.37, .63), c(.4, .4)),
               line(c(.37, .63), c(.6, .6)), line(c(.3, .7), c(.78, .78))),
    matrix = lapply(seq_len(12), function(i) {
      rectangle(((i-1) %% 4 + .5)/4, ((i-1) %/% 4 + .5)/3, .2, .26,
                c('#e6f0ed', '#9fc9bf', color)[(i %% 3)+1])
    }),
    clock = list(circle(.5, .5, .42), line(c(.5, .5, .74), c(.8, .5, .38)),
                 line(c(.5, .5), c(.1, .19)), line(c(.1, .19), c(.5, .5)),
                 line(c(.81, .9), c(.5, .5))),
    cells = list(circle(.3, .65, .21, '#dfebe6'), circle(.72, .6, .23, '#edf3e1'),
                 circle(.5, .2, .18, '#dce9f1'), circle(.3, .65, .07, color),
                 circle(.72, .6, .09, color), circle(.5, .2, .05, color)),
    benchmark = list(line(c(.12, .12, .9), c(.9, .1, .1)),
                     rectangle(.3, .3, .14, .4, color),
                     rectangle(.53, .43, .14, .66, color),
                     rectangle(.76, .55, .14, .9, color)),
    arena = list(rectangle(.5, .58, .9, .66),
                 line(c(.5, .5), c(.25, .08)), line(c(.25, .75), c(.08, .08)),
                 circle(.5, .58, .23), line(c(.5, .5, .65), c(.74, .58, .49))),
    integrate = list(rectangle(.5, .5, .85, .75),
                     line(c(.1, .9), c(.65, .65)), line(c(.37, .37), c(.12, .88)),
                     line(c(.64, .64), c(.12, .88)), line(c(.1, .9), c(.38, .38)))
  )
  grid::gTree(children = do.call(grid::gList, parts))
}
# A database cylinder uses three oval outlines joined by vertical sides.
# Keeping it separate makes the geometry easy to edit.
database_icon <- function(color) {
  oval <- function(y) grid::polygonGrob(
    .5 + .4 * cos(seq(0, 2*pi, length.out = 60)),
    y + .12 * sin(seq(0, 2*pi, length.out = 60)),
    gp = grid::gpar(fill = '#eef3f8', col = color, lwd = 1.5))
  grid::grobTree(grid::rectGrob(.5, .48, .8, .6,
                  gp = grid::gpar(fill = '#eef3f8', col = color, lwd = 1.5)),
                oval(.18), oval(.48), oval(.78))
}
boxes <- tribble(
  ~label,                ~icon,        ~x,  ~y,  ~fill,      ~color,
  'PubMed',              'paper',     .11, .82, '#eef3f8', '#587c9b',
  'Clock\ncuration',     'curate',    .37, .82, '#eef3f8', '#587c9b',
  'Clock\napplication',  'clock',     .64, .82, '#eaf3ef', '#287d80',
  'GEO',                 'database',  .11, .49, '#eef3f8', '#587c9b',
  'Data\nprocessing',    'dna',       .37, .49, '#eaf3ef', '#287d80',
  'Betas',               'matrix',    .64, .49, '#eaf3ef', '#287d80',
  'Cell\ncomposition',   'cells',     .88, .49, '#f3f2e8', '#969151',
  'Metadata\ncuration',  'curate',    .37, .17, '#eef3f8', '#587c9b',
  'Data integration',    'integrate', .64, .17, '#d7ebe4', '#287d80',
  'Benchmark',           'benchmark', .64, -.16, '#eef3f8', '#587c9b',
  'ClockArena\nwebsite', 'arena',     .90, -.16, '#eaf3ef', '#287d80'
)
edges <- tribble(
  ~path, ~x, ~y,
  1, .21, .82, 1, .27, .82,                    # PubMed -> curation
  2, .47, .82, 2, .54, .82,                    # curation -> application
  3, .21, .49, 3, .27, .49,                    # GEO -> processing
  4, .11, .375, 4, .11, .17, 4, .27, .17,      # GEO -> metadata
  5, .47, .49, 5, .54, .49,                    # processing -> betas
  6, .64, .605, 6, .64, .705,                  # betas -> application
  7, .74, .49, 7, .78, .49,                    # betas -> cells
  8, .74, .82, 8, 1.02, .82, 8, 1.02, .17, 8, .74, .17,
  9, .47, .17, 9, .54, .17,                    # metadata -> integration
  10, .88, .375, 10, .88, .23, 10, .74, .23,   # cells -> integration
  11, .64, .055, 11, .64, -.045,               # integration -> benchmark
  12, .74, -.16, 12, .80, -.16                # benchmark -> ClockArena
)
A <- ggplot() +
  geom_path(data = edges, aes(x, y, group = path), linewidth = .45,
            lineend = 'round', arrow = arrow(length = unit(2, 'mm'), type = 'closed'),
            color = '#6b7c83') +
  coord_cartesian(xlim = c(-.01, 1.05), ylim = c(-.30, .97), expand = FALSE) +
  theme_pubr() + no_axes
for (i in seq_len(nrow(boxes))) {
  node <- boxes[i, ]
  card <- grid::roundrectGrob(r = unit(3, 'mm'),
              gp = grid::gpar(fill = node$fill, col = NA))
  icon <- if (node$icon == 'database') database_icon(node$color) else
    workflow_icon(node$icon, node$color)
  A <- A + annotation_custom(card, node$x-.10, node$x+.10, node$y-.115, node$y+.115) +
    annotation_custom(icon, node$x-.032, node$x+.032, node$y+.005, node$y+.085)
}
A <- A + geom_text(data = boxes, aes(x, y = y-.06, label = label), color = '#263b43')
A
# B. All clocks in the local copy of the online Clocks Table -------------------
# Use the saved current Clocks_Table.csv under paper/input.
clock_table <- read_csv(file.path(data_dir, 'Clocks_Table.csv'),
                        col_types = cols(.default = col_character())) %>%
  filter(!is.na(`Catalogue ID`)) %>%
  mutate(year = as.integer(`Publication year`),
         tissue_group = factor(Tissue, levels = names(tissue_colors)))
a <- clock_table %>% count(year, tissue_group, name = 'clocks')
B <- ggplot(a, aes(year, clocks, fill = tissue_group)) +
  theme_pubr() +
  geom_col(width = .8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = tissue_colors) +
  scale_x_continuous(breaks = seq(min(a$year), max(a$year), 2)) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x = 'Publication year', y = 'Number of epigenetic clocks', fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
B
# C. Year versus model size for the selected clocks ----------------------------
# Fit uses actual years; x_display only separates points with identical values.
selected_catalogue <- read_csv('input/06_data_integration/clock_names/selected_clock_catalogue.csv', show_col_types = FALSE)
calculable_clocks <- clock_table %>%
  filter(`Catalogue ID` %in% selected_catalogue$catalogue_id) %>%
  arrange(year, `Internal calculation ID`)

# Supplementary Table 1: source columns B, F, G, I, L, N, O, P, Q, R, S, T.
supplementary_table_1 <- calculable_clocks %>%
  select(`PubMed ID / reference`, Author, `Publication year`,
         `Benchmark / platform clock name`, Tissue, `Number CpGs`, Outcome,
         `Age range (years)`, `Training cohort`, `Training sample size`,
         Algorithm, `Methylation platform`) %>%
  mutate(across(c(`Publication year`, `Number CpGs`), as.numeric))

b <- calculable_clocks %>%
  transmute(catalogue_id = `Catalogue ID`, model_id = `Internal calculation ID`,
            clock = `Benchmark / platform clock name`, year,
            resolved_cpgs = as.numeric(`Number CpGs`)) %>%
  arrange(year, model_id) %>%
  group_by(year, resolved_cpgs) %>%
  mutate(x_display = year + if (n() > 1) seq(-.18, .18, length.out = n()) else 0) %>%
  ungroup()
fit <- lm(log10(resolved_cpgs) ~ year, data = b)
fit_line <- tibble(year = range(b$year)) %>%
  mutate(cpgs = 10^predict(fit, newdata = .))
panel_c_test <- cor.test(b$year, log10(b$resolved_cpgs), method = "pearson")
r <- unname(panel_c_test$estimate)
panel_c_p_label <- sprintf("p = %.2f", panel_c_test$p.value)
C <- ggplot(b, aes(x_display, resolved_cpgs)) +
  theme_pubr() +
  geom_smooth(method = "lm", color = "black", fill = palette_3[["grey"]])+
  geom_point(color = teal, alpha = 0.7, size = 3) +
  scale_y_log10(limits = c(1, 1e6), breaks = c(1, 100, 10000, 1e6), labels = label_comma(),
                expand = expansion(mult = .05)) +
  scale_x_continuous(limits = c(2010, 2027), breaks = c(2012, 2016, 2020, 2024),
                     expand = expansion(mult = .05)) +
  labs(x = 'Publication year', y = 'Required CpGs (log scale)') +
  # Convert panel fractions to data coordinates, including 5% axis expansion.
  # The y position is converted on the log10 scale.
  annotate('text', x = 2010 + 17 * (1.1 * .2 - .05),
           y = 10^(6 * (1.1 * .9 - .05)),
           label = sprintf('R = %.2f', r), hjust = .5, vjust = .5) +
  annotate('text', x = 2010 + 17 * (1.1 * .2 - .05),
           y = 10^(6 * (1.1 * .83 - .05)),
           label = panel_c_p_label, hjust = .5, vjust = .5)
C
# Supplementary Figure 1. Publication year versus dataset-median MAE -----------
# One point per selected clock; MAE is in years, without a log transform.
# Use the current paper benchmark snapshot, with equal dataset weighting in
# each clock's median. Fit and correlation use actual publication years.
age_mae_summary <- read_csv(file.path(data_dir, 'benchmark/clock_summary.csv'),
                            show_col_types = FALSE) %>%
  select(clock, age_datasets, median_age_MAE)
stopifnot(!anyDuplicated(age_mae_summary$clock))
supplementary_figure_1_data <- calculable_clocks %>%
  transmute(catalogue_id = `Catalogue ID`, model_id = `Internal calculation ID`,
            clock = `Benchmark / platform clock name`, year) %>%
  left_join(age_mae_summary, by = 'clock') %>%
  filter(is.finite(year), is.finite(median_age_MAE)) %>%
  arrange(year, model_id) %>%
  group_by(year, median_age_MAE) %>%
  mutate(x_display = year + if (n() > 1) seq(-.18, .18, length.out = n()) else 0) %>%
  ungroup()
supplementary_figure_1_test <- with(supplementary_figure_1_data,
  cor.test(year, median_age_MAE, method = 'pearson'))
mae_y_limits <- c(0, max(supplementary_figure_1_data$median_age_MAE) * 1.05)
Supplementary_Figure_1 <- ggplot(supplementary_figure_1_data,
                                aes(year, median_age_MAE)) +
  theme_pubr() +
  geom_smooth(method = 'lm', color = 'black', fill = palette_3[['grey']]) +
  geom_point(aes(x = x_display), color = teal, alpha = 0.7, size = 3) +
  scale_y_continuous(limits = mae_y_limits, expand = expansion(mult = .05)) +
  scale_x_continuous(limits = c(2010, 2027), breaks = c(2012, 2016, 2020, 2024),
                     expand = expansion(mult = .05)) +
  labs(x = 'Publication year', y = 'Dataset-median MAE (years)') +
  annotate('text', x = 2010 + 17 * (1.1 * .2 - .05),
           y = mae_y_limits[2] * (1.1 * .9 - .05),
           label = sprintf('R = %.2f', unname(supplementary_figure_1_test$estimate)),
           hjust = .5, vjust = .5) +
  annotate('text', x = 2010 + 17 * (1.1 * .2 - .05),
           y = mae_y_limits[2] * (1.1 * .83 - .05),
           label = sprintf('p = %.2f', supplementary_figure_1_test$p.value),
           hjust = .5, vjust = .5)
Supplementary_Figure_1
# D. Most reused CpGs, with an aligned CpG / nearest-gene table -----------------
# Matrix rows are CpGs, columns are internal clock IDs; entries are 0 or 1.
# Nearest-gene annotations are retained as an attribute of the matrix.
cpg_matrix <- readRDS(file.path(data_dir, 'CpG_clock_matrix.rds'))
# Additional top-15 labels use the same nearest curated RefSeq TSS definition.
panel_d_genes <- bind_rows(
  attr(cpg_matrix, 'nearest_gene_annotation') %>% select(cpg, nearest_gene),
  read_csv(file.path(data_dir, 'additional_nearest_genes.csv'),
           show_col_types = FALSE) %>% select(cpg, nearest_gene)) %>%
  distinct(cpg, .keep_all = TRUE)
c <- tibble(cpg = rownames(cpg_matrix), clock_count = rowSums(cpg_matrix)) %>%
  arrange(desc(clock_count), cpg) %>%
  slice_head(n = 15) %>%
  left_join(panel_d_genes, by = 'cpg') %>%
  mutate(row = rev(row_number()))
stopifnot(!anyNA(c$nearest_gene), all(nzchar(c$nearest_gene)))
y_range <- c(.4, nrow(c) + 1.3)
D_table <- ggplot(c, aes(y = row)) +
  geom_text(aes(x = 0, label = cpg), hjust = 0) +
  geom_text(aes(x = .45, label = nearest_gene), hjust = 0) +
  annotate('text', x = c(0, .45), y = nrow(c) + 1,
           label = c('CpG', 'Nearest gene'), hjust = c(0, 0), fontface = 'bold') +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = y_range, expand = c(0, 0)) +
  theme_pubr() + no_axes + theme(plot.margin = margin(10, 0, 5, 5))
D_bars <- ggplot(c, aes(clock_count, row)) +
  theme_pubr() +
  geom_col(aes(fill = clock_count), orientation = 'y', width = .66) +
  scale_fill_gradient(low = colors$blue[1], high = colors$blue[3], guide = 'none') +
  geom_text(aes(label = clock_count), hjust = -.3) +
  scale_x_continuous(limits = c(0, nrow(calculable_clocks)), breaks = seq(0, nrow(calculable_clocks), 20), expand = c(0, 0)) +
  scale_y_continuous(limits = y_range, expand = c(0, 0)) +
  labs(x = 'Number of epigenetic clocks', y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(10, 5, 5, 0))
D <- plot_grid(D_table, D_bars, nrow = 1, rel_widths = c(1, 1.15), align = 'h', axis = 'tb') + theme_pubr() + no_axes
D
# Supplementary Table 2. All CpGs used by the selected clocks ------------------
# Two columns: CpG identifier and number of clocks using that probe.
supplementary_table_2 <- tibble(cpg = rownames(cpg_matrix),
                               clocks_using_probe = rowSums(cpg_matrix)) %>%
  filter(clocks_using_probe > 0) %>%
  arrange(desc(clocks_using_probe), cpg)
write_csv(supplementary_table_2, file.path(out_dir, 'Supplementary_Table_2.csv'),
          na = '')
# E. Mirrored age histograms, stacked by the agreed display health groups -------
# Blank disease is grouped as healthy/normal ONLY for this display, as requested.
# This is not a verified healthy-control classification.
metadata <- read_csv(file.path(data_dir, 'integrated_metadata.csv'),
                     col_types = cols(.default = col_character())) %>%
  mutate(age = as.numeric(age), across(c(disease, treatment, country), str_squish))
e <- metadata %>%
  filter(is.finite(age), sex %in% c('female', 'male')) %>%
  mutate(age_left = floor(age / 5) * 5,
         health_group = if_else(is.na(disease) | disease == '',
                                'Healthy / normal', 'Disease documented')) %>%
  count(age_left, sex, health_group, name = 'records') %>%
  mutate(value = if_else(sex == 'female', records, -records),
         health_group = factor(health_group, c('Healthy / normal', 'Disease documented')))
sex_counts <- e %>% group_by(sex) %>% summarise(n = sum(records), .groups = 'drop')
# Disease starts at zero for both sexes; healthy/normal extends outward.
e <- e %>%
  group_by(age_left, sex) %>%
  arrange(desc(health_group), .by_group = TRUE) %>%
  mutate(outer = cumsum(records), inner = outer - records,
         ymin = if_else(sex == 'male', -outer, inner),
         ymax = if_else(sex == 'male', -inner, outer)) %>% ungroup()
E <- ggplot(e, aes(fill = health_group)) +
  theme_pubr() +
  geom_rect(aes(xmin = age_left, xmax = age_left + 5, ymin = ymin, ymax = ymax),
            color = 'white', linewidth = .1) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_fill_manual(values = c(teal, coral), labels = c('Healthy/normal', 'Disease')) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  scale_y_continuous(labels = function(x) comma(abs(x)), expand = expansion(mult = .2)) +
  coord_cartesian(xlim = c(0, 105)) +
  annotate('text', x = 103, y = Inf, vjust = 1.5, hjust = 1,
           label = paste0('Female\nn = ', comma(sex_counts$n[sex_counts$sex == 'female']))) +
  annotate('text', x = 103, y = -Inf, vjust = -1, hjust = 1,
           label = paste0('Male\nn = ', comma(sex_counts$n[sex_counts$sex == 'male']))) +
  labs(x = 'Chronological age (years)', y = 'Number of samples', fill = NULL)
E
# F. Country counts from integrated metadata; map geometry from the maps package.
# Curation v2 uses Congo as the requested display name for DR Congo.
# Match the generalizability country exclusions, only within this panel.
panel_f_excluded_countries <- c('Mexico', 'Gambia')
f <- metadata %>% filter(!is.na(country), country != '',
                          !country %in% panel_f_excluded_countries) %>%
  count(country, name = 'records') %>%
  mutate(map_country = recode(country, 'United States' = 'USA',
                             'United Kingdom' = 'UK',
                             'Congo' = 'Democratic Republic of the Congo',
                             'Democratic Republic of the Congo' = 'Democratic Republic of the Congo'))
# Wrap polygons at the date line so western fragments are not clipped away.
world <- map_data('world', wrap = c(-180, 180)) %>%
  left_join(f %>% select(map_country, records), by = c('region' = 'map_country'))
F_map <- ggplot(world, aes(long, lat, group = group, fill = records)) +
  geom_polygon(color = palette_3[['grey']], linewidth = .1) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-58, 85), expand = FALSE) +
  scale_fill_gradientn(colors = c('white', palette_3[['orange']]), trans = 'log10',
                       breaks = c(10, 100, 1000, 10000), labels = label_comma(),
                       na.value = '#eeeeee',
                       name = 'Number of samples (log scale)') +
  guides(fill = guide_colorbar(title.position = 'bottom', title.hjust = .5,
                               barwidth = unit(70, 'mm'), barheight = unit(2, 'mm'))) +
  theme_pubr() + no_axes + theme(legend.position = 'bottom', legend.title = element_text(),
                       legend.margin = margin(0))
f_labels <- f %>% arrange(desc(records)) %>%
  mutate(short = recode(country, 'United States' = 'US', 'United Kingdom' = 'UK',
                        'Democratic Republic of the Congo' = 'DR Congo'),
         label = paste0(short, ' (', comma(records), ')'),
         column = (row_number()-1) %/% ceiling(n()/5),
         row = ceiling(n()/5)-(row_number()-1) %% ceiling(n()/5))
# Tighten only the label table; preserve the map and its allocated height.
f_labels <- f_labels %>%
  mutate(column = .70 + .90 * column,
         row = mean(range(row)) + .80 * (row - mean(range(row))))
F_table <- ggplot(f_labels, aes(column, row, label = label)) +
  geom_text(hjust = 0, size = 2.7) +
  coord_cartesian(xlim = c(-.03, 4.95), ylim = c(.4, ceiling(nrow(f)/5) + .5), expand = FALSE) + theme_pubr() + no_axes
F <- plot_grid(F_map, F_table, ncol = 1, rel_heights = c(3, 1)) + theme_pubr() + no_axes
F
# Shared layout for G and H: exact counts, log axis, dotted row separators ------
paired_counts <- function(data, x_limits, ticks, zero_position, rank_by = 'total') {
  rows <- data %>% arrange(desc(.data[[rank_by]]), label) %>% mutate(row = rev(row_number()))
  points <- rows %>% pivot_longer(c(cases, controls), names_to = 'group', values_to = 'n') %>%
    mutate(y = row + if_else(group == 'cases', .23, -.23), hjust = 0)
  # Axis labels sit directly beside the dots; no separate label panel.
  # A zero count is shown at a separate labelled position left of the log range.
  points <- points %>% mutate(x = if_else(n == 0, zero_position, as.numeric(n)),
                             text_x = x * 1.3)
  ggplot(points, aes(x, y, color = group)) +
    geom_hline(yintercept = seq(1.5, nrow(rows)-.5), linetype = 'dotted', color = palette_3[['grey']]) +
    geom_point(size = 3) +
    geom_text(aes(x = text_x, label = comma(n), hjust = hjust), show.legend = FALSE) +
    scale_color_manual(values = c(cases = coral, controls = blue), labels = c('Cases', 'Controls')) +
    scale_x_log10(limits = x_limits,
                 breaks = ticks, labels = function(x) ifelse(x == zero_position, '0', comma(x))) +
    scale_y_continuous(breaks = rows$row, labels = str_wrap(rows$label, 24),
                       limits = c(.4, nrow(rows)+.6), expand = c(0, 0)) +
    labs(x = 'Number of samples (log scale)', y = NULL) + theme_pubr() +
    theme(legend.title = element_blank())
}

# G. Cases carry the disease label; controls have blank disease in the same study.
# Each study's blank-disease records count once per disease present in that study.
# Select and display the top ten disease labels by case count, highest first.
disease_cases <- metadata %>% filter(!is.na(disease), disease != '') %>%
  count(GSE_ID, label = disease, name = 'cases')
disease_controls <- metadata %>% filter(is.na(disease) | disease == '') %>%
  count(GSE_ID, name = 'controls')
g <- disease_cases %>%
  left_join(disease_controls, by = 'GSE_ID') %>%
  mutate(controls = replace_na(controls, 0L)) %>%
  group_by(label) %>% summarise(cases = sum(cases), controls = sum(controls), .groups = 'drop') %>%
  mutate(total = cases + controls) %>% arrange(desc(cases), label) %>% slice_head(n = 10)

# H. Top ten recorded intervention/exposure labels, ranked by case count.
# Keep the integrated labels; do not reuse earlier manually grouped totals.
h <- metadata %>%
  filter(!is.na(treatment), treatment != '', case_control %in% c('case', 'ctrl')) %>%
  mutate(group = if_else(case_control == 'case', 'cases', 'controls')) %>%
  count(label = treatment, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  mutate(total = cases + controls) %>%
  arrange(desc(cases), label) %>% slice_head(n = 10) %>%
  mutate(label = str_wrap(label, 33))
# Use the same log axis for both panels, including any zero-count position.
shared_counts <- c(g$cases, g$controls, h$cases, h$controls)
zero_position <- min(shared_counts[shared_counts > 0]) / 3
shared_x <- if_else(shared_counts == 0, zero_position, as.numeric(shared_counts))
x_limits <- c(min(shared_x) / 1.5, max(shared_x) * 1.8)
x_ticks <- c(1, 10, 100, 1000, 10000)
if (any(shared_counts == 0)) x_ticks <- c(zero_position, x_ticks[x_ticks > zero_position])
G <- paired_counts(g, x_limits, x_ticks, zero_position, rank_by = 'cases')
H <- paired_counts(h, x_limits, x_ticks, zero_position, rank_by = 'cases')
G
H
# Assemble A–H in reading order with cowplot ----------------------------------
figure <- plot_grid(A, B, C, D, E, F, G, H, ncol = 2,
                    labels = LETTERS[1:8], label_size = 12,
                    label_x = .01, label_y = 1,
                    rel_heights = c(1, 1.05, 1.05, 1.3)) + theme_pubr() + no_axes


# The larger canvas accommodates theme_pubr default typography without per-axis sizing.
# Save the combined figure, Supplementary Figure 1 and Supplementary Table 1 ---
invisible(dev.off())
pdf(file.path(out_dir, 'Figure_1.pdf'),
    width = 366 / 25.4, height = 550 / 25.4)
print(figure)
dev.off()

# Export the same assembled figure at the PDF's physical size, with a white background.
ragg::agg_png(file.path(out_dir, 'Figure_1.png'),
    width = 366 / 25.4, height = 550 / 25.4,
    units = 'in', res = 300, background = 'white')
print(figure)
dev.off()

ggsave(file.path(out_dir, 'Supplementary_Figure_1.pdf'), Supplementary_Figure_1,
       width = 183, height = 140, units = 'mm')
write_csv(supplementary_table_1, file.path(out_dir, 'Supplementary_Table_1.csv'))
message('Saved Figure 1 (PDF and PNG), Supplementary Figure 1 and supplementary table values; XLSX is refreshed separately preserving its formatting.')
