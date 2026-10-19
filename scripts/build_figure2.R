#!/usr/bin/env Rscript
# Figure 2: clock agreement and benchmark behaviour.
# Run from the project root: Rscript paper/scripts/build_figure2.R
# Generates Figure 2 (PDF and PNG), Supplementary Figures 2–6 and Supplementary Tables 3–15.
# Scientific inputs come from paper/input. All ages are included.
# Panels E/F exclude Gambia and Mexico, matching the supplied generalizability scope.
# Panels G/H use source-timed biological intervals and direct reliability MAD.

# Libraries ------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(writexl)
dir.create('paper/output', recursive=TRUE, showWarnings=FALSE)

# Exact palette from build_figure1.R; copied here to keep this script standalone.
colors <- list(
  tone = c('#f6f2ed','#dedbca','#c4c0a5','#a59f83','#878264','#5e5948'),
  grey = c('#e5e4e8','#c5c9d6','#959fb3','#6e788c','#425468','#1b2942'),
  olive = c('#f2edb2','#dbdc63','#c4c400','#95a008','#637314','#304215'),
  green = c('#d7e5c5','#9fc978','#5db342','#41912f','#1d6e29','#0e3716'),
  teal = c('#c9e4ef','#96ced3','#48bcbc','#00959f','#006479','#003547'),
  blue = c('#c5e4fb','#9bc9e8','#5495ce','#006eae','#01478c','#002259'),
  purple = c('#e9d3e7','#d1a9ce','#b678b3','#a44990','#792373','#430b4d'),
  red = c('#f5cfc9','#e9a0a4','#db6463','#c5373c','#9c241b','#730c0d'),
  orange = c('#fbdcbc','#f9bd7b','#f29741','#e96700','#b34a00','#832a00'),
  yellow = c('#ffedc1','#f7dc86','#e8c54d','#ca9a23','#9b730a','#685409'),
  skin = c('#f6e5d3','#dcbc9f','#bc9778','#906852','#734d3d','#422a17'))
# Role-specific colours from the shared palette, following Paul Tol's guidance.
# https://sronpersonalpages.nl/~pault/ — qualitative, diverging and sequential roles.
paper_clock_palette <- vapply(colors, function(x) x[3], character(1))
blue <- colors$blue[3]
teal <- colors$teal[3]
coral <- colors$red[3]

#### Panel A ####
set.seed(20260922)
# Shared inputs are loaded at their first use and reused by later panels.
meta <- read_csv('paper/input/integrated_metadata.csv', col_types = cols(.default = col_character()), na = c('', 'NA', 'NaN'), progress = FALSE)
# Read the current analysis-facing names without relying on column positions.
clocks <- read_csv('paper/input/clock_field_registry.csv', show_col_types = FALSE)$clock %>% sort()
stopifnot(length(clocks) == 100, all(clocks %in% names(meta)))
cells <- c('B', 'NK', 'CD4T', 'CD8T', 'Mono', 'Neutro', 'Eosino')

meta <- meta %>%
  mutate(across(all_of(c('age', clocks, cells)), as.numeric), sex = factor(sex), patient_id = na_if(as.character(patient_id), '')) %>%
  filter(is.finite(age))
datasets <- meta %>% group_by(GSE_ID) %>% group_split(.keep = TRUE)

# Shared lm/residuals helper for age/sex/cell adjustment, permitting aliased terms.
# A full-cell model is accepted only with at least 3 residual df.
fit_cell_adjustment <- function(d) {
  y <- d %>% select(all_of(clocks)) %>% as.matrix()
  full_residuals <- matrix(NA_real_, nrow(d), length(clocks), dimnames = list(NULL, clocks))
  shift_rows <- tibble(dataset = character(), clock = character(), n = integer(), mean_abs_shift = double())
  patterns <- apply(is.finite(y), 2, paste0, collapse = '')
  complete_cells <- d %>% transmute(complete = if_all(all_of(cells), is.finite)) %>% pull(complete)

  for (pattern in unique(patterns)) {
    ix <- which(patterns == pattern)
    ok <- is.finite(d$age) & is.finite(y[, ix[1]]) & complete_cells
    s <- d[ok, , drop = FALSE]
    yy <- y[ok, ix, drop = FALSE]
    if (nrow(s) < 10) next

    use_sex <- !anyNA(s$sex) && sum(table(s$sex) >= 2) >= 2
    for (j in seq_along(ix)) {
      s$predicted_age <- yy[, j]
      if (use_sex) {
        reduced_fit <- lm(predicted_age ~ age + sex, data = s)
        full_fit <- lm(predicted_age ~ age + sex + B + NK + CD4T + CD8T + Mono + Neutro + Eosino,
                       data = s)
      } else {
        reduced_fit <- lm(predicted_age ~ age, data = s)
        full_fit <- lm(predicted_age ~ age + B + NK + CD4T + CD8T + Mono + Neutro + Eosino,
                       data = s)
      }
      if (df.residual(full_fit) < 3) next
      residual_reduced <- residuals(reduced_fit)
      residual_full <- residuals(full_fit)
      full_residuals[ok, ix[j]] <- residual_full
      shift_rows <- bind_rows(shift_rows,
        tibble(dataset = d$GSE_ID[1], clock = clocks[ix[j]], n = nrow(s),
               mean_abs_shift = mean(abs(residual_full - residual_reduced))))
    }
  }
  list(residuals = full_residuals, shifts = shift_rows)
}

# Panels A and J share the same blue-red correlation scale.
# Scratch graphics device for plot assembly; does not create a file.
pdf(NULL)
cor_colors <- colorRamp2(c(-1, -.5, 0, .5, 1),
  c(colors$blue[4], colors$blue[2], 'white', colors$red[2], colors$red[4]))

# Dataset-adjusted clock agreement: correlations of full age/sex/cell residuals.
# Exclude constant-age datasets directly; no age-performance eligibility mask.
# Average valid Pearson correlations arithmetically with equal dataset weights.
a_cor_sum <- matrix(0, length(clocks), length(clocks), dimnames = list(clocks, clocks))
a_cor_n <- a_cor_sum

for (d in datasets) {
  if (n_distinct(d$age) < 2) next
  adjusted <- fit_cell_adjustment(d)
  residuals <- adjusted$residuals
  dataset_cor <- suppressWarnings(cor(residuals, use = 'pairwise.complete.obs'))
  pair_counts <- crossprod(is.finite(residuals))
  usable <- is.finite(dataset_cor) & pair_counts >= 10
  a_cor_sum[usable] <- a_cor_sum[usable] + dataset_cor[usable]
  a_cor_n[usable] <- a_cor_n[usable] + 1
}
a_cor <- a_cor_sum / ifelse(a_cor_n > 0, a_cor_n, NA_real_)
# Keep unsupported diagonal entries missing, like unsupported clock pairs.
diag(a_cor)[diag(a_cor_n) > 0] <- 1
a_cluster_matrix <- a_cor
a_cluster_matrix[!is.finite(a_cluster_matrix)] <- 0
a_cluster <- hclust(as.dist(pmax(1 - a_cluster_matrix, 0)), method = 'average')
# Annotation values are copied verbatim from columns L (Tissue) and O (Outcome)
# of the clock table, matched using its column I benchmark clock names.
a_annotations <- read_csv('paper/input/clock_heatmap_annotations.csv',
                          show_col_types = FALSE)
stopifnot(!anyDuplicated(a_annotations$clock),
          all(colnames(a_cor) %in% a_annotations$clock))
a_annotations <- a_annotations[match(colnames(a_cor), a_annotations$clock), ]
stopifnot(!anyNA(a_annotations[, c('Tissue', 'Outcome')]))
a_annotation_colors <- list(
  Tissue = c('Blood' = colors$grey[2], 'Multiple tissues' = colors$blue[3]),
  Outcome = c('Chronological age' = colors$teal[3],
              'Mortality-related age' = colors$red[3],
              'Other age estimate' = colors$purple[3]))
stopifnot(all(a_annotations$Tissue %in% names(a_annotation_colors$Tissue)),
          all(a_annotations$Outcome %in% names(a_annotation_colors$Outcome)))
a_top_annotation <- HeatmapAnnotation(
  Tissue = a_annotations$Tissue, Outcome = a_annotations$Outcome,
  col = a_annotation_colors, simple_anno_size = grid::unit(2.5, 'mm'),
  annotation_name_side = 'left', annotation_name_gp = grid::gpar(fontsize = 9),
  show_legend = FALSE)
# One ordered legend row with matching typography and colour-key heights.
a_legend_style <- list(
  title_position = 'topleft', title_gp = grid::gpar(fontsize = 10),
  labels_gp = grid::gpar(fontsize = 9),
  grid_height = grid::unit(2.5, 'mm'))
a_pearson_legend <- do.call(Legend, c(list(
  title = 'Pearson r', col_fun = cor_colors,
  at = c(-1, -.5, 0, .5, 1), direction = 'horizontal',
  legend_width = grid::unit(28, 'mm')), a_legend_style))
a_tissue_legend <- do.call(Legend, c(list(
  title = 'Tissue', labels = names(a_annotation_colors$Tissue),
  legend_gp = grid::gpar(fill = a_annotation_colors$Tissue),
  grid_width = grid::unit(2.5, 'mm'), nrow = 1), a_legend_style))
a_outcome_legend <- do.call(Legend, c(list(
  title = 'Outcome', labels = names(a_annotation_colors$Outcome),
  legend_gp = grid::gpar(fill = a_annotation_colors$Outcome),
  grid_width = grid::unit(2.5, 'mm'), nrow = 1), a_legend_style))
a_legends <- packLegend(a_pearson_legend, a_tissue_legend, a_outcome_legend,
                        direction = 'horizontal', gap = grid::unit(4, 'mm'))

a_heatmap <- Heatmap(a_cor, name = 'Pearson R', col = cor_colors, na_col = colors$grey[2],
  top_annotation = a_top_annotation,
  width = grid::unit(7.4, 'in'), height = grid::unit(7.4, 'in'),
  show_row_dend = FALSE, show_column_dend = FALSE,
  cluster_rows = as.dendrogram(a_cluster), cluster_columns = as.dendrogram(a_cluster),
  row_dend_side = 'right', column_dend_side = 'top',
  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
  show_column_names = TRUE, column_names_rot = 90,
  column_names_max_height = grid::unit(2.2, 'in'), row_names_side = 'left',
  row_names_gp = grid::gpar(fontsize = 5.8),
  row_names_max_width = grid::unit(2.8, 'in'),
  column_names_gp = grid::gpar(fontsize = 5.8),
  show_heatmap_legend = FALSE)
a_grob <- grid::grid.grabExpr(draw(a_heatmap, heatmap_legend_side = 'top',
  heatmap_legend_list = list(a_legends)), width = 8, height = 8)
A <- ggdraw() + draw_grob(a_grob, scale = 1.3) +
  theme_pubr() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank())
A

# Supplementary Table 3: the exact clock-agreement matrix shown in panel A.
# Pearson residual correlations are averaged equally across eligible datasets.
# Use the heatmap's shared dendrogram order; missing pairs remain blank.
a_display_order <- a_cluster$labels[a_cluster$order]
supplementary_table_3_matrix <- a_cor[a_display_order, a_display_order, drop = FALSE]
supplementary_table_3 <- supplementary_table_3_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column('clock')
write_csv(supplementary_table_3, 'paper/output/Supplementary_Table_3.csv', na = '')
# Supplementary Table 4: mean absolute agreement with the other 99 clocks.
# Take absolute values of the panel A correlations, then average off-diagonal
# entries per clock. Exclude self-correlations and preserve missing pairs.
a_absolute_cor <- abs(a_cor)
diag(a_absolute_cor) <- NA_real_
supplementary_table_4 <- tibble(
  clock = rownames(a_absolute_cor),
  mean_absolute_correlation = rowMeans(a_absolute_cor, na.rm = TRUE),
  other_clocks = rowSums(is.finite(a_absolute_cor))) %>%
  mutate(mean_absolute_correlation = if_else(other_clocks > 0,
    mean_absolute_correlation, NA_real_)) %>%
  arrange(mean_absolute_correlation, clock) %>%
  mutate(rank = if_else(is.finite(mean_absolute_correlation), row_number(), NA_integer_)) %>%
  select(rank, clock, mean_absolute_correlation, other_clocks)
supplementary_table_4 <- bind_rows(supplementary_table_4,
  tibble(rank = NA_integer_, clock = 'Average across all clocks',
         mean_absolute_correlation = mean(supplementary_table_4$mean_absolute_correlation,
                                          na.rm = TRUE),
         other_clocks = NA_integer_))
write_csv(supplementary_table_4 %>% select(-other_clocks),
          'paper/output/Supplementary_Table_4.csv', na = '')
#### Panel B ####
# Rank clocks by the mean absolute change from reduced to cell-adjusted residuals.
b_adjusted <- datasets %>% map(fit_cell_adjustment)
b_data <- read_csv('paper/input/benchmark/cell_influence_summary.csv', show_col_types = FALSE) %>%
  transmute(clock, shift = cell_influence_years, datasets = cell_influence_datasets) %>%
  filter(is.finite(shift)) %>% arrange(desc(shift), clock) %>% mutate(rank = row_number())
B <- ggplot(b_data, aes(rank, shift)) +
  theme_pubr(base_size = 13) +
  theme(plot.margin = margin(3, 4, 3, 3), legend.title = element_blank()) +
  geom_segment(aes(xend = rank, yend = 0), color = colors$grey[3], linewidth = .35) +
  geom_point(shape = 21, fill = colors$yellow[4], color = 'black', stroke = .25, size = 2.4) +
  geom_text(aes(y = shift + .30, label = clock), angle = 90, hjust = 0, vjust = .5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(.02, .55))) +
  scale_x_continuous(expand = expansion(mult = c(.025, .025))) +
  labs(x = 'Clock rank', y = 'Absolute residual shift (years)') +
  theme()
B

#### Supplementary Figure 2: cell-composition influence characteristics ####
# Exploratory clock-level associations; annotations match panel A.
# Higher residual shifts indicate greater influence, not better performance.
s2_catalogue <- read_csv('paper/input/Clocks_Table.csv', show_col_types = FALSE) %>%
  transmute(clock_id = `Catalogue ID`, cpg_count = as.numeric(`Number CpGs`))
stopifnot(!anyDuplicated(s2_catalogue$clock_id), !anyDuplicated(a_annotations$clock))
s2_data <- b_data %>%
  left_join(a_annotations %>% select(clock_id, clock, Tissue, Outcome), by = 'clock') %>%
  left_join(s2_catalogue, by = 'clock_id')
stopifnot(nrow(s2_data) == length(clocks), !anyDuplicated(s2_data$clock),
          all(is.finite(s2_data$shift)), all(is.finite(s2_data$cpg_count)),
          all(s2_data$cpg_count > 0), !anyNA(s2_data$Tissue), !anyNA(s2_data$Outcome))
# Panel A excludes other age estimates from both the plot and the test.
s2_target_data <- s2_data %>%
  filter(Outcome %in% c('Chronological age', 'Mortality-related age'))
s2_target_test <- wilcox.test(shift ~ Outcome, data = s2_target_data, exact = FALSE)
s2_tissue_test <- wilcox.test(shift ~ Tissue, data = s2_data, exact = FALSE)
s2_cpg_test <- cor.test(s2_data$cpg_count, s2_data$shift, method = 'spearman', exact = FALSE)
s2_tests <- tibble(characteristic = c('Training target', 'Training tissue', 'CpG count'),
  test = c('Wilcoxon rank-sum', 'Wilcoxon rank-sum', 'Spearman'),
  n = c(nrow(s2_target_data), nrow(s2_data), nrow(s2_data)), rho = c(NA_real_, NA_real_, unname(s2_cpg_test$estimate)),
  p = c(s2_target_test$p.value, s2_tissue_test$p.value, s2_cpg_test$p.value)) %>%
  mutate(p_BH = p.adjust(p, method = 'BH'))
# Related clocks need not be independent; these tests do not establish independent
# effects of training characteristics or a biological mechanism.
s2_group_plot <- function(variable, title, adjusted_p, data = s2_data) {
  ggplot(data, aes(x = .data[[variable]], y = shift)) +
    geom_boxplot(fill = colors$blue[1], outlier.shape = NA, width = .55) +
    geom_point(size = 2, alpha = .7,
                color = colors$blue[4], position = position_jitter(width = .12, height = 0, seed = 42)) +
    theme_pubr(base_size = 12) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 17)) +
    labs(x = title, y = 'Cell-composition influence (years)',
         subtitle = sprintf('Wilcoxon; BH-adjusted p = %.2g', adjusted_p))
}
s2_target_plot <- s2_group_plot('Outcome', 'Training target', s2_tests$p_BH[1],
                                data = s2_target_data)
s2_tissue_plot <- s2_group_plot('Tissue', 'Training tissue', s2_tests$p_BH[2])
s2_cpg_plot <- ggplot(s2_data, aes(cpg_count, shift)) +
  geom_point(color = colors$blue[4], alpha = .7, size = 2) +
  scale_x_log10() + theme_pubr(base_size = 12) +
  labs(x = 'CpG count (log scale)', y = 'Cell-composition influence (years)',
       subtitle = sprintf('Spearman rho = %.2f; BH-adjusted p = %.2g',
                          s2_cpg_test$estimate, s2_tests$p_BH[3]))
supplementary_figure_2 <- plot_grid(s2_target_plot, s2_tissue_plot, s2_cpg_plot,
                                  nrow = 1, labels = c('A', 'B', 'C'))
pdf('paper/output/Supplementary_Figure_2.pdf', width = 15, height = 5)
print(supplementary_figure_2)
dev.off()

#### Supplementary Figure 3: individual cell-type adjustment ####
# Same complete samples and eligibility as panel B. Single-cell shifts overlap;
# they are not independent contributions and must not be summed.
s3_by_dataset <- map_dfr(datasets, function(d) {
  y <- as.matrix(d %>% select(all_of(clocks)))
  complete_cells <- d %>% transmute(ok = if_all(all_of(cells), is.finite)) %>% pull(ok)
  patterns <- apply(is.finite(y), 2, paste0, collapse = '')
  map_dfr(unique(patterns), function(pattern) {
    ix <- which(patterns == pattern)
    ok <- is.finite(d$age) & complete_cells & is.finite(y[, ix[1]])
    s <- d[ok, , drop = FALSE]
    shifts <- matrix(NA_real_, length(ix), length(cells))
    full_shift <- rep(NA_real_, length(ix))
    status <- 'fewer_than_10_complete_samples'
    use_sex <- FALSE
    if (nrow(s) >= 10) {
      use_sex <- !anyNA(s$sex) && sum(table(s$sex) >= 2) >= 2
      yy <- y[ok, ix, drop = FALSE]
      status <- 'fewer_than_3_full_model_residual_df'
      for (i in seq_along(ix)) {
        s$predicted_age <- yy[, i]
        if (use_sex) {
          reduced_fit <- lm(predicted_age ~ age + sex, data = s)
          full_fit <- lm(predicted_age ~ age + sex + B + NK + CD4T + CD8T + Mono + Neutro + Eosino,
                         data = s)
        } else {
          reduced_fit <- lm(predicted_age ~ age, data = s)
          full_fit <- lm(predicted_age ~ age + B + NK + CD4T + CD8T + Mono + Neutro + Eosino,
                         data = s)
        }
        if (df.residual(full_fit) < 3) next
        reduced <- residuals(reduced_fit)
        full_shift[i] <- mean(abs(reduced - residuals(full_fit)))
        for (j in seq_along(cells)) {
          s$cell_fraction <- s[[cells[j]]]
          if (use_sex) {
            single_fit <- lm(predicted_age ~ age + sex + cell_fraction, data = s)
          } else {
            single_fit <- lm(predicted_age ~ age + cell_fraction, data = s)
          }
          shifts[i, j] <- mean(abs(reduced - residuals(single_fit)))
        }
        status <- 'ok'
      }
    }
    tibble(dataset = as.character(d$GSE_ID[1]),
      clock = rep(clocks[ix], each = length(cells)), cell = rep(cells, length(ix)),
      n = nrow(s), adjustment = if (use_sex) 'age + sex' else 'age',
      shift_years = as.vector(t(shifts)), full_shift_years = rep(full_shift, each = length(cells)),
      status = status)
  })
})
s3_summary <- s3_by_dataset %>% group_by(clock, cell) %>%
  summarise(mean_shift_years = if (any(status == 'ok')) mean(shift_years[status == 'ok']) else NA_real_,
            datasets = sum(status == 'ok'), .groups = 'drop')
s3_full_check <- s3_by_dataset %>% distinct(dataset, clock, full_shift_years, status) %>%
  filter(status == 'ok') %>% group_by(clock) %>%
  summarise(recalculated = mean(full_shift_years), n_datasets = n(), .groups = 'drop') %>%
  left_join(b_data, by = 'clock')
stopifnot(nrow(s3_summary) == length(clocks) * length(cells),
          nrow(s3_full_check) == length(clocks),
          all(s3_full_check$n_datasets == s3_full_check$datasets),
          max(abs(s3_full_check$recalculated - s3_full_check$shift)) < 1e-8)
s3_order <- b_data %>% arrange(desc(shift), clock) %>% pull(clock)
s3_matrix <- matrix(NA_real_, length(clocks), length(cells), dimnames = list(s3_order, cells))
s3_matrix[cbind(match(s3_summary$clock, s3_order), match(s3_summary$cell, cells))] <- s3_summary$mean_shift_years
s3_annotation <- rowAnnotation('All cells (years)' = anno_barplot(
  b_data$shift[match(s3_order, b_data$clock)], width = grid::unit(25, 'mm'),
  gp = grid::gpar(fill = colors$blue[3], col = NA)))
supplementary_figure_3 <- Heatmap(s3_matrix, name = 'Shift (years)',
  col = colorRamp2(c(0, max(s3_matrix) / 2, max(s3_matrix)), c('#fff7ec', '#fc8d59', '#7f0000')),
  cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = grid::gpar(fontsize = 8),
  column_labels = c('B cells', 'NK cells', 'CD4 T cells', 'CD8 T cells', 'Monocytes', 'Neutrophils', 'Eosinophils'),
  column_names_rot = 45, right_annotation = s3_annotation, na_col = '#cccccc',
  column_title = 'Residual shift from adjusting for each cell type individually')
pdf('paper/output/Supplementary_Figure_3.pdf', width = 10, height = 22)
draw(supplementary_figure_3)
dev.off()

# Supplementary Table 5: panel B residual shifts, ranked lowest to highest.
# Shift is the mean absolute residual change in years, equally weighted across
# eligible datasets. This table's ascending rank does not change panel B order.
supplementary_table_5 <- b_data %>%
  arrange(shift, clock) %>%
  mutate(rank = row_number()) %>%
  transmute(rank, clock, shift_years = shift)
write_csv(supplementary_table_5, 'paper/output/Supplementary_Table_5.csv', na = '')
#### Panel C: MAE by chronological age group ####
# Same selection and aggregation as the approved standalone preview.
# Each clock uses its available samples; this is not a matched-sample comparison.
age_mae_by_dataset <- read_csv('paper/input/benchmark/age_range_by_dataset.csv', show_col_types = FALSE)
age_mae_curves <- read_csv('paper/input/benchmark/age_range_summary.csv', show_col_types = FALSE) %>%
  filter(mae_status == 'ok', is.finite(mae)) %>%
  transmute(age_bin, clock, mae, datasets = mae_datasets, samples = mae_samples) %>%
  arrange(age_bin, mae, clock) %>%
  group_by(age_bin) %>% mutate(rank = row_number()) %>% ungroup()
age_mae_selected <- age_mae_curves %>% filter(rank <= 3) %>% pull(clock) %>% unique()
age_mae_plot_data <- age_mae_curves %>%
  filter(clock %in% age_mae_selected) %>%
  mutate(age_midpoint = age_bin + 5) %>%
  arrange(clock, age_bin) %>% group_by(clock) %>%
  mutate(segment = cumsum(c(TRUE, diff(age_bin) != 10))) %>% ungroup()
age_mae_ends <- age_mae_plot_data %>% group_by(clock) %>% slice_tail(n = 1) %>% ungroup()
# Explicit one-to-one clock colours: no cycling or repeated hues at the same shade.
# All swatches belong to the project's palette; retain stable assignments across bins.
age_mae_palette <- c(
  'ABEC' = colors$purple[5], 'BLUP' = colors$orange[4],
  'Childhood CpG clock' = colors$blue[5], 'Childhood PCA clock' = colors$teal[3],
  'Childhood SL PCA clock' = colors$red[3], 'DNAmFitAge' = colors$grey[5],
  'GOLD-R DNAmAge' = colors$yellow[4], 'IntrinClock' = colors$green[4],
  'MAPLE' = colors$olive[4], 'PCGrimAge' = colors$skin[3],
  'Wu (111 CpG)' = colors$purple[3], 'Zhang (514 CpG)' = colors$blue[3],
  'eABEC' = colors$red[5])
stopifnot(all(age_mae_selected %in% names(age_mae_palette)),
          !anyDuplicated(unname(age_mae_palette[age_mae_selected])))
age_mae_ends <- age_mae_ends %>% mutate(label_color = unname(age_mae_palette[clock]))
C <- ggplot(age_mae_plot_data, aes(age_midpoint, mae, color = clock)) +
  geom_line(aes(group = interaction(clock, segment)), linewidth = .75, alpha = 1) +
  geom_point(shape = 16, size = 2.1) +
  ggrepel::geom_text_repel(data = age_mae_ends, aes(label = clock, color = I(label_color)), nudge_x = 1.5,
    direction = 'y', hjust = 0, size = 3.4, box.padding = .15,
    segment.size = .25, min.segment.length = 0, max.overlaps = Inf,
    max.time = 3, seed = 42) +
  scale_color_manual(values = age_mae_palette) +
  scale_x_continuous(breaks = sort(unique(age_mae_plot_data$age_midpoint)),
    labels = function(x) paste0(x - 5, '-', x + 4),
    expand = expansion(mult = c(.025, .16))) +
  scale_y_continuous(expand = expansion(mult = c(.03, .08))) +
  labs(x = 'Chronological age group (years)', y = 'Mean absolute error (years)') +
  theme_pubr(base_size = 13) +
  theme(legend.position = 'none', axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        plot.margin = margin(3, 4, 3, 3))
C
# Supplementary Table 6: all-age R-squared and MAE for every selected clock.
# Keep the saved benchmark's unweighted dataset medians and eligibility scope.
# Separate ranks: higher R-squared is better; lower MAE (years) is better.
supplementary_table_6 <- read_csv('paper/input/benchmark/clock_summary.csv',
                                show_col_types = FALSE) %>%
  filter(clock %in% clocks) %>%
  transmute(clock,
            median_r_squared = median_age_r_squared,
            median_mae_years = median_age_MAE) %>%
  mutate(r_squared_rank = min_rank(desc(median_r_squared)),
         mae_rank = min_rank(median_mae_years)) %>%
  arrange(r_squared_rank, clock) %>%
  select(clock, median_r_squared, r_squared_rank, median_mae_years, mae_rank)
stopifnot(nrow(supplementary_table_6) == length(clocks),
          !anyDuplicated(supplementary_table_6$clock))
write_csv(supplementary_table_6, 'paper/output/Supplementary_Table_6.csv', na = '')
# Supplementary Table 7: R-squared and MAE by ten-year age group.
# These are equal-dataset means (not the all-age medians in Table 6).
# Preserve metric-specific ranks; unsupported groups remain blank in the export.
# Eligibility statuses remain available in the benchmark source table.
supplementary_table_7 <- read_csv('paper/input/benchmark/age_range_summary.csv',
                                show_col_types = FALSE) %>%
  filter(clock %in% clocks) %>%
  arrange(age_lower, r_squared_rank, clock) %>%
  transmute(age_group = paste0(age_lower, '-', age_upper_exclusive - 1),
            clock, mean_r_squared = pearson_r_squared, r_squared_rank,
            mean_mae_years = mae, mae_rank)
write_csv(supplementary_table_7, 'paper/output/Supplementary_Table_7.csv', na = '')
#### Panel D ####
# Highest-mean-AUC dataset per disease; show the top 15 diseases.
d_auc <- read_csv('paper/input/benchmark/disease_auc_by_dataset.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
d_auc <- d_auc %>% filter(eligible, is.finite(auc))
d_top <- d_auc %>%
  group_by(dataset_id, disease) %>%
  summarise(mean_auc = mean(auc), .groups = 'drop') %>%
  arrange(desc(mean_auc), disease, dataset_id) %>%
  distinct(disease, .keep_all = TRUE) %>%
  slice_head(n = 15) %>%
  mutate(label = paste(str_wrap(disease, 38), dataset_id, sep = '\n'))
d_display <- d_auc %>% inner_join(d_top, by = c('dataset_id', 'disease'))
d_levels <- rev(d_top$label)
D <- ggplot(d_display, aes(auc, factor(label, levels = d_levels))) +
  theme_pubr(base_size = 13, border = TRUE) +
  geom_vline(xintercept = .5, color = colors$grey[3], linetype = 2, linewidth = .4) +
  geom_point(position = position_jitter(width = 0, height = .14, seed = 42), color = colors$orange[3], alpha = .65, size = 1.1) +
  geom_point(data = d_top, aes(x = mean_auc), color = 'black', shape = 18, size = 2.2) +
  geom_text(data = d_top, aes(x = .985, label = sprintf('%.3f', mean_auc)), hjust = 1, size = 3.1, color = colors$grey[5]) +
  scale_x_continuous(breaks = seq(0, 1, .2), limits = c(0, 1),
                     expand = expansion(mult = c(.03, .05))) +
  scale_y_discrete(expand = expansion(add = .6)) +
  labs(x = 'AUC', y = NULL) +
  theme(axis.text.y = element_text(size = 9), legend.position = 'none',
        axis.title.x = element_text(margin = margin(t = 7)),
        plot.margin = margin(3, 4, 10, 3))
D

# Supplementary Table 8: all eligible clock-dataset-disease AUC pairs.
# Use the benchmark's adjusted-age-acceleration AUC, highest first.
# Include all valid pairs, beyond the top 15 disease examples in panel D.
supplementary_table_8 <- d_auc %>%
  select(clock, dataset_id, disease, auc) %>%
  arrange(desc(auc), clock, dataset_id, disease)
write_csv(supplementary_table_8, 'paper/output/Supplementary_Table_8.csv', na = '')
# Supplementary Table 9: benchmark Cox HRs per 1 SD for all three mortality cohorts.
# Dataset estimates retain their 95% CIs, deaths, censoring and diagnostic flags.
# The benchmark summary is a sample-count-weighted median of valid dataset HRs.
# At exact half weight, adjacent HRs are averaged; no pooled CI is attached.
mortality_dataset_results <- read_csv('paper/input/benchmark/mortality_hr_by_dataset.csv',
                                     show_col_types = FALSE)
mortality_clock_results <- read_csv('paper/input/benchmark/mortality_summary.csv',
                                   show_col_types = FALSE)
stopifnot(all(mortality_dataset_results$hr_unit == 'per 1 SD within dataset and clock'))
supplementary_table_9 <- mortality_dataset_results %>%
  filter(clock %in% clocks) %>%
  left_join(mortality_clock_results %>%
              select(clock, weighted_median_hazard_ratio, weighting_samples, mortality_rank),
            by = 'clock') %>%
  arrange(mortality_rank, clock, dataset_id) %>%
  select(clock, dataset_id, hazard_ratio, conf_low, conf_high, n, dead, alive,
         weighted_median_hazard_ratio, weighting_samples, mortality_rank,
         p_value, ph_p_value, ph_status, adjustment_model, status)
write_csv(supplementary_table_9, 'paper/output/Supplementary_Table_9.csv', na = '')
#### Panel E ####
general_excluded_countries <- c('Gambia', 'Mexico')
# Dataset-country MAE heatmap. A black border marks each dataset-country minimum.
# Eligibility: known country, age 0-120, finite clock estimate within +/-500 years,
# and at least three usable samples per dataset-country-clock. MAE needs no age variation.
e_dataset_country <- datasets %>%
  map_dfr(function(d) {
    d %>%
      select(country, age, all_of(clocks)) %>%
      filter(!is.na(country), !country %in% general_excluded_countries, between(age, 0, 120)) %>%
      pivot_longer(all_of(clocks), names_to = 'clock', values_to = 'estimate') %>%
      filter(is.finite(estimate), abs(estimate) <= 500) %>%
      group_by(country, clock) %>%
      summarise(mae = mean(abs(estimate - age)), n = n(), .groups = 'drop') %>%
      filter(n >= 3) %>%
      mutate(dataset = d$GSE_ID[1])
  })
e_comparisons <- e_dataset_country %>%
  group_by(country, dataset) %>%
  mutate(minimum_mae = min(mae), winner = abs(mae - minimum_mae) <= 1e-10) %>%
  ungroup()
e_wins <- e_comparisons %>%
  group_by(clock) %>%
  summarise(wins = sum(winner), mean_dataset_mae = mean(mae), .groups = 'drop') %>%
  arrange(desc(wins), mean_dataset_mae, clock)
e_winning_clocks <- e_wins %>% filter(wins > 0) %>% pull(clock)
e_country_order <- e_dataset_country %>%
  group_by(country, clock) %>%
  summarise(mae = mean(mae), .groups = 'drop') %>%
  group_by(country) %>%
  summarise(mean_mae = mean(mae), .groups = 'drop') %>%
  arrange(mean_mae) %>%
  pull(country)
e_pair_order <- e_comparisons %>%
  distinct(country, dataset) %>%
  mutate(country_order = match(country, e_country_order), pair = paste(dataset, country, sep = ' | ')) %>%
  arrange(country_order, dataset)
e_display <- e_comparisons %>%
  filter(clock %in% e_winning_clocks) %>%
  mutate(pair = paste(dataset, country, sep = ' | '))
e_matrix <- e_display %>%
  select(clock, pair, mae) %>%
  pivot_wider(names_from = pair, values_from = mae) %>%
  column_to_rownames('clock') %>%
  as.matrix()
e_matrix <- e_matrix[e_winning_clocks, e_pair_order$pair, drop = FALSE]
e_winner_matrix <- e_display %>%
  select(clock, pair, winner) %>%
  pivot_wider(names_from = pair, values_from = winner, values_fill = FALSE) %>%
  column_to_rownames('clock') %>%
  as.matrix()
e_winner_matrix <- e_winner_matrix[e_winning_clocks, e_pair_order$pair, drop = FALSE]
e_winner_matrix[is.na(e_winner_matrix)] <- FALSE
e_row_labels <- paste0( rownames(e_matrix), ' (', e_wins$wins[match(rownames(e_matrix), e_wins$clock)], ')')
e_ticks <- c(0, 2, 5, 10, 20)
# colorRamp2 clamps values above the final break to the darkest colour.
e_colors <- colorRamp2(seq(0, log1p(20), length.out = 4), colors$blue[1:4])
# Capture panel E at its final two-column, one-row size.
# Flexible heatmap dimensions fill the available space after names and legend;
# row and column cell sizes are intentionally independent.
figure_width <- 22
square_side <- figure_width / 3
second_row_height <- 6.25
e_heatmap <- Heatmap(log1p(e_matrix), name = 'MAE', col = e_colors, na_col = colors$grey[2],
  cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = 'left',
  height = grid::unit(3.0, 'in'),
  row_labels = e_row_labels, row_names_gp = grid::gpar(fontsize = 9),
  row_names_max_width = grid::unit(2.8, 'in'),
  column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 90,
  column_names_max_height = grid::unit(2.4, 'in'),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (e_winner_matrix[i, j]) {
      grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = NA, col = 'black', lwd = 1.1))
    }
  },
  heatmap_legend_param = list(title = 'MAE (years)', at = log1p(e_ticks),
    labels = c('0', '2', '5', '10', '>20'), direction = 'horizontal',
    title_position = 'topcenter', legend_width = grid::unit(2.5, 'in'),
    title_gp = grid::gpar(fontsize = 10), labels_gp = grid::gpar(fontsize = 9)))
e_grob <- local({
  previous_padding <- ComplexHeatmap::ht_opt('HEATMAP_LEGEND_PADDING')
  on.exit(ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = previous_padding))
  ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = grid::unit(10, 'mm'))
  grid::grid.grabExpr(draw(e_heatmap, heatmap_legend_side = 'top'),
                     width = 2 * square_side, height = second_row_height)
})
E <- ggdraw(e_grob) + theme_pubr() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank())
E
# Supplementary Table 10: median dataset MAE for every country-clock pair.
# Rank clocks within each country, and retain the overall country-balanced rank.
country_median_mae <- read_csv('paper/input/benchmark/country_clock_median_mae.csv',
                              show_col_types = FALSE)
country_median_summary <- read_csv('paper/input/benchmark/clock_generalizability_summary.csv',
                                  show_col_types = FALSE)
supplementary_table_10 <- country_median_mae %>%
  filter(clock %in% clocks) %>%
  left_join(country_median_summary %>%
              select(clock, median_country_median_mae, generalizability_rank), by = 'clock') %>%
  arrange(median_dataset_mae, country, clock) %>%
  transmute(country, clock, median_mae_years = median_dataset_mae,
            rank_within_country = country_rank,
            overall_median_mae_years = median_country_median_mae,
            overall_clock_rank = generalizability_rank)
write_csv(supplementary_table_10, 'paper/output/Supplementary_Table_10.csv', na = '')
#### Panel F ####
# Show dataset-country errors for the clock that wins the most minimum-MAE cells.
# Ties in win count are resolved by lower mean MAE, then alphabetical clock name.
f_clock <- e_wins$clock[1]
f_points <- e_dataset_country %>% filter(clock == f_clock)
f_country <- f_points %>%
  group_by(country) %>%
  summarise(mean_mae = mean(mae), .groups = 'drop') %>%
  arrange(mean_mae)
f_country_order <- f_country$country
F <- ggplot(f_points, aes(mae, factor(country, levels = rev(f_country_order)))) +
  theme_pubr(base_size = 13) +
  theme(plot.margin = margin(3, 4, 3, 3), legend.title = element_blank()) +
  geom_point(position = position_jitter(width = 0, height = .15, seed = 42), color = colors$green[3], size = 1.7, alpha = .8) +
  geom_segment(data = f_country, inherit.aes = FALSE, aes(x = mean_mae, xend = mean_mae, y = as.numeric(factor(country, levels = rev(f_country_order))) - .28,
                   yend = as.numeric(factor(country, levels = rev(f_country_order))) + .28),
               color = colors$green[5], linewidth = .8) +
  labs(x = 'Dataset MAE (years)', y = NULL) +
  theme(axis.text.y = element_text(size = 10)) +
  theme()
F
#### Panel G ####
# Use the benchmark's source-timed biological intervals and direct MAD formula.
g_source <- read_csv('paper/input/benchmark/all_datasets_clock_deviations.csv',
                     show_col_types = FALSE, progress = FALSE)
g_pairs <- g_source %>%
  transmute(dataset = dataset_id, participant = participant_id, clock,
            baseline = sample_1, followup = sample_2,
            elapsed = interval_years, deviation, timing_certainty)
g_levels <- c('<3 h', '3-6 h', '6-48 h', '2 d-2 mo', '2-6 mo', '6-12 mo',
              '1-2 yr', '2-3 yr', '3-5 yr', '5-7 yr', '7-10 yr', '10-14 yr')
g_breaks <- c(0, 3/(24*365.25), 6/(24*365.25), 48/(24*365.25),
              2/12, 6/12, 1, 2, 3, 5, 7, 10, 14)
g_pairs <- g_pairs %>%
  filter(is.finite(deviation), elapsed > 0) %>%
  # Rounding removes CSV precision differences at exact scheduled boundaries.
  mutate(time_bin = cut(round(elapsed,12), breaks = round(g_breaks,12),
    labels = g_levels, right = FALSE, include.lowest = TRUE))
stopifnot(!anyNA(g_pairs$time_bin))
g_clock_coverage <- g_pairs %>%
  group_by(time_bin,clock) %>%
  summarise(participants=n_distinct(dataset,participant), pairs=n(),
            datasets=n_distinct(dataset), .groups='drop') %>%
  mutate(included=participants>=6 & clock!='Schmidt (7 CpG)',
         exclusion_reason=case_when(
           clock=='Schmidt (7 CpG)' ~ 'excluded from Panel G by user instruction (2 October 2026)',
           participants<6 ~ 'fewer than six people with valid paired clock estimates',
           TRUE ~ ''))
g_data <- g_pairs %>%
  semi_join(g_clock_coverage %>% filter(included),by=c('time_bin','clock')) %>%
  group_by(time_bin, dataset, participant, clock) %>%
  summarise(mad = mean(abs(deviation)), .groups = 'drop') %>%
  group_by(time_bin, dataset, clock) %>%
  summarise(mad = mean(mad), .groups = 'drop') %>%
  group_by(time_bin, clock) %>%
  summarise(mad = mean(mad), .groups = 'drop')
g_counts <- g_pairs %>%
  group_by(time_bin, .drop = FALSE) %>%
  summarise(participants = n_distinct(dataset, participant), datasets = n_distinct(dataset), .groups = 'drop')
stopifnot(nrow(g_counts)==12,all(g_counts$participants>=6))
G <- ggplot(g_data, aes(time_bin, mad)) +
  theme_pubr(base_size = 13) +
  theme(plot.margin = margin(3, 4, 12, 3), legend.title = element_blank()) +
  geom_boxplot(fill = colors$purple[2], color = colors$grey[5], outlier.shape = NA, width = .55) +
  geom_point(position = position_jitter(width = .17, seed = 42), color = colors$purple[3], size = .9, alpha = .4) +
  scale_y_continuous(expand=expansion(mult=c(.03,.10))) +
  scale_x_discrete(drop = FALSE, labels = function(x) {
    counts <- g_counts[match(x, g_counts$time_bin), ]
    paste0(x, '\nn=', counts$participants, '\nk=', counts$datasets)
  }) +
  labs(x = 'Collection interval (n = people; k = datasets)', y = 'Mean absolute difference (years)') +
  theme(axis.text.x = element_text(size = 8,angle=45,hjust=1),
        axis.title.x = element_text(margin = margin(t = 8))) +
  theme()
G
# Supplementary Figure 4: across-clock spread of signed deviations by dataset/bin.
# Keep Panel G's clock eligibility and participant-within-dataset averaging.
# LOESS uses ordered bins; its approximate band does not adjust for repeated
# dataset contributions. This separate export does not change Panel G.
g_signed_dataset_clock <- g_pairs %>%
  semi_join(g_clock_coverage %>% filter(included), by = c('time_bin', 'clock')) %>%
  group_by(time_bin, dataset, participant, clock) %>%
  summarise(mean_signed_deviation = mean(deviation), .groups = 'drop') %>%
  group_by(time_bin, dataset, clock) %>%
  summarise(mean_signed_deviation = mean(mean_signed_deviation), .groups = 'drop')
g_signed_spread <- g_signed_dataset_clock %>%
  group_by(time_bin, dataset) %>%
  summarise(iqr_years = IQR(mean_signed_deviation), clocks = n(), .groups = 'drop') %>%
  mutate(interval_order = match(as.character(time_bin), g_levels))
stopifnot(all(is.finite(g_signed_spread$iqr_years)),
          !anyNA(g_signed_spread$interval_order),
          !anyDuplicated(g_signed_spread[c('time_bin', 'dataset')]))
g_signed_spread_plot <- ggplot(g_signed_spread, aes(interval_order, iqr_years)) +
  geom_smooth(method = 'loess', formula = y ~ x, span = .75,
    se = TRUE, level = .95, color = colors$purple[5],
    fill = colors$purple[2], alpha = .25, linewidth = .9) +
  geom_point(position = position_jitter(width = .15, height = 0, seed = 42),
    color = colors$purple[3], size = 3, alpha = .85) +
  scale_x_continuous(breaks = seq_along(g_levels), labels = g_levels) +
  theme_classic(base_size = 13) +
  labs(x = 'Collection interval (ordered bins)',
    y = 'Across-clock interquartile range (years)',
    subtitle = 'Dataset-interval points; LOESS fit with approximate 95% confidence band') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
pdf('paper/output/Supplementary_Figure_4.pdf', width = 9, height = 6)
print(g_signed_spread_plot)
dev.off()

#### Panel H ####
# GSE227809 one-week reliability.
# Reference: comment_clinical_trials/draft/v3/Figure1.pdf, Panel A.
# Select the lowest 15 MADs within the one-week cross-session T1 comparison.
h_points <- read_csv(
  'paper/input/benchmark/all_datasets_clock_deviations.csv',
  show_col_types = FALSE) %>%
  filter(dataset_id == 'GSE227809', session=='cross_session_T1',
         is.finite(deviation),
         round(interval_years,12)==round(7/365.25,12)) %>%
  mutate(absolute_difference = abs(deviation))
h_all_summary <- h_points %>% group_by(clock) %>%
  summarise(mean_absolute_deviation=mean(absolute_difference),
            reliability_coefficient_95=1.96*sd(deviation),
            mean_absolute_error=mean(c(abs(clock_1-age_1),abs(clock_2-age_2))),
            mean_error=mean(c(clock_1-age_1,clock_2-age_2)),
            age_estimates=2L*n(),
            participants=n_distinct(participant_id),.groups='drop') %>%
  arrange(mean_absolute_deviation,clock)
h_summary <- h_all_summary %>% slice_head(n=15)
stopifnot(nrow(h_summary)==15,all(h_summary$participants>=6))
h_points <- h_points %>% filter(clock %in% h_summary$clock)
h_mean_age <- h_points %>% distinct(participant_id,age_1) %>%
  summarise(mean_age=mean(age_1)) %>% pull(mean_age)
# Subtract the approximate seven-day interval; A first by user instruction.
stopifnot(all(h_points$interval_years > 0),
          all(abs(h_points$expected_clock_change - h_points$interval_years) < 1e-12))
h_check <- h_points %>% group_by(clock) %>%
  summarise(recomputed_mad = mean(absolute_difference), .groups = 'drop') %>%
  inner_join(h_summary, by = 'clock')
stopifnot(all(abs(h_check$recomputed_mad - h_check$mean_absolute_deviation) < 1e-8))
h_levels <- rev(h_summary$clock)
h_summary <- h_summary %>% mutate(clock = factor(clock, levels = h_levels))
h_points <- h_points %>% mutate(clock = factor(clock, levels = h_levels))
h_table <- h_summary %>%
  transmute(clock, MAD = mean_absolute_deviation,
            RC95 = reliability_coefficient_95, ME=mean_error) %>%
  pivot_longer(c(MAD, RC95, ME), names_to = 'measure', values_to = 'value') %>%
  mutate(table_x = case_when(measure=='MAD' ~ 3.85,measure=='RC95' ~ 4.45,TRUE ~ 5.10),
         display = sprintf('%.2f', value))
h_colors <- setNames(rep(colors$skin[4], length(h_levels)), h_levels)
H <- ggplot(h_points, aes(absolute_difference, clock, color = clock)) +
  geom_point(position = position_jitter(width = 0, height = .14, seed = 227809),
             alpha = .18, size = 1.5, show.legend = FALSE) +
  geom_segment(data = h_summary,
    aes(x = mean_absolute_deviation, xend = reliability_coefficient_95,
        y = clock, yend = clock), linewidth = .5, show.legend = FALSE) +
  geom_point(data = h_summary,
    aes(x = reliability_coefficient_95, y = clock, shape = 'RC95'),
    inherit.aes = FALSE, color = 'black', size = 3.4) +
  geom_point(data = h_summary,
    aes(x = mean_absolute_deviation, y = clock, fill = clock, shape = 'MAD'),
    inherit.aes = FALSE, color = 'black', stroke = .4, size = 2.5) +
  scale_color_manual(values = h_colors, guide = 'none') +
  scale_fill_manual(values = h_colors, guide = 'none') +
  scale_shape_manual(values = c(MAD = 21, RC95 = 124), guide = 'none') +
  geom_text(data = h_table, aes(x = table_x, y = clock, label = display),
            inherit.aes = FALSE, size = 3.1, color = 'black') +
  annotate('point', x = 3.65, y = 15.75, shape = 21, fill = 'white',
           color = 'black', size = 2.5) +
  annotate('text', x = 3.85, y = 15.75, label = 'MAD', size = 3.2) +
  annotate('point', x = 4.22, y = 15.75, shape = 124,
           color = 'black', size = 3.4) +
  annotate('text', x = 4.45, y = 15.75, label = 'RC95', size = 3.2) +
  annotate('text', x = 5.10, y = 15.75, label = 'ME', size = 3.2) +
  scale_x_continuous(breaks=seq(0,3.5,.5), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand=expansion(add=c(.6,1.3))) +
  coord_cartesian(xlim=c(0,3.5), clip = 'off') +
  labs(x = 'Reliability error (years)', y = NULL,
       title=sprintf('GSE227809: 1 week (n=11; mean age %.1f years)',h_mean_age)) +
  theme_pubr(base_size = 13) +
  theme(axis.text.y = element_text(size = 9),plot.title=element_text(size=12),
        plot.margin = margin(10, 145, 3, 3))

# Supplementary Table 11: all eligible GSE227809 one-week clocks, lowest MAD first.
# MAD = mean(abs(predicted-age change - elapsed years)); RC95 = 1.96*SD(signed error).
# ME = mean(predicted age - recorded chronological age) across both visits.
# All metrics are in years; seven days is approximate and A first is user-directed.
supplementary_table_11 <- h_all_summary %>%
  transmute(clock, MAD=mean_absolute_deviation, RC95=reliability_coefficient_95,
            ME=mean_error)
stopifnot(!anyDuplicated(supplementary_table_11$clock),
          all(diff(supplementary_table_11$MAD)>=0))
write_csv(supplementary_table_11, 'paper/output/Supplementary_Table_11.csv', na='')

#### Panel I ####
response_saved <- read_csv('paper/input/benchmark/clock_sensitivity_ranking.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
# Select the ten highest and lowest primary responsiveness scores.
i_response_order <- response_saved %>%
  filter(clock %in% clocks, is.finite(responsiveness_score))
i_selected <- bind_rows(
  i_response_order %>% arrange(desc(responsiveness_score), clock) %>% slice_head(n = 10) %>% mutate(response_group = 'Responsive'),
  i_response_order %>% arrange(responsiveness_score, clock) %>% slice_head(n = 10) %>% mutate(response_group = 'Non-responsive'))

# Retain independent contrasts with at least four samples per group and finite positive pooled SD.
# Reuse panel B's metadata-derived residuals.
i_examples <- map2_dfr(datasets, b_adjusted, function(d, adjusted) {
  bind_cols(d %>% select(GSE_ID, patient_id, treatment, case_control), as_tibble(adjusted$residuals))
}) %>%
  filter(!is.na(treatment), case_control %in% c('case', 'ctrl')) %>%
  group_by(GSE_ID, treatment) %>%
  filter(!anyDuplicated(patient_id[!is.na(patient_id)]), sum(case_control == 'case') >= 4, sum(case_control == 'ctrl') >= 4) %>%
  ungroup() %>%
  select(GSE_ID, treatment, case_control, all_of(i_selected$clock)) %>%
  pivot_longer(all_of(i_selected$clock), names_to = 'clock', values_to = 'acceleration') %>%
  filter(is.finite(acceleration)) %>%
  group_by(GSE_ID, treatment, clock, case_control) %>%
  summarise(n = n(), mean = mean(acceleration), variance = var(acceleration), .groups = 'drop') %>%
  pivot_wider(names_from = case_control, values_from = c(n, mean, variance)) %>%
  filter(n_case >= 4, n_ctrl >= 4) %>%
  mutate(pooled_sd = sqrt(((n_case - 1) * variance_case + (n_ctrl - 1) * variance_ctrl) / (n_case + n_ctrl - 2)),
         acceleration_difference = mean_case - mean_ctrl) %>%
  filter(is.finite(pooled_sd), pooled_sd > 0, is.finite(acceleration_difference)) %>%
  left_join(i_selected %>% select(clock, response_group), by = 'clock')

# Three displayed examples per category; labels and category assignments are
# explicit display mappings. Selection uses sample coverage, not effect size.
i_categories <- tribble(
  ~treatment, ~category,
  'Psychosocial stress exposure', 'Exposure',
  'Sexual/physical abuse exposure', 'Exposure',
  'Early life adversity', 'Exposure',
  'Prenatal alcohol exposure', 'Exposure',
  'In utero famine exposure', 'Exposure',
  'Environmental chromium exposure', 'Exposure',
  'Methylmercury exposure', 'Exposure',
  'Intrauterine Zika virus exposure', 'Exposure',
  'Cancer therapy', 'Treatment',
  'Levodopa treatment', 'Treatment',
  'Chemotherapy and radiotherapy', 'Treatment',
  'Radiotherapy', 'Treatment',
  'Vitamin K supplementation', 'Treatment',
  'Leflunomide', 'Treatment',
  'Transcranial magnetic stimulation', 'Treatment',
  'Doxorubicin-based chemotherapy', 'Treatment',
  'DCCT conventional versus intensive diabetes therapy', 'Treatment',
  'Cigarette smoking exposure', 'Lifestyle',
  'Athletic participation', 'Lifestyle',
  'Aerobic/endurance exercise training', 'Lifestyle')
i_show <- i_examples %>%
  group_by(GSE_ID, treatment) %>%
  summarise(clocks = n_distinct(clock), response_groups = n_distinct(response_group),
            n = min(n_case + n_ctrl), .groups = 'drop') %>%
  left_join(i_categories, by = 'treatment') %>%
  filter(response_groups == 2, !is.na(category)) %>%
  arrange(desc(n), GSE_ID) %>%
  distinct(treatment, .keep_all = TRUE) %>%
  group_by(category) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(category = factor(category, c('Exposure', 'Treatment', 'Lifestyle'))) %>%
  arrange(category, desc(n))
i_plot_points <- i_examples %>%
  inner_join(i_show %>% select(GSE_ID, treatment, category),
             by = c('GSE_ID', 'treatment')) %>%
  mutate(label = factor(treatment, levels = rev(i_show$treatment)), response_group = factor(response_group, c('Responsive', 'Non-responsive')))
i_plot_data <- i_plot_points %>%
  group_by(GSE_ID, treatment, category, label, response_group) %>%
  summarise(acceleration_difference = mean(acceleration_difference), .groups = 'drop')

I <- ggplot(i_plot_data, aes(acceleration_difference, label, color = response_group, shape = response_group)) +
  theme_pubr(base_size = 13, border = TRUE) +
  theme(plot.margin = margin(3, 4, 3, 3), legend.title = element_blank()) +
  geom_vline(xintercept = 0, color = colors$grey[3], linetype = 2) +
  geom_point(data = i_plot_points, aes(group = response_group), position = position_jitterdodge( jitter.width = 0, jitter.height = .12,
               dodge.width = .55, seed = 42),
             color = colors$grey[4], alpha = .65, size = 1.5) +
  geom_point(position = position_dodge(width = .55), size = 2) +
  facet_grid(category ~ ., scales = 'free_y', space = 'free_y') +
  scale_y_discrete(labels = function(x) str_wrap(x, 25)) +
  scale_color_manual( values = c(Responsive = colors$orange[4], 'Non-responsive' = colors$teal[4]), labels = c(Responsive = 'Responsive clocks',
               'Non-responsive' = 'Non-responsive clocks')) +
  scale_shape_manual(values = c(Responsive = 16, 'Non-responsive' = 17),
    labels = c(Responsive = 'Responsive clocks', 'Non-responsive' = 'Non-responsive clocks')) +
  labs(x = 'Mean age-acceleration difference (case - control, years)', y = NULL) +
  theme(axis.text.y = element_text(size = 9), strip.text.y = element_text(size = 10), legend.position = 'top', legend.text = element_text(size = 9)) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme()
I

# Supplementary Table 12: primary responsiveness score for every selected clock.
# Mean direction-aligned adjusted age differences, median reference sign and equal comparison weights.
# Units: years; higher scores rank first. Retain missing scores as blank.
supplementary_table_12 <- tibble(clock=clocks) %>%
  left_join(response_saved %>% select(clock, responsiveness_score), by='clock') %>%
  arrange(desc(responsiveness_score), clock) %>%
  mutate(responsiveness_rank=if_else(is.finite(responsiveness_score), row_number(), NA_integer_)) %>%
  select(responsiveness_rank, clock, responsiveness_score)
stopifnot(nrow(supplementary_table_12)==length(clocks),
          !anyDuplicated(supplementary_table_12$clock))
write_csv(supplementary_table_12, 'paper/output/Supplementary_Table_12.csv', na='')

#### Panel J ####
summary_saved <- read_csv('paper/input/benchmark/clock_summary.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
mortality_saved <- read_csv('paper/input/benchmark/mortality_summary.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
disease_saved <- read_csv('paper/input/benchmark/disease_auc_summary.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
general_saved <- read_csv('paper/input/benchmark/clock_generalizability_summary.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
reliability_saved <- read_csv('paper/input/benchmark/overall_clock_reliability_metrics.csv', na = c('', 'NA', 'NaN'), show_col_types = FALSE, progress = FALSE)
# Compare benchmark ranks: rank 1 is preferable, including least cell influence.
# Spearman rho is calculated on benchmark ranks with pairwise complete clocks.
j_values <- list(summary_saved %>% select(clock, `Age prediction`=age_rank),
  disease_saved %>% select(clock, Disease=disease_rank),
  mortality_saved %>% select(clock, Mortality=mortality_rank),
  reliability_saved %>% select(clock, Reliability=reliability_rank),
  general_saved %>% select(clock, Generalizability=generalizability_rank),
  response_saved %>% select(clock, Responsiveness=responsiveness_rank),
  b_data %>% transmute(clock, `Cell composition influence`=n()+1L-rank)
) %>% reduce(left_join, by='clock')
j_cor <- j_values %>% select(-clock) %>% as.matrix() %>%
  cor(method='spearman', use='pairwise.complete.obs')
j_labels <- c('Age prediction', 'Disease\ndiscrimination', 'Mortality\nassociation',
              'Within-person\nreliability', 'Geographic\ngeneralizability',
              'Responsiveness\nto interventions', 'Cell-composition\ninfluence')
j_cluster_matrix <- j_cor
j_cluster_matrix[!is.finite(j_cluster_matrix)] <- 0
j_cluster <- hclust(as.dist(pmax(1 - j_cluster_matrix, 0)), method = 'average')
j_heatmap <- Heatmap(j_cor, name = 'Spearman rho', col = cor_colors, na_col = colors$grey[2],
  width = grid::unit(5.7, 'in'), height = grid::unit(5.2, 'in'),
  cluster_rows = as.dendrogram(j_cluster),
  cluster_columns = as.dendrogram(j_cluster),
  row_dend_side = 'right', column_dend_side = 'top',
  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_labels = j_labels, column_labels = j_labels,
  row_names_side = 'left', column_names_side = 'bottom', column_names_rot = 90,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  rect_gp = grid::gpar(type = 'none'),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill, col = '#d9dde3', lwd = .5))
    grid::grid.text(sprintf('%.2f', j_cor[i, j]), x, y, gp = grid::gpar(fontsize = 11, col = 'black'))
  },
  heatmap_legend_param = list(direction = 'horizontal', title_position = 'topcenter',
    legend_width = grid::unit(1.8, 'in'), title_gp = grid::gpar(fontsize = 10),
    labels_gp = grid::gpar(fontsize = 9)))
# Capture displayed J at its final square panel size so fixed-size heatmap
# cells occupy the panel rather than sitting in a smaller nested canvas.
j_grob <- local({
  previous_padding <- ComplexHeatmap::ht_opt('HEATMAP_LEGEND_PADDING')
  on.exit(ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = previous_padding))
  ComplexHeatmap::ht_opt(HEATMAP_LEGEND_PADDING = grid::unit(1, 'mm'))
  grid::grid.grabExpr(draw(j_heatmap, heatmap_legend_side = 'top',
    padding = grid::unit(rep(1, 4), 'mm')),
    width = square_side, height = square_side, wrap.grobs = TRUE)
})
J <- ggdraw(j_grob) +
  theme_pubr() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank())
J

# Supplementary Table 13: the seven domain ranks used in panel J.
# Rank 1 is preferable; cell composition rank 1 means the least influence.
# Mean rank is the unweighted arithmetic mean of all seven domain ranks.
# Overall rank orders those means ascending; equal means share the minimum rank.
# Lower values indicate a lower mean rank across domains.
# Keep both summary columns blank if any domain rank is unavailable.
# These summary columns do not enter panel J's correlations or panel K's group comparison.
supplementary_table_13 <- j_values %>%
  mutate(`Mean rank` = rowMeans(across(-clock), na.rm = FALSE),
         `Overall rank` = min_rank(`Mean rank`)) %>%
  relocate(`Overall rank`, `Mean rank`, .after = clock) %>%
  arrange(`Mean rank`, clock)
stopifnot(nrow(supplementary_table_13)==length(clocks),
          !anyDuplicated(supplementary_table_13$clock))
write_csv(supplementary_table_13, 'paper/output/Supplementary_Table_13.csv', na='')

#### Supplementary Figure 5 ####
# Clock-level reliability deviation versus responsiveness score.
k_data <- reliability_saved %>%
  select(clock, reliability_deviation = mean_absolute_deviation) %>%
  left_join(response_saved %>% select(clock, responsiveness_score), by = 'clock') %>%
  filter(is.finite(reliability_deviation), is.finite(responsiveness_score))
supplementary_figure_5 <- ggplot(k_data, aes(reliability_deviation, responsiveness_score)) +
  theme_pubr(base_size = 13) +
  theme(plot.margin = margin(3, 4, 3, 3), legend.title = element_blank()) +
  geom_smooth(method = 'lm', formula = y ~ x, color = colors$grey[6], fill = colors$red[1], linewidth = .55) +
  geom_point(color = colors$red[3], alpha = .8, size = 2) +
  ggrepel::geom_text_repel(data = bind_rows(
      k_data %>% slice_min(reliability_deviation, n = 5, with_ties = FALSE),
      k_data %>% slice_max(reliability_deviation, n = 5, with_ties = FALSE),
      k_data %>% slice_min(responsiveness_score, n = 5, with_ties = FALSE),
      k_data %>% slice_max(responsiveness_score, n = 5, with_ties = FALSE)) %>%
      distinct(clock, .keep_all = TRUE), aes(label = clock),
    color = colors$grey[6], size = 3.1, box.padding = .25,
    min.segment.length = 0, max.overlaps = Inf, seed = 42) +
  labs(x = 'Reliability MAD (years; lower is better)', y = 'Responsiveness score (years; higher is better)') +
  theme(aspect.ratio = NULL)
supplementary_figure_5

#### Panel K ####
# Group comparison from Supplementary Tables 1 and 13 only.
# Original 100-clock ranks, with rank 1 meaning least cell influence.
k_domains <- c('Age prediction', 'Generalizability', 'Reliability',
  'Cell composition influence', 'Disease', 'Mortality', 'Responsiveness')
k_annotations <- read_csv('paper/output/Supplementary_Table_1.csv', show_col_types = FALSE) %>%
  transmute(clock = `Benchmark / platform clock name`, Outcome)
k_ranks <- read_csv('paper/output/Supplementary_Table_13.csv', show_col_types = FALSE)
stopifnot(!anyDuplicated(k_annotations$clock), !anyDuplicated(k_ranks$clock),
          setequal(k_annotations$clock, k_ranks$clock))
k_groups <- c('Chronological age', 'Mortality-related age')
k_clock_points <- k_ranks %>% left_join(k_annotations, by = 'clock') %>%
  filter(Outcome %in% k_groups) %>%
  pivot_longer(all_of(k_domains), names_to = 'domain', values_to = 'rank') %>%
  mutate(y = 8 - match(domain, k_domains) + if_else(Outcome == k_groups[1], .12, -.12),
         Outcome = factor(Outcome, levels = k_groups))
k_plot_data <- k_clock_points %>% group_by(Outcome, domain) %>%
  summarise(n = n(), median = median(rank), q25 = quantile(rank, .25),
            q75 = quantile(rank, .75), .groups = 'drop') %>%
  mutate(y = 8 - match(domain, k_domains) + if_else(Outcome == k_groups[1], .12, -.12),
         Outcome = factor(Outcome, levels = k_groups))
stopifnot(nrow(k_plot_data) == 14L,
          all(k_plot_data$n[k_plot_data$Outcome == k_groups[1]] == 62L),
          all(k_plot_data$n[k_plot_data$Outcome == k_groups[2]] == 26L),
          all(is.finite(k_plot_data$median)))
K <- ggplot(k_plot_data, aes(color = Outcome, shape = Outcome)) +
  geom_vline(xintercept = c(1, 25, 50, 75, 100), color = '#e7e7e7', linewidth = .35) +
  geom_hline(yintercept = 3.5, linetype = 'dashed', color = '#aaaaaa', linewidth = .4) +
  geom_point(data = k_clock_points, aes(x = rank, y = y),
    position = position_jitter(width = 0, height = .065, seed = 20261016),
    size = 1.25, alpha = .35, show.legend = FALSE) +
  geom_segment(aes(x = q25, xend = q75, y = y, yend = y), linewidth = 1) +
  geom_point(aes(x = median, y = y), size = 2.8) +
  scale_color_manual(values = c('#2374AB', '#C65A28'),
    labels = c('Age predictors (n = 62)', 'Mortality-related (n = 26)')) +
  scale_shape_manual(values = c(16, 17),
    labels = c('Age predictors (n = 62)', 'Mortality-related (n = 26)')) +
  scale_x_continuous(limits = c(1, 100), breaks = c(1, 25, 50, 75, 100),
    expand = expansion(mult = c(.02, .025))) +
  scale_y_continuous(limits = c(.5, 7.5), breaks = 7:1,
    labels = c('Age prediction', 'Geographic\ngeneralizability',
      'Within-person\nreliability', 'Cell-composition\ninfluence',
      'Disease\ndiscrimination', 'Mortality\nassociation', 'Responsiveness\nto interventions'),
    expand = expansion(mult = 0)) +
  labs(x = 'Benchmark rank', y = NULL, color = NULL, shape = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = 'top', legend.direction = 'vertical',
    legend.text = element_text(size = 11),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(18, 9, 8, 5)) +
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
K
# Figure assembly and display -------------------------------------------------
# Plot objects and section headings use the displayed panel letters.
# Top: A occupies half the width, B above age curves C in the other half.
# Lower rows: D | E (two columns); F | G | H; I | J | K.
A <- A + theme(aspect.ratio = 1, plot.margin = margin(0, 0, 0, 0))
C <- C + theme(aspect.ratio = NULL)
D <- D + theme(aspect.ratio = NULL)
F <- F + theme(aspect.ratio = 1)
G <- G + theme(aspect.ratio = NULL)
H <- H + theme(aspect.ratio = NULL)
J <- J + theme(aspect.ratio = 1, plot.margin = margin(0, 0, 0, 0))
K <- K + theme(aspect.ratio = NULL)
# Give C a wide plotting area like B and remove excess height around A.
first_row_height <- 9.5
b_height <- first_row_height / 2
figure_height <- first_row_height + second_row_height + 2 * square_side
right_top <- plot_grid(B, C, ncol = 1,
  rel_heights = c(b_height, first_row_height - b_height),
  labels = c('B', 'C'), label_size = 18, label_x = .012, label_y = .99) +
  theme_pubr() + theme_void()
row1 <- plot_grid(A, right_top, ncol = 2, rel_widths = c(9.5, 12.5),
  labels = c('A', ''), label_size = 18, label_x = .012, label_y = .99) +
  theme_pubr() + theme_void()
row2 <- plot_grid(D, E, ncol = 2, rel_widths = c(1, 2),
  labels = c('D', 'E'), label_size = 18, label_x = .012, label_y = .99) +
  theme_pubr() + theme_void()
row3 <- plot_grid(F, G, H, ncol = 3, rel_widths = c(1, 1, 1),
  labels = c('F', 'G', 'H'), label_size = 18, label_x = .012, label_y = .99) +
  theme_pubr() + theme_void()
row4 <- plot_grid(I, J, K, ncol = 3, rel_widths = c(1, 1, 1),
  labels = c('I', 'J', 'K'), label_size = 18, label_x = .012, label_y = .99) +
  theme_pubr() + theme_void()
combined <- plot_grid(row1, row2, row3, row4, ncol = 1,
  rel_heights = c(first_row_height, second_row_height, square_side, square_side)) +
  theme_pubr() + theme_void()
figure <- ggdraw(combined) + theme_pubr() + theme_void()

invisible(dev.off())
# Display interactively; avoid Rscript creating an automatic Rplots.pdf.
if (interactive()) print(figure)

# Physical heights keep the lower two rows square while D/E is more compact.
pdf('paper/output/Figure_2.pdf', width = figure_width, height = figure_height)
print(figure)
dev.off()
# Export the same assembled figure at the PDF's physical size, with a white background.
ragg::agg_png('paper/output/Figure_2.png', width = figure_width, height = figure_height,
    units = 'in', res = 300, background = 'white')
print(figure)
dev.off()
message('Saved Figure 2 PDF and PNG.')

pdf('paper/output/Supplementary_Figure_5.pdf', width = square_side, height = square_side)
print(supplementary_figure_5 + theme(aspect.ratio = 1))
dev.off()
message('Saved Supplementary Figure 5 PDF.')

# Supplementary Table 14: required CpG coverage before/after reference filling.
# This block can be run independently without generating figures. It reads the
# current Stage 06 selection, rather than the historical paper/input snapshot.
# Percentages use the required CpGs as denominator. After-fill coverage includes
# reference-imputed probes, but does not count residual zero-filled probes as
# measured/reference-covered. Retained means at least one finite clock estimate
# in the integrated dataset, not eligibility for every individual benchmark.
# If retained samples use multiple source accessions, report the minimum coverage
# among the calculation sources actually contributing finite estimates.
coverage_audit <- readr::read_csv(
  'output/06_data_integration/clock_pair_eligibility.csv',
  show_col_types = FALSE, progress = FALSE)
coverage_selection <- readr::read_csv(
  'output/06_data_integration/integrated_metadata.csv',
  col_types = readr::cols(.default = readr::col_character()),
  na = c('', 'NA', 'NaN'), progress = FALSE)
coverage_sources <- readr::read_csv(
  'output/06_data_integration/parent_sample_provenance.csv',
  col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
coverage_clocks <- readr::read_csv(
  'input/06_data_integration/clock_names/selected_clock_catalogue.csv',
  show_col_types = FALSE, progress = FALSE)$benchmark_name %>% sort()
stopifnot(!anyDuplicated(coverage_selection$GSM_ID),
          !anyDuplicated(coverage_sources$GSM_ID),
          !anyDuplicated(coverage_audit[c('GSE_ID', 'clock')]),
          all(coverage_clocks %in% names(coverage_selection)),
          all(coverage_selection$GSM_ID %in% coverage_sources$GSM_ID))
coverage_selected_sources <- coverage_selection %>%
  select(GSE_ID, GSM_ID) %>%
  left_join(coverage_sources %>% select(GSM_ID, estimates_source_dataset_id),
            by = 'GSM_ID')
coverage_used_sources <- bind_rows(lapply(coverage_clocks, function(clock_name) {
  finite <- is.finite(suppressWarnings(as.numeric(coverage_selection[[clock_name]])))
  coverage_selected_sources[finite, ] %>%
    distinct(GSE_ID, estimates_source_dataset_id) %>%
    mutate(clock = clock_name)
}))
coverage_values <- coverage_audit %>%
  transmute(estimates_source_dataset_id = GSE_ID, clock,
    before = 100 * (required_cpgs - absent_before_reference) / required_cpgs,
    after = 100 * (required_cpgs - missing_cpgs) / required_cpgs,
    eligible, clock_id,
    coverage_pass = if_else(clock_id == 'CLK002',
      5 * missing_cpgs <= required_cpgs,
      10 * missing_cpgs < required_cpgs))
stopifnot(all(is.finite(coverage_values$before)),
          all(is.finite(coverage_values$after)),
          all(coverage_values$before >= 0),
          all(coverage_values$after <= 100),
          all(coverage_values$before <= coverage_values$after),
          all(coverage_audit$reference_imputed_cpgs[
            coverage_audit$clock_id %in% c('CLK137', 'CLK152')] == 0),
          all(coverage_audit$absent_before_reference -
              coverage_audit$reference_imputed_cpgs == coverage_audit$missing_cpgs))
coverage_contributors <- coverage_used_sources %>%
  left_join(coverage_values, by = c('estimates_source_dataset_id', 'clock'))
stopifnot(all(is.finite(coverage_contributors$before)),
          all(is.finite(coverage_contributors$after)),
          all(coverage_contributors$eligible),
          all(coverage_contributors$coverage_pass))
coverage_retained <- coverage_contributors %>%
  group_by(GSE_ID, clock) %>%
  summarise(before_retained = min(before), after_retained = min(after),
            .groups = 'drop')
supplementary_table_14 <- tidyr::expand_grid(
  GSE_ID = sort(unique(coverage_selection$GSE_ID)), clock = coverage_clocks) %>%
  left_join(coverage_values %>% rename(GSE_ID = estimates_source_dataset_id),
            by = c('GSE_ID', 'clock')) %>%
  left_join(coverage_retained, by = c('GSE_ID', 'clock')) %>%
  transmute(dataset = GSE_ID, clock,
    `CpG coverage before filling (%)` = coalesce(before_retained, before),
    `CpG coverage after reference filling (%)` = coalesce(after_retained, after),
    status = if_else(is.finite(before_retained), 'Retained', 'Excluded'))
stopifnot(nrow(supplementary_table_14) ==
            length(unique(coverage_selection$GSE_ID)) * length(coverage_clocks),
          !anyDuplicated(supplementary_table_14[c('dataset', 'clock')]),
          all(is.finite(supplementary_table_14[[3]])),
          all(is.finite(supplementary_table_14[[4]])))
readr::write_csv(supplementary_table_14,
  'paper/output/Supplementary_Table_14.csv', na = '')
message('Saved Supplementary Table 14: pre/post-fill CpG coverage and selection status.')

# Supplementary Figure 6: coverage versus primary score in all seven domains.
# One point per clock. Coverage counts the units used in the primary aggregation:
# datasets, disease/responsiveness comparisons, or generalizability countries.
# evaluable_studies in the responsiveness export counts comparisons, not datasets.
coverage_plot_specs <- list(
  list(domain = 'Age', file = 'clock_summary.csv', coverage = 'age_datasets',
       score = 'median_age_r_squared', unit = 'datasets',
       label = 'Median squared Pearson correlation'),
  list(domain = 'Disease', file = 'disease_auc_summary.csv', coverage = 'comparisons',
       score = 'median_auc', unit = 'comparisons', label = 'Median AUC'),
  list(domain = 'Mortality', file = 'mortality_summary.csv', coverage = 'datasets',
       score = 'weighted_median_hazard_ratio', unit = 'datasets',
       label = 'Sample-weighted median hazard ratio'),
  list(domain = 'Cell influence', file = 'cell_influence_summary.csv',
       coverage = 'cell_influence_datasets', score = 'cell_influence_years',
       unit = 'datasets', label = 'Mean absolute residual difference (years)'),
  list(domain = 'Reliability', file = 'overall_clock_reliability_metrics.csv',
       coverage = 'n_datasets', score = 'mean_absolute_difference',
       unit = 'datasets', label = 'Mean dataset MAD (years)'),
  list(domain = 'Responsiveness', file = 'clock_sensitivity_ranking.csv',
       coverage = 'evaluable_studies', score = 'responsiveness_score',
       unit = 'comparisons', label = 'Responsiveness score (years)'),
  list(domain = 'Generalizability', file = 'clock_generalizability_summary.csv',
       coverage = 'median_country_mae_n_countries', score = 'median_country_median_mae',
       unit = 'countries', label = 'Median country median MAE (years)'))
coverage_plots <- list()
coverage_correlation_rows <- list()
for (coverage_spec in coverage_plot_specs) {
  coverage_source <- readr::read_csv(
    file.path('paper/input/benchmark', coverage_spec$file),
    show_col_types = FALSE, progress = FALSE)
  coverage_data <- data.frame(
    clock = coverage_source$clock,
    coverage = coverage_source[[coverage_spec$coverage]],
    score = coverage_source[[coverage_spec$score]])
  stopifnot(nrow(coverage_data) == 100L,
            !anyDuplicated(coverage_data$clock),
            all(is.finite(coverage_data$coverage)),
            all(is.finite(coverage_data$score)),
            all(coverage_data$coverage > 0),
            all(coverage_data$coverage == floor(coverage_data$coverage)))
  coverage_plot <- ggplot2::ggplot(coverage_data, ggplot2::aes(coverage, score)) +
    ggplot2::geom_point(size = 2.8, alpha = .65, colour = '#5495ce') +
    ggplot2::labs(x = paste('Number of valid', coverage_spec$unit),
                 y = paste(coverage_spec$domain, coverage_spec$label, sep = '\n')) +
    ggplot2::theme_classic(base_size = 13)
  if (coverage_spec$domain == 'Mortality') {
    coverage_plot <- coverage_plot + ggplot2::scale_x_continuous(
      breaks = sort(unique(coverage_data$coverage)))
  }
  if (coverage_spec$domain %in% c('Reliability', 'Generalizability')) {
    coverage_plot <- coverage_plot + ggplot2::scale_x_continuous(
      breaks = seq(min(coverage_data$coverage), max(coverage_data$coverage), 2))
  }
  if (coverage_spec$domain == 'Responsiveness') {
    coverage_plot <- coverage_plot + ggplot2::scale_x_continuous(
      breaks = seq(25, 65, 5))
  }
  coverage_plots[[coverage_spec$domain]] <- coverage_plot + ggplot2::theme(
    axis.title = ggplot2::element_text(size = 10),
    axis.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(8, 10, 8, 8))
  # Exploratory two-sided Pearson tests: related clocks are not fully independent.
  coverage_test <- stats::cor.test(coverage_data$coverage, coverage_data$score,
                                  method = 'pearson', alternative = 'two.sided')
  coverage_correlation_rows[[coverage_spec$domain]] <- data.frame(
    domain = coverage_spec$domain,
    coverage_unit = coverage_spec$unit,
    coverage_field = coverage_spec$coverage,
    score_field = coverage_spec$score,
    n_clocks = nrow(coverage_data),
    min_coverage = min(coverage_data$coverage),
    max_coverage = max(coverage_data$coverage),
    pearson_r = unname(coverage_test$estimate),
    ci_low = coverage_test$conf.int[1],
    ci_high = coverage_test$conf.int[2],
    p_value = coverage_test$p.value)
}
coverage_correlations <- do.call(rbind, coverage_correlation_rows)
coverage_correlations$q_value_BH <- stats::p.adjust(
  coverage_correlations$p_value, method = 'BH')
coverage_correlations$p_value_Holm <- stats::p.adjust(
  coverage_correlations$p_value, method = 'holm')
supplementary_figure_6 <- cowplot::plot_grid(
  plotlist = c(coverage_plots, list(NULL)), ncol = 2, align = 'none')
pdf('paper/output/Supplementary_Figure_6.pdf', width = 12, height = 15)
print(supplementary_figure_6)
dev.off()
message('Saved Supplementary Figure 6: seven-domain coverage plots.')

# Supplementary Table 15: cell-reference marker coverage before imputation.
# Coverage counts probe rows present, not sample-level completeness of beta values.
# Include integrated datasets retained under the benchmark dataset-selection policy.
cell_marker_processing <- readr::read_csv(
  'output/05_cell_composition/no_external_imputation_20261007/processing_summary.csv',
  show_col_types = FALSE)
cell_reference_environment <- new.env()
utils::data('centDHSbloodDMC.m', package = 'EpiDISH',
            envir = cell_reference_environment)
cell_reference_marker_count <- nrow(cell_reference_environment$centDHSbloodDMC.m)
stopifnot(cell_reference_marker_count == 333L,
          !anyDuplicated(cell_marker_processing$dataset))
supplementary_table_15 <- tibble::tibble(
  dataset = sort(unique(coverage_selection$GSE_ID))) %>%
  left_join(cell_marker_processing %>% select(dataset, cohort_probes), by = 'dataset') %>%
  transmute(dataset,
    `Cell-reference CpGs` = cell_reference_marker_count,
    `CpGs present before imputation` = cohort_probes,
    `Cell-reference coverage before imputation (%)` =
      100 * cohort_probes / cell_reference_marker_count)
stopifnot(nrow(supplementary_table_15) == length(unique(coverage_selection$GSE_ID)),
          !anyNA(supplementary_table_15[[3]]),
          all(supplementary_table_15[[3]] >= 0),
          all(supplementary_table_15[[3]] <= cell_reference_marker_count))
readr::write_csv(supplementary_table_15,
  'paper/output/Supplementary_Table_15.csv', na = '')
message('Saved Supplementary Table 15: cell-reference marker coverage before imputation.')
