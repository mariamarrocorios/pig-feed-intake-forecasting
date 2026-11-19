# Figures of the paper and the supplementary material.
# Nothing is refitted here: the figures are drawn from the saved results, so
# this can be run on its own after 02_forecasting_study.R to change how a
# figure looks. When 02 sources it, the objects already in memory are used.
#
# Figure numbers in the manuscript:
#   Fig1 -> Figure 1     Fig3 -> Figure 2     Fig7 -> Figure 3
#   Fig8 -> Figure 4     Fig5 -> Figure 5     Fig9 -> Figure 6
#   Fig4 -> Figure S1    Fig6 -> Figure S2    Fig10 -> Figure S3
#   FigS1 -> Figure S4   (Fig2 is not used)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

if (!exists("out_path")) out_path <- "outputs/results"
if (!exists("theme_paper")) source("R/functions.R")
if (!exists("out_tables")) out_tables <- file.path(out_path, "tables")

if (!dir.exists(out_tables)) {
  stop("Not found: ", out_tables, ". Run 02_forecasting_study.R first.")
}

# use the object in memory if it exists, otherwise read it from disk
need <- function(name, loader) {
  if (exists(name, envir = globalenv())) return(get(name, envir = globalenv()))
  loader()
}
tbl <- function(f) read_csv(file.path(out_tables, f), show_col_types = FALSE)

relabel_fset <- function(d) {
  d %>% mutate(fset_label = ifelse(fset == "none", fset_label_none,
                                   unname(fset_labels[fset])))
}

metrics_main   <- need("metrics_main",   function() tbl("T1_metrics_by_fold.csv")) %>% relabel_fset()
summary_main   <- need("summary_main",   function() tbl("T2_metrics_summary.csv")) %>% relabel_fset()
added_value    <- need("added_value",    function() tbl("T4_added_value_vs_feed_history.csv"))
bottom_up      <- need("bottom_up",      function() tbl("T5_direct_vs_bottom_up.csv"))
by_phase       <- need("by_phase",       function() tbl("T8_error_by_phase_pen.csv"))
perm_imp       <- need("perm_imp",       function() tbl("T12_grouped_permutation_importance.csv"))
ale_tbl        <- need("ale_tbl",        function() tbl("T13_ale_pen_full_model.csv"))
ale_range      <- need("ale_range",      function() tbl("T13b_ale_range_pen_full_model.csv"))
learning_curve <- need("learning_curve", function() tbl("S4_learning_curve_incl_linear.csv"))
lm_vs_rf       <- need("lm_vs_rf",       function() tbl("S2_linear_vs_rf_by_fold.csv"))

all_preds <- need("all_preds", function() {
  f <- file.path(out_path, "predictions", "predictions_pen_batch.csv")
  if (!file.exists(f)) {
    warning("Not found: ", f, " - Fig9 is skipped.", call. = FALSE)
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE)
})

ale_rug <- need("ale_rug", function() {
  f <- file.path(out_tables, "T13c_ale_rug_sample.csv")
  if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
})

horizons <- if (exists("horizons")) horizons else sort(unique(summary_main$horizon))
h1       <- min(horizons)
if (!exists("final_fset")) final_fset <- "F+B"

# predictor sets actually present in part A (without climate data part A
# uses the same sets as part B)
setsA <- intersect(c(fsetA_order, "F+G+B"), unique(summary_main$fset[summary_main$part == "A"]))
fsetA_levels <- rev(unname(fset_labels[setsA]))

message("Drawing figures into ", out_figures, " (", fig_dpi, " dpi)")

# ============================================================
# Figures
# ============================================================

# Fig 1 - RQ1 x RQ3: nRMSE by scale and horizon (feed history)
fig1 <- summary_main %>%
  filter(part == "A", fset %in% c("none", "F")) %>%
  prep_scale() %>% prep_model() %>%
  ggplot(aes(horizon, nrmse_mean, colour = model_label, shape = model_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_x_continuous(breaks = horizons) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = "Forecast horizon (days)", y = "Normalised RMSE (% of mean intake)") +
  theme_paper
save_fig(fig1, "Fig1_nrmse_by_scale_horizon.png", 9, 4.2)

# Fig 2 - RQ2: skill vs persistence by predictor set (first horizon)
# not in the paper (the kg model goes down to about -250% at batch level)
fig2 <- summary_main %>%
  filter(part == "A", horizon == h1, model %in% c("rf_level", "rf_rel")) %>%
  prep_scale() %>% prep_model() %>%
  mutate(fset_label = factor(fset_label, levels = fsetA_levels)) %>%
  ggplot(aes(100 * skill_mean, fset_label, colour = model_label, shape = model_label)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = 100 * skill_min, xmax = 100 * skill_max),
                width = 0, linewidth = 0.4, orientation = "y",
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ scale_lab, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = "RMSE reduction vs persistence (%) - mean and range over test batches",
       y = NULL) +
  theme_paper
save_fig(fig2, "Fig2_skill_by_predictor_set.png", 10, 4.8)

# Fig 3 -> manuscript Figure 2
# RQ2 x RQ3: change in RMSE when predictors are added (RF relative target)
fig3 <- added_value %>%
  filter(part == "A", model == "rf_rel") %>%
  prep_scale() %>%
  mutate(fset_label = factor(unname(fset_labels[fset]), levels = fsetA_levels),
         horizon_lab = paste0("h = ", horizon, " d")) %>%
  ggplot(aes(pct_change_vs_F, fset_label, shape = test_batch)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(size = 2.2, colour = model_colours[["RF (relative target)"]],
             position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  facet_grid(horizon_lab ~ scale_lab, scales = "free_x") +
  scale_shape_manual(values = batch_shapes, labels = batch_labels, name = "Test batch") +
  labs(x = "Change in RMSE vs feed history only (%)  (negative = better)", y = NULL) +
  theme_paper
save_fig(fig3, "Fig3_added_value_of_predictors.png", 10, 2.6 + 2 * length(horizons))

# Fig 4 -> supplementary Figure S1 : direct vs bottom-up
fig4_data <- bottom_up %>%
  filter(part == "A", model %in% c("persistence", "mean5", "rf_level", "rf_rel"),
         fset %in% c("none", "F"))
if (nrow(fig4_data) > 0) {
  fig4 <- fig4_data %>%
    prep_model() %>%
    mutate(target_lab = factor(ifelse(target == "pen", "Pen target (from animal forecasts)",
                                      "Batch target (from pen forecasts)"),
                               levels = c("Pen target (from animal forecasts)",
                                          "Batch target (from pen forecasts)"))) %>%
    ggplot(aes(pct_diff_bottom_up, model_label, colour = model_label, shape = test_batch)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_point(size = 2.2, position = position_jitter(height = 0.15, width = 0, seed = 1)) +
    facet_grid(paste0("h = ", horizon, " d") ~ target_lab) +
    scale_colour_manual(values = model_colours, guide = "none") +
    scale_shape_manual(values = batch_shapes, labels = batch_labels, name = "Test batch") +
    labs(x = "RMSE of bottom-up minus direct forecast (%)  (positive = direct is better)",
         y = NULL) +
    theme_paper
  save_fig(fig4, "Fig4_direct_vs_bottom_up.png", 10, 2.2 + 1.8 * length(horizons))
}

# Fig 5 -> manuscript Figure 5 : error by fattening phase (pen, first horizon)
fig5 <- by_phase %>%
  prep_model() %>%
  mutate(test_lab = batch_labels[test_batch]) %>%
  ggplot(aes(phase, rmse, colour = model_label, shape = model_label, group = model_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ test_lab, nrow = 1) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = NULL, y = "RMSE (kg/animal/day)") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(fig5, "Fig5_error_by_phase_pen.png", 11, 4.2)

# Fig 6 -> supplementary Figure S2 : external test on batch 4
fig6 <- metrics_main %>%
  filter(part == "B") %>%
  prep_scale() %>% prep_model() %>%
  mutate(fset_label = factor(fset_label, levels = fsetB_levels)) %>%
  ggplot(aes(nrmse, fset_label, colour = model_label, shape = model_label)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.5)) +
  facet_grid(paste0("h = ", horizon, " d") ~ scale_lab, scales = "free_x") +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_shape_manual(values = model_shapes, name = NULL) +
  labs(x = "Normalised RMSE on batch 4 (% of mean intake)", y = NULL) +
  theme_paper
save_fig(fig6, "Fig6_external_test_batch4.png", 10, 2.4 + 2 * length(horizons))

# Fig 7 -> manuscript Figure 3 : grouped permutation importance
fig7 <- perm_imp %>%
  prep_scale() %>%
  mutate(test_lab = batch_labels[test_batch],
         group = factor(group, levels = rev(names(group_colours)))) %>%
  ggplot(aes(pct_increase, group, shape = test_batch)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = group), size = 2.4,
             position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  facet_wrap(~ scale_lab, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_shape_manual(values = batch_shapes, labels = batch_labels, name = "Test batch") +
  scale_x_continuous(trans = "pseudo_log", breaks = c(0, 5, 20, 100, 500)) +
  labs(x = "Increase in RMSE when the group is permuted (%, pseudo-log scale)", y = NULL) +
  theme_paper
save_fig(fig7, "Fig7_grouped_permutation_importance.png", 10, 3.8)

# Fig 8 -> manuscript Figure 4 : ALE curves (pen, RF relative target, full set)
if (is.null(ale_rug)) {
  message("  Fig8 (ALE) skipped: T13c_ale_rug_sample.csv not found.")
} else {
  ale_order <- ale_range$feature_label
  fig8 <- ale_tbl %>%
    mutate(feature_label = factor(feature_label, levels = ale_order)) %>%
    ggplot(aes(x, ale, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.8) +
    geom_rug(data = ale_rug %>%
               mutate(feature_label = factor(unname(feature_labels[feature]),
                                             levels = ale_order)),
             aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.25, linewidth = 0.2) +
    facet_wrap(~ feature_label, scales = "free_x", ncol = 4) +
    scale_colour_manual(values = group_colours, name = NULL) +
    labs(x = NULL, y = "Accumulated local effect (kg/animal/day)") +
    theme_paper
  save_fig(fig8, "Fig8_ale_pen_full_model.png", 10, 8.5)
}

# Fig 9 -> manuscript Figure 6 : prediction intervals, batch 4, pen level
if (!is.null(all_preds)) {
  fig9_data <- all_preds %>%
    filter(part == "B", scale == "pen", horizon == h1,
           model == "rf_rel", fset == final_fset) %>%
    mutate(status = case_when(y < q05 ~ "Below 90% interval",
                              y > q95 ~ "Above 90% interval",
                              TRUE ~ "Inside interval"),
           pen_lab = factor(paste("Pen", pen), levels = paste("Pen", sort(unique(pen)))))
  if (nrow(fig9_data) > 0) {
    fig9 <- ggplot(fig9_data, aes(day)) +
      geom_ribbon(aes(ymin = q05, ymax = q95),
                  fill = model_colours[["RF (relative target)"]], alpha = 0.15) +
      geom_line(aes(y = pred), colour = model_colours[["RF (relative target)"]],
                linewidth = 0.5) +
      geom_point(aes(y = y, shape = status, colour = status), size = 1.1) +
      facet_wrap(~ pen_lab, ncol = 2) +
      scale_colour_manual(values = status_colours, name = NULL) +
      scale_shape_manual(values = status_shapes, name = NULL) +
      labs(x = "Fattening day", y = "Feed intake (kg/animal/day)") +
      theme_paper
    save_fig(fig9, "Fig9_prediction_intervals_batch4_pen.png", 10, 11)
  }
}

# Fig 10 -> supplementary Figure S3 : learning curve (pen, first horizon)
# kg and relative target in separate rows, each with its own y scale
fig10 <- learning_curve %>%
  filter(model %in% c("rf_level", "rf_rel")) %>%
  prep_model() %>%
  mutate(test_lab = batch_labels[test_batch]) %>%
  ggplot(aes(n_train_batches, 100 * skill_vs_persistence)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = model_label, shape = model_label), size = 2.2) +
  stat_summary(aes(colour = model_label), fun = mean, geom = "line",
               linewidth = 0.6) +
  facet_grid(model_label ~ test_lab, scales = "free_y") +
  scale_x_continuous(breaks = 1:3) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_shape_manual(values = model_shapes, guide = "none") +
  labs(x = "Number of training batches",
       y = "RMSE reduction vs persistence (%)") +
  theme_paper +
  theme(strip.text.y = element_text(angle = 0))
save_fig(fig10, "Fig10_learning_curve_pen.png", 10, 5.4)

# Fig S1 -> supplementary Figure S4 : linear regression vs RF
figS1 <- lm_vs_rf %>%
  prep_scale() %>%
  mutate(fset_label = factor(unname(fset_labels[fset]),
                             levels = rev(unname(fset_labels))),
         horizon_lab = paste0("h = ", horizon, " d")) %>%
  ggplot(aes(pct_diff_lm_vs_rf, fset_label, shape = test_batch)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(size = 2.2, colour = model_colours[["Linear regression"]],
             position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  facet_grid(horizon_lab ~ scale_lab, scales = "free_x") +
  scale_shape_manual(values = batch_shapes, labels = batch_labels, name = "Test batch") +
  labs(x = "RMSE of linear regression relative to RF with relative target (%)  (positive = RF better)",
       y = NULL) +
  theme_paper
save_fig(figS1, "FigS1_linear_vs_rf.png", 10, 2.8 + 2 * length(horizons))
