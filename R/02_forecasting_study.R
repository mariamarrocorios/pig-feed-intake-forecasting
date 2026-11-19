# Forecasting daily feed intake of growing-finishing pigs with Random Forest:
# aggregation scale (animal, pen, batch), predictor groups and forecast horizon.
#
# Code for: Effects of aggregation scale, predictor groups and forecast horizon
# on feed intake forecasting in growing-finishing pigs (Marroco et al.).
#
# Validation
#   Part A: leave-one-batch-out on batches 1-3, all predictor groups.
#   Part B: train on batches 1-3, test on batch 4 (no environmental data).
#   Within each scale, horizon and test batch, every model and predictor set is
#   evaluated on the same rows. Hyperparameters are tuned by an inner
#   leave-one-batch-out on the training batches.
#
# Input:  outputs/tables/ (made by 01_daily_tables.R)
# Output: outputs/results/ (tables/, figures/, predictions/, study_summary.txt)
#
# Quick test (horizon 1 only, smaller grids, fewer trees): set quick <- TRUE
# in run_all.R. The results of a quick run are only to check that everything works.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ranger); library(ggplot2)
})

# ============================================================
# 0. Settings
# ============================================================

if (!exists("tables_path")) tables_path <- "outputs/tables"
if (!exists("out_path"))    out_path    <- "outputs/results"
if (!exists("quick"))       quick       <- FALSE

out_tables  <- file.path(out_path, "tables")
out_figures <- file.path(out_path, "figures")
out_preds   <- file.path(out_path, "predictions")
for (p in c(out_tables, out_figures, out_preds)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

final_fset <- "F+B"     # final model of the paper: feed history + feeding behaviour
fig_dpi    <- 600
fig_pdf    <- FALSE     # TRUE also writes a pdf copy of each figure

horizons        <- if (quick) 1 else c(1, 3, 7)
h1              <- min(horizons)
scales          <- c("animal", "pen", "batch")
batches_A       <- c("batch1", "batch2", "batch3")   # leave-one-batch-out
batches_B_train <- c("batch1", "batch2", "batch3")
batch_B_test    <- "batch4"                          # external test batch

trim_last_days  <- 2      # last days of each batch removed (pigs leaving for slaughter)
ratio_floor_kg  <- 0.10   # lower bound of the denominator of the relative target
n_trees_final   <- if (quick) 200 else 500
n_trees_tune    <- if (quick) 100 else 250
n_threads       <- max(1, parallel::detectCores())
rf_max_rows     <- 10000  # max rows per tree (only matters at animal level)
stability_seeds <- if (quick) 1:2 else 1:5
perm_reps       <- if (quick) 2 else 10
ale_bins        <- 20
interval_probs  <- c(0.025, 0.05, 0.5, 0.95, 0.975)
hot_day_temp_c  <- 26

set.seed(123)

# Predictor groups
feat_F  <- c("f_lag_a", "f_lag_b", "f_lag_c", "f_mean5")      # feed intake history
feat_G  <- c("g_day", "g_weight", "g_nanimals")              # growth
feat_B  <- c("b_visits", "b_visits3", "b_duration", "b_rate") # feeding behaviour
feat_E  <- c("e_temp", "e_hum", "e_co2", "e_nh3")            # barn environment
feat_EO <- c("eo_temp", "eo_hum")                            # climate of the target day

feature_sets_A <- list(
  "F"      = feat_F,
  "F+G"    = c(feat_F, feat_G),
  "F+B"    = c(feat_F, feat_B),
  "F+E"    = c(feat_F, feat_E),
  "F+E+EO" = c(feat_F, feat_E, feat_EO),
  "Full"   = c(feat_F, feat_G, feat_B, feat_E)
)
feature_sets_B <- list(
  "F"     = feat_F,
  "F+G"   = c(feat_F, feat_G),
  "F+B"   = c(feat_F, feat_B),
  "F+G+B" = c(feat_F, feat_G, feat_B)
)
importance_set <- c(A = "Full", B = "F+G+B")    # set used for importance and ALE
predictor_groups <- list(
  "Feed history" = feat_F,
  "Growth"       = feat_G,
  "Behaviour"    = feat_B,
  "Environment"  = c(feat_E, feat_EO)
)

baseline_models <- c("persistence", "mean5")
rf_models       <- c("rf_level", "rf_rel")       # RF on kg and on the relative target
supp_models     <- "lm"                          # linear regression (supplementary)
fitted_models   <- c(rf_models, supp_models)
main_models     <- c(baseline_models, rf_models)

source("R/functions.R")

# Without climate data, the sets with the environment group are dropped
if (nrow(env_pen) == 0) {
  message("No environmental data: the environment group is left out.")
  feature_sets_A <- feature_sets_B
  importance_set["A"] <- "F+G+B"
}

# ============================================================
# 1-5. Main loop
# ============================================================

folds_A <- lapply(batches_A, function(b) list(train = setdiff(batches_A, b), test = b))
folds_B <- list(list(train = batches_B_train, test = batch_B_test))

t0 <- Sys.time()
results <- list()
for (sc in scales) {
  results[[paste0("A_", sc)]] <- run_part("A", sc, feature_sets_A, folds_A)
  results[[paste0("B_", sc)]] <- run_part("B", sc, feature_sets_B, folds_B)
}
message("Main loop finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

all_preds <- bind_rows(lapply(results, `[[`, "preds")) %>% mutate(error = y - pred)
tuned_params <- bind_rows(lapply(results, `[[`, "tuned"))
perm_imp <- bind_rows(lapply(results, `[[`, "imp"))

saveRDS(all_preds, file.path(out_preds, "predictions_all.rds"))
write_csv(all_preds %>% filter(scale != "animal"),
          file.path(out_preds, "predictions_pen_batch.csv"))
write_csv(tuned_params, file.path(out_tables, "T0_tuned_hyperparameters.csv"))

# ============================================================
# 6. Accuracy metrics
# ============================================================

metrics <- all_preds %>%
  group_by(part, scale, horizon, test_batch, fset, model) %>%
  summarise(n = n(), mean_obs = mean(y),
            rmse = rmse_fun(y, pred), mae = mae_fun(y, pred),
            mape = mape_fun(y, pred), bias = mean(pred - y), .groups = "drop") %>%
  mutate(nrmse = 100 * rmse / mean_obs) %>%
  group_by(part, scale, horizon, test_batch) %>%
  mutate(rmse_persistence = rmse[model == "persistence"][1],
         skill_vs_persistence = 1 - rmse / rmse_persistence) %>%
  ungroup() %>%
  mutate(model_label = model_labels[model],
         fset_label = ifelse(fset == "none", fset_label_none, fset_labels[fset]))

summarise_metrics <- function(m) {
  m %>%
    group_by(part, scale, horizon, fset, fset_label, model, model_label) %>%
    summarise(n_folds = n(),
              rmse_mean = mean(rmse), rmse_min = min(rmse), rmse_max = max(rmse),
              mae_mean = mean(mae),
              nrmse_mean = mean(nrmse), nrmse_min = min(nrmse), nrmse_max = max(nrmse),
              mape_mean = mean(mape), bias_mean = mean(bias),
              skill_mean = mean(skill_vs_persistence),
              skill_min = min(skill_vs_persistence),
              skill_max = max(skill_vs_persistence), .groups = "drop") %>%
    arrange(part, scale, horizon, rmse_mean)
}

metrics_main <- metrics %>% filter(model %in% main_models)
summary_main <- summarise_metrics(metrics_main)
write_csv(metrics_main, file.path(out_tables, "T1_metrics_by_fold.csv"))
write_csv(summary_main, file.path(out_tables, "T2_metrics_summary.csv"))

best_rf <- summary_main %>%
  filter(model %in% rf_models) %>%
  group_by(part, scale, horizon) %>%
  slice_min(rmse_mean, n = 1, with_ties = FALSE) %>%
  ungroup()
write_csv(best_rf, file.path(out_tables, "T3_best_rf_by_scale_horizon.csv"))

added_value <- metrics %>%
  filter(model %in% fitted_models) %>%
  group_by(part, scale, horizon, test_batch, model) %>%
  mutate(rmse_F = rmse[fset == "F"][1],
         pct_change_vs_F = 100 * (rmse - rmse_F) / rmse_F) %>%
  ungroup() %>%
  filter(fset != "F")
write_csv(added_value %>% filter(model %in% rf_models),
          file.path(out_tables, "T4_added_value_vs_feed_history.csv"))

# ============================================================
# 7. RQ1 - direct vs bottom-up forecasts (animal -> pen, pen -> batch)
# ============================================================

pen_obs <- load_scale("pen") %>% select(batch, pen, day, n_animals)
join_keys <- c("part", "horizon", "test_batch", "fset", "model")

bottom_up_pen <- all_preds %>%
  filter(scale == "animal", single_pen) %>%
  group_by(across(all_of(join_keys)), batch, pen, day) %>%
  summarise(pred_bu = mean(pred), n_pred = n(), .groups = "drop") %>%
  inner_join(pen_obs, by = c("batch", "pen", "day")) %>%
  filter(n_pred == n_animals) %>%                    # every pig of the pen forecast
  inner_join(all_preds %>% filter(scale == "pen") %>%
               select(all_of(join_keys), batch, pen, day, y, pred_direct = pred),
             by = c(join_keys, "batch", "pen", "day")) %>%
  mutate(target = "pen", source = "animal")

bottom_up_batch <- all_preds %>%
  filter(scale == "pen") %>%
  inner_join(pen_obs, by = c("batch", "pen", "day")) %>%
  group_by(across(all_of(join_keys)), batch, day) %>%
  summarise(pred_bu = sum(pred * n_animals) / sum(n_animals), n_pens = n(),
            .groups = "drop") %>%
  group_by(batch) %>% filter(n_pens == max(n_pens)) %>% ungroup() %>%   # every pen forecast that day
  inner_join(all_preds %>% filter(scale == "batch") %>%
               select(all_of(join_keys), batch, day, y, pred_direct = pred),
             by = c(join_keys, "batch", "day")) %>%
  mutate(target = "batch", source = "pen")

bottom_up <- bind_rows(bottom_up_pen, bottom_up_batch) %>%
  group_by(part, target, source, horizon, test_batch, fset, model) %>%
  summarise(n = n(), rmse_direct = rmse_fun(y, pred_direct),
            rmse_bottom_up = rmse_fun(y, pred_bu), .groups = "drop") %>%
  mutate(pct_diff_bottom_up = 100 * (rmse_bottom_up - rmse_direct) / rmse_direct,
         model_label = model_labels[model])
write_csv(bottom_up %>% filter(model %in% main_models),
          file.path(out_tables, "T5_direct_vs_bottom_up.csv"))

# ============================================================
# 8. Diebold-Mariano tests
# ============================================================

run_dm <- function(pairs) {
  res <- list()
  for (i in seq_len(nrow(pairs))) {
    pr <- pairs[i, ]
    a <- all_preds %>% filter(model == pr$m1, fset == pr$f1) %>%
      select(part, scale, horizon, test_batch, entity, day, e1 = error)
    b <- all_preds %>% filter(model == pr$m2, fset == pr$f2) %>%
      select(part, scale, horizon, test_batch, entity, day, e2 = error)
    ab <- inner_join(a, b, by = c("part", "scale", "horizon", "test_batch", "entity", "day"))
    if (nrow(ab) == 0) next
    res[[length(res) + 1]] <- ab %>%
      group_by(part, scale, horizon, test_batch) %>%
      group_modify(~ {
        t <- dm_test(.x$e1, .x$e2, .x$day, .y$horizon)
        tibble(n = nrow(.x), rmse_1 = sqrt(mean(.x$e1^2)), rmse_2 = sqrt(mean(.x$e2^2)),
               dm_stat = t[["dm_stat"]], p_value = t[["p_value"]], n_days = t[["n_days"]])
      }) %>%
      ungroup() %>%
      mutate(question = pr$question,
             model_1 = paste(model_labels[pr$m1], if (pr$f1 != "none") pr$f1 else ""),
             model_2 = paste(model_labels[pr$m2], if (pr$f2 != "none") pr$f2 else ""))
  }
  bind_rows(res) %>%
    mutate(model_1 = trimws(model_1), model_2 = trimws(model_2),
           result = case_when(
             is.na(p_value) ~ NA_character_,
             p_value >= 0.05 ~ "no significant difference",
             rmse_1 < rmse_2 ~ paste(model_1, "better"),
             TRUE ~ paste(model_2, "better")))
}

rq2_sets <- c("F+G", "F+B", "F+E", "F+E+EO", "Full", "F+G+B")

dm_main <- run_dm(bind_rows(
  tibble(m1 = "rf_rel", f1 = "F", m2 = "persistence", f2 = "none",
         question = "RF (relative) vs persistence"),
  tibble(m1 = "rf_rel", f1 = "F", m2 = "mean5", f2 = "none",
         question = "RF (relative) vs 5-day mean"),
  tibble(m1 = "rf_rel", f1 = "F", m2 = "rf_level", f2 = "F",
         question = "RF relative vs kg target"),
  tibble(m1 = "rf_rel", f1 = rq2_sets, m2 = "rf_rel", f2 = "F",
         question = paste("RQ2", rq2_sets, "vs F")),
  # The final model of Section 3.4 against the two naive references.
  # The three comparisons above use feed intake history alone; these use
  # the predictor set reported as the final model, and are the ones the
  # manuscript quotes for the external test batch.
  tibble(m1 = "rf_rel", f1 = final_fset, m2 = "persistence", f2 = "none",
         question = "Final model vs persistence"),
  tibble(m1 = "rf_rel", f1 = final_fset, m2 = "mean5", f2 = "none",
         question = "Final model vs 5-day mean")
))
write_csv(dm_main, file.path(out_tables, "T6_diebold_mariano.csv"))

# ============================================================
# 9. Which factor matters most for the RF? (share of variance in nRMSE)
# ============================================================

vd_data <- metrics %>%
  filter(part == "A", model == "rf_rel") %>%
  mutate(across(c(scale, fset, test_batch), factor), horizon = factor(horizon))

vd_formula <- if (length(horizons) > 1) {
  nrmse ~ scale + horizon + fset + test_batch +
    scale:test_batch + scale:horizon + scale:fset + horizon:test_batch
} else {
  nrmse ~ scale + fset + test_batch + scale:test_batch + scale:fset
}
vd_tab <- function(fit) {
  a <- anova(fit)
  tibble(term = trimws(rownames(a)), df = a$Df, sum_sq = a$`Sum Sq`) %>%
    mutate(pct_of_total = 100 * sum_sq / sum(sum_sq)) %>%
    arrange(desc(pct_of_total))
}
variance_decomposition <- vd_tab(lm(vd_formula, data = vd_data))
write_csv(variance_decomposition, file.path(out_tables, "T7_variance_decomposition_nrmse.csv"))

vd_within <- bind_rows(lapply(scales, function(sc) {
  f <- if (length(horizons) > 1) {
    nrmse ~ horizon + fset + test_batch + horizon:test_batch + fset:test_batch
  } else {
    nrmse ~ fset + test_batch
  }
  vd_tab(lm(f, data = vd_data %>% filter(scale == sc) %>% droplevels())) %>%
    mutate(scale = sc, .before = 1)
}))
write_csv(vd_within, file.path(out_tables, "T7b_variance_decomposition_within_scale.csv"))

# ============================================================
# 10. Diagnostics: phase, pen, hot days, seed stability
# ============================================================

by_phase <- all_preds %>%
  filter(scale == "pen", horizon == h1, fset %in% c("none", final_fset),
         model %in% main_models) %>%
  mutate(phase = cut(day, c(0, 30, 60, 90, Inf),
                     labels = c("Days 1-30", "Days 31-60", "Days 61-90", "Days > 90"))) %>%
  group_by(part, test_batch, phase, model) %>%
  summarise(n = n(), rmse = rmse_fun(y, pred), bias = mean(pred - y),
            nrmse = 100 * rmse / mean(y), .groups = "drop")
write_csv(by_phase, file.path(out_tables, "T8_error_by_phase_pen.csv"))

by_pen <- all_preds %>%
  filter(scale == "pen", horizon == h1, fset %in% c("none", final_fset),
         model %in% main_models) %>%
  group_by(part, test_batch, pen, model) %>%
  summarise(n = n(), rmse = rmse_fun(y, pred), bias = mean(pred - y), .groups = "drop")
write_csv(by_pen, file.path(out_tables, "T9_error_by_pen.csv"))

hot_days <- all_preds %>%
  filter(part == "A", scale == "pen", model == "rf_rel", fset %in% c("F", "F+E", "F+E+EO")) %>%
  mutate(day_type = ifelse(pmax(e_temp, eo_temp) > hot_day_temp_c,
                           paste0("Hot (> ", hot_day_temp_c, " C)"), "Other")) %>%
  group_by(horizon, test_batch, day_type, fset) %>%
  summarise(n = n(), rmse = rmse_fun(y, pred), .groups = "drop") %>%
  group_by(horizon, test_batch, day_type) %>%
  mutate(pct_change_vs_F = 100 * (rmse - rmse[fset == "F"][1]) / rmse[fset == "F"][1]) %>%
  ungroup()
write_csv(hot_days, file.path(out_tables, "T10_hot_days_environment_check.csv"))

pen_base <- load_scale("pen") %>% trim_batch_end()
pen_h1 <- build_features(pen_base, h1) %>%
  filter(!is.na(y), if_all(all_of(feat_F), ~ !is.na(.x)))

get_params <- function(part_, scale_, fset_, model_, batch_ = NULL) {
  t <- tuned_params %>% filter(part == part_, scale == scale_, fset == fset_, model == model_)
  if (!is.null(batch_)) t <- t %>% filter(test_batch == batch_)
  # most frequent combination across folds
  t <- t %>% count(mtry, min.node.size, num.trees, sort = TRUE) %>% slice(1)
  list(mtry = t$mtry, min.node.size = t$min.node.size, num.trees = t$num.trees)
}

seed_stab <- list()
for (fo in folds_A) {
  tr <- pen_h1 %>% filter(batch %in% fo$train)
  te <- pen_h1 %>% filter(batch == fo$test)
  for (m in rf_models) {
    pr <- get_params("A", "pen", "F", m, fo$test)
    for (s in stability_seeds) {
      seed_stab[[length(seed_stab) + 1]] <- tibble(
        test_batch = fo$test, model = m, seed = s,
        rmse = rmse_fun(te$y, fit_model(m, tr, feat_F, pr, seed = s)$predict(te)))
    }
  }
}
seed_stab <- bind_rows(seed_stab) %>%
  group_by(test_batch, model) %>%
  summarise(n_seeds = n(), rmse_mean = mean(rmse), rmse_sd = sd(rmse), .groups = "drop")
write_csv(seed_stab, file.path(out_tables, "T11_seed_stability_pen_F.csv"))

# Same seed check on the final model (feed history + behaviour), for the
# three leave-one-batch-out folds and for batch 4
seed_final <- list()
for (pt in c("A", "B")) {
  fsets <- if (pt == "A") feature_sets_A else feature_sets_B
  folds <- if (pt == "A") folds_A else folds_B
  needed <- unique(unlist(lapply(fsets, scale_features, scale = "pen")))
  dh <- build_features(pen_base, h1)
  vars <- feature_sets_B[[final_fset]]
  for (fo in folds) {
    d <- dh %>% filter(batch %in% c(fo$train, fo$test)) %>%
      filter(!is.na(y), if_all(all_of(needed), ~ !is.na(.x)))
    tr <- d %>% filter(batch %in% fo$train)
    te <- d %>% filter(batch == fo$test)
    pr <- get_params(pt, "pen", final_fset, "rf_rel", fo$test)
    for (s in stability_seeds) {
      seed_final[[length(seed_final) + 1]] <- tibble(
        part = pt, test_batch = fo$test, seed = s,
        rmse = rmse_fun(te$y, fit_model("rf_rel", tr, vars, pr, seed = s)$predict(te)))
    }
  }
}
seed_final <- bind_rows(seed_final) %>%
  group_by(part, test_batch) %>%
  summarise(n_seeds = n(), rmse_mean = mean(rmse), rmse_sd = sd(rmse),
            rmse_range = max(rmse) - min(rmse), .groups = "drop")
write_csv(seed_final, file.path(out_tables, "T11b_seed_stability_pen_final.csv"))


# ============================================================
# 11. Interpretation: grouped permutation importance and ALE
# ============================================================

write_csv(perm_imp, file.path(out_tables, "T12_grouped_permutation_importance.csv"))

# ALE for the pen-level RF (relative target, first horizon, full predictor
# set), refitted on all of batches 1-3.
pen_full_vars <- feature_sets_A[[importance_set[["A"]]]]
pen_A <- build_features(pen_base, h1) %>%
  filter(batch %in% batches_A, !is.na(y), if_all(all_of(pen_full_vars), ~ !is.na(.x)))
ale_fit <- fit_model("rf_rel", pen_A, pen_full_vars, get_params("A", "pen", importance_set[["A"]], "rf_rel"))
ale_tbl <- bind_rows(lapply(pen_full_vars, function(v) ale_1d(ale_fit$predict, pen_A, v))) %>%
  mutate(feature_label = feature_labels[feature],
         group = case_when(feature %in% feat_F ~ "Feed history",
                           feature %in% feat_G ~ "Growth",
                           feature %in% feat_B ~ "Behaviour",
                           TRUE ~ "Environment"))
ale_range <- ale_tbl %>%
  group_by(group, feature, feature_label) %>%
  summarise(ale_range_kg = max(ale) - min(ale), .groups = "drop") %>%
  arrange(desc(ale_range_kg))
write_csv(ale_tbl, file.path(out_tables, "T13_ale_pen_full_model.csv"))
write_csv(ale_range, file.path(out_tables, "T13b_ale_range_pen_full_model.csv"))

# sample of predictor values for the rug under the ALE curves (03_figures.R)
set.seed(1)
ale_rug <- pen_A %>%
  select(all_of(pen_full_vars)) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "x") %>%
  group_by(feature) %>% slice_sample(n = 200) %>% ungroup()
write_csv(ale_rug, file.path(out_tables, "T13c_ale_rug_sample.csv"))

# ============================================================
# 12. Prediction intervals (quantile regression forest)
# ============================================================

intervals <- all_preds %>%
  filter(model == "rf_rel", !is.na(q05)) %>%
  group_by(part, scale, horizon, test_batch, fset) %>%
  summarise(
    n = n(),
    coverage_90 = mean(y >= q05 & y <= q95),
    coverage_95 = mean(y >= q025 & y <= q975),
    width_90_kg = mean(q95 - q05),
    width_95_kg = mean(q975 - q025),
    width_90_pct = 100 * mean((q95 - q05) / pred),
    below_90 = mean(y < q05),                  # potential "low intake" alarms
    above_90 = mean(y > q95),
    interval_score_90 = mean((q95 - q05) + (2 / 0.10) * (q05 - y) * (y < q05) +
                               (2 / 0.10) * (y - q95) * (y > q95)),
    .groups = "drop")
write_csv(intervals, file.path(out_tables, "T14_prediction_intervals.csv"))

# ============================================================
# 13. Learning curve: how many training batches are needed? (pen, first horizon)
# ============================================================

lc_params <- get_params("A", "pen", "F", "rf_rel")
lc_params_kg <- get_params("A", "pen", "F", "rf_level")
all_batches <- c(batches_A, batch_B_test)
lc <- list()
for (tb in all_batches) {
  pool <- if (tb == batch_B_test) batches_B_train else setdiff(batches_A, tb)
  te <- pen_h1 %>% filter(batch == tb)
  for (k in seq_along(pool)) {
    for (cmb in combn(pool, k, simplify = FALSE)) {
      tr <- pen_h1 %>% filter(batch %in% cmb)
      p_list <- list(
        persistence = te$f_lag_a,
        rf_rel   = fit_model("rf_rel", tr, feat_F, lc_params)$predict(te),
        rf_level = fit_model("rf_level", tr, feat_F, lc_params_kg)$predict(te),
        lm       = fit_model("lm", tr, feat_F, list())$predict(te)
      )
      for (m in names(p_list)) {
        lc[[length(lc) + 1]] <- tibble(test_batch = tb, n_train_batches = k,
                                       train_batches = paste(cmb, collapse = "+"),
                                       model = m, n = nrow(te),
                                       rmse = rmse_fun(te$y, p_list[[m]]))
      }
    }
  }
}
learning_curve <- bind_rows(lc) %>%
  group_by(test_batch, train_batches) %>%
  mutate(skill_vs_persistence = 1 - rmse / rmse[model == "persistence"][1]) %>%
  ungroup()
write_csv(learning_curve %>% filter(model != "lm"),
          file.path(out_tables, "T15_learning_curve_pen.csv"))

# ============================================================
# 14. Supplementary material: linear regression
# ============================================================

supp_metrics <- metrics %>% filter(model %in% c("lm", "rf_rel"))
write_csv(summarise_metrics(supp_metrics), file.path(out_tables, "S1_linear_regression_metrics.csv"))

lm_vs_rf <- supp_metrics %>%
  select(part, scale, horizon, test_batch, fset, model, rmse) %>%
  pivot_wider(names_from = model, values_from = rmse) %>%
  mutate(pct_diff_lm_vs_rf = 100 * (lm - rf_rel) / rf_rel) %>%
  rename(rmse_lm = lm, rmse_rf_rel = rf_rel)
write_csv(lm_vs_rf, file.path(out_tables, "S2_linear_vs_rf_by_fold.csv"))

dm_supp <- run_dm(bind_rows(
  tibble(m1 = "lm", f1 = unique(c(names(feature_sets_A), names(feature_sets_B))),
         m2 = "rf_rel", f2 = unique(c(names(feature_sets_A), names(feature_sets_B))),
         question = "Linear regression vs RF (relative), same predictors"),
  tibble(m1 = "lm", f1 = rq2_sets, m2 = "lm", f2 = "F",
         question = paste("RQ2 (linear)", rq2_sets, "vs F"))
))
write_csv(dm_supp, file.path(out_tables, "S3_diebold_mariano_linear.csv"))
write_csv(learning_curve, file.path(out_tables, "S4_learning_curve_incl_linear.csv"))
write_csv(added_value %>% filter(model == "lm"),
          file.path(out_tables, "S5_added_value_linear.csv"))

# ============================================================
# 15. Figures (all defined in 03_figures.R)
# ============================================================

source("R/03_figures.R")

# ============================================================
# 16. Text summary
# ============================================================

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
lines <- c(
  "Forecasting feed intake: aggregation scale, predictor groups and horizon",
  "=======================================================================",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M"), if (quick) "(QUICK MODE - only to check that the code runs)" else ""),
  paste("Horizons:", paste(horizons, collapse = ", ")),
  "Part A: leave-one-batch-out on batches 1-3 (all predictor groups)",
  "Part B: train batches 1-3, external test on batch 4 (no environment)",
  "",
  "Best RF per scale and horizon (mean RMSE over test batches):"
)
for (i in seq_len(nrow(best_rf))) {
  r <- best_rf[i, ]
  lines <- c(lines, sprintf("  %s | %-6s | h=%d | %-22s | %-32s | RMSE %s | nRMSE %s%% | skill %s%%",
                            r$part, r$scale, r$horizon, r$model_label, r$fset_label,
                            fmt(r$rmse_mean), fmt(r$nrmse_mean, 1), fmt(100 * r$skill_mean, 1)))
}
lines <- c(lines, "", "Share of variance in nRMSE of the RF (relative target), Part A:")
for (i in seq_len(min(8, nrow(variance_decomposition)))) {
  r <- variance_decomposition[i, ]
  lines <- c(lines, sprintf("  %-22s %5.1f%%", r$term, r$pct_of_total))
}
lines <- c(lines, "", "Within each scale:")
for (sc in scales) {
  v <- vd_within %>% filter(scale == sc)
  lines <- c(lines, paste0("  ", sc, ": ",
                           paste(sprintf("%s %.1f%%", v$term, v$pct_of_total), collapse = " | ")))
}
lines <- c(lines, "", "Mean % change in RMSE vs feed history, RF relative target (Part A):")
av <- added_value %>% filter(part == "A", model == "rf_rel") %>%
  group_by(scale, fset) %>%
  summarise(mean_pct = mean(pct_change_vs_F), share_better = mean(pct_change_vs_F < 0),
            .groups = "drop")
for (i in seq_len(nrow(av))) {
  r <- av[i, ]
  lines <- c(lines, sprintf("  %-6s %-7s %+6.1f%%  (better in %3.0f%% of fold x horizon)",
                            r$scale, r$fset, r$mean_pct, 100 * r$share_better))
}
lines <- c(lines, "", "90% prediction interval coverage (RF relative, feed history):")
iv <- intervals %>% filter(fset == "F") %>%
  group_by(part, scale, horizon) %>%
  summarise(cov = mean(coverage_90), width = mean(width_90_pct), .groups = "drop")
for (i in seq_len(nrow(iv))) {
  r <- iv[i, ]
  lines <- c(lines, sprintf("  %s | %-5s | h=%d | coverage %.1f%% | width %.1f%% of forecast",
                            r$part, r$scale, r$horizon, 100 * r$cov, r$width))
}
lines <- c(lines, "", "Supplementary - linear regression vs RF (relative), mean % RMSE difference:")
lv <- lm_vs_rf %>% group_by(part, scale) %>%
  summarise(mean_pct = mean(pct_diff_lm_vs_rf), .groups = "drop")
for (i in seq_len(nrow(lv))) {
  lines <- c(lines, sprintf("  %s | %-6s | %+5.1f%% (positive = RF better)",
                            lv$part[i], lv$scale[i], lv$mean_pct[i]))
}
writeLines(lines, file.path(out_path, "study_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_path, "sessionInfo.txt"))
message(paste(lines, collapse = "\n"))
message("Total time: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
