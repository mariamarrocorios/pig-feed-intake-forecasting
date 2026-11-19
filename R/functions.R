# Functions used by 02_forecasting_study.R and 03_figures.R:
# metrics, data loading, predictors, model fitting, the main loop,
# and the labels, colours and theme of the figures.
# Settings (paths, horizons, number of trees...) are set in 02_forecasting_study.R.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ranger); library(ggplot2)
})

if (!exists("tables_path")) tables_path <- "outputs/tables"
if (!exists("out_path"))    out_path    <- "outputs/results"
if (!exists("out_tables"))  out_tables  <- file.path(out_path, "tables")
if (!exists("out_figures")) out_figures <- file.path(out_path, "figures")
if (!exists("out_preds"))   out_preds   <- file.path(out_path, "predictions")
if (!exists("quick"))       quick       <- FALSE


# ============================================================
# Metrics
# ============================================================

rmse_fun <- function(o, p) sqrt(mean((o - p)^2))
mae_fun  <- function(o, p) mean(abs(o - p))
mape_fun <- function(o, p) {
  k <- o >= 0.1                       # avoid dividing by ~0 at animal level
  100 * mean(abs((o[k] - p[k]) / o[k]))
}

# Diebold-Mariano test with the Harvey-Leybourne-Newbold correction.
# Losses of all pens (or animals) of the same day are averaged first, so the
# test runs on one daily series. HAC variance (Bartlett), lag (h - 1) + floor(T^(1/3)).
dm_test <- function(e1, e2, time, h) {
  d <- tapply(e1^2 - e2^2, time, mean)
  d <- d[order(as.numeric(names(d)))]
  n <- length(d)
  if (n < 10) return(c(dm_stat = NA, p_value = NA, n_days = n))
  dbar <- mean(d)
  L <- min(n - 1, (h - 1) + floor(n^(1 / 3)))
  gam <- sapply(0:L, function(k) {
    sum((d[(k + 1):n] - dbar) * (d[1:(n - k)] - dbar)) / n
  })
  w <- 1 - (0:L) / (L + 1)
  lrv <- gam[1] + 2 * sum(w[-1] * gam[-1])
  if (!is.finite(lrv) || lrv <= 0) return(c(dm_stat = NA, p_value = NA, n_days = n))
  stat <- dbar / sqrt(lrv / n) * sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  c(dm_stat = stat, p_value = 2 * pt(-abs(stat), df = n - 1), n_days = n)
}


# ============================================================
# Data at the three scales (one row per entity and day)
# ============================================================

trim_batch_end <- function(d) {
  d %>% group_by(batch) %>% filter(day <= max(day) - trim_last_days) %>% ungroup()
}

env_pen <- read_csv(file.path(tables_path, "environment_pen_day.csv"),
                    col_types = cols(batch = "c", date = "D", .default = "d")) %>%
  transmute(batch, location = as.integer(location), date = as.Date(date),
            temp = inside_temp_c_day_mean, hum = inside_humidity_pct_day_mean,
            co2 = co2_ppm_day_mean, nh3 = nh3_ppm_day_mean)

env_batch <- env_pen %>%
  group_by(batch, date) %>%
  summarise(across(c(temp, hum, co2, nh3), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

load_scale <- function(scale) {
  if (scale == "pen") {
    read_csv(file.path(tables_path, "pen_day.csv"), show_col_types = FALSE) %>%
      transmute(batch, location = as.integer(location), date = as.Date(date),
                entity = paste0(batch, "_P", location), pen = location,
                day = batch_day_fattening, y = mean_feed_kg,
                weight = mean_weight_kg, n_animals = n_animals_pen_day,
                visits = mean_visits, duration = mean_duration_h, rate = mean_rate_g_min) %>%
      left_join(env_pen, by = c("batch", "location", "date")) %>%
      select(-location)

  } else if (scale == "batch") {
    read_csv(file.path(tables_path, "batch_day.csv"), show_col_types = FALSE) %>%
      transmute(batch, date = as.Date(date), entity = batch, pen = NA_integer_,
                day = batch_day_fattening, y = mean_feed_kg,
                weight = mean_weight_kg, n_animals = n_animals_day,
                visits = mean_visits, duration = mean_duration_h, rate = mean_rate_g_min) %>%
      left_join(env_batch, by = c("batch", "date"))

  } else if (scale == "animal") {
    read_csv(file.path(tables_path, "animal_day.csv"), show_col_types = FALSE) %>%
      filter(!is.na(animal_number)) %>%
      transmute(batch, date = as.Date(date),
                entity = paste0(batch, "_A", animal_number),
                # if a pig was recorded in two pens on the same day, take the first one
                pen = suppressWarnings(as.integer(sub("[^0-9].*$", "", as.character(pens_day)))),
                single_pen = n_pens_day == 1,
                day = batch_day_fattening, y = daily_feed_kg,
                weight = mean_weight_kg, n_animals = NA_real_,
                visits = daily_visits, duration = daily_duration_sec / 3600,
                rate = mean_rate_g_min) %>%
      left_join(env_pen, by = c("batch", "pen" = "location", "date"))
  }
}


# ============================================================
# Predictors for each horizon (lags by calendar day)
# ============================================================

roll_mean <- function(x, k, min_obs) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    w <- x[max(1, i - k + 1):i]
    if (sum(!is.na(w)) >= min_obs) out[i] <- mean(w, na.rm = TRUE)
  }
  out
}

build_features <- function(d, h) {
  d %>%
    group_by(entity) %>%
    complete(day = full_seq(day, 1)) %>%      # missing days become NA, lags stay correct
    arrange(day, .by_group = TRUE) %>%
    mutate(
      observed = !is.na(y),
      m5 = roll_mean(y, 5, 4),
      v3 = roll_mean(visits, 3, 2),
      f_lag_a = lag(y, h), f_lag_b = lag(y, h + 1), f_lag_c = lag(y, h + 2),
      f_mean5 = lag(m5, h),
      g_day = day, g_weight = lag(weight, h), g_nanimals = lag(n_animals, h),
      b_visits = lag(visits, h), b_visits3 = lag(v3, h),
      b_duration = lag(duration, h), b_rate = lag(rate, h),
      e_temp = lag(temp, h), e_hum = lag(hum, h), e_co2 = lag(co2, h), e_nh3 = lag(nh3, h),
      eo_temp = temp, eo_hum = hum             # climate of the target day
    ) %>%
    ungroup() %>%
    filter(observed) %>%
    mutate(horizon = h)
}

scale_features <- function(vars, scale) {
  if (scale == "animal") setdiff(vars, "g_nanimals") else vars
}


# ============================================================
# Models
# ============================================================

rf_fraction <- function(n) min(0.632, rf_max_rows / n)
den_fun <- function(d) pmax(d$f_mean5, ratio_floor_kg)

# relative target: intake / 5-day mean. The lags are divided by the same
# denominator; the 5-day mean stays in kg as a proxy of the stage of fattening.
rel_frame <- function(d, vars) {
  den <- den_fun(d)
  x <- as.data.frame(d[, vars, drop = FALSE])
  for (v in intersect(c("f_lag_a", "f_lag_b", "f_lag_c"), vars)) x[[v]] <- d[[v]] / den
  x
}

predict_baseline <- function(model, d) {
  if (model == "persistence") d$f_lag_a else d$f_mean5
}

fit_model <- function(model, train, vars, params, seed = 1, quantiles = FALSE) {
  if (model == "lm") {
    fit <- lm(as.formula(paste("y ~", paste(vars, collapse = " + "))), data = train)
    return(list(predict = function(nd) as.numeric(predict(fit, newdata = nd)),
                predict_q = NULL))
  }
  rel <- model == "rf_rel"
  xf <- function(nd) if (rel) rel_frame(nd, vars) else as.data.frame(nd[, vars, drop = FALSE])
  sc <- function(nd) if (rel) den_fun(nd) else 1
  y  <- if (rel) train$y / den_fun(train) else train$y
  rf <- ranger(x = xf(train), y = y,
               num.trees = params$num.trees, mtry = min(params$mtry, length(vars)),
               min.node.size = params$min.node.size,
               replace = FALSE, sample.fraction = rf_fraction(nrow(train)),
               quantreg = quantiles, num.threads = n_threads, seed = seed)
  list(
    predict = function(nd) {
      sc(nd) * predict(rf, data = xf(nd), num.threads = n_threads)$predictions
    },
    predict_q = if (quantiles) function(nd, probs) {
      sc(nd) * predict(rf, data = xf(nd), type = "quantiles", quantiles = probs,
                       num.threads = n_threads)$predictions
    } else NULL
  )
}

param_grid <- function(model, p, scale) {
  if (!model %in% rf_models) return(NULL)
  nodes <- switch(scale, animal = c(20, 50, 150), pen = c(5, 20, 50), batch = c(3, 10, 25))
  mtrys <- unique(pmax(1, round(p * c(0.25, 0.5, 0.8))))
  g <- expand.grid(mtry = mtrys, min.node.size = nodes)
  if (quick) g <- g[c(1, nrow(g)), , drop = FALSE]
  g$num.trees <- n_trees_tune
  g
}

# tuning by inner leave-one-batch-out, using the training batches only
tune_model <- function(model, train, vars, scale) {
  grid <- param_grid(model, length(vars), scale)
  if (is.null(grid)) return(list())
  inner <- unique(train$batch)
  score <- sapply(seq_len(nrow(grid)), function(i) {
    pr <- as.list(grid[i, , drop = FALSE])
    mean(sapply(inner, function(b) {
      tr <- train[train$batch != b, ]
      va <- train[train$batch == b, ]
      rmse_fun(va$y, fit_model(model, tr, vars, pr)$predict(va))
    }))
  })
  best <- as.list(grid[which.min(score), , drop = FALSE])
  best$num.trees <- n_trees_final
  best
}

# grouped permutation importance: all the variables of a group are shuffled together
perm_importance <- function(fit, test, vars, reps) {
  base <- rmse_fun(test$y, fit$predict(test))
  bind_rows(lapply(names(predictor_groups), function(g) {
    gv <- intersect(predictor_groups[[g]], vars)
    if (length(gv) == 0) return(NULL)
    inc <- sapply(seq_len(reps), function(r) {
      set.seed(r)
      tp <- test
      tp[gv] <- test[sample(nrow(test)), gv]
      rmse_fun(test$y, fit$predict(tp)) - base
    })
    tibble(group = g, n_vars = length(gv), rmse_full = base,
           rmse_increase = mean(inc), rmse_increase_sd = sd(inc),
           pct_increase = 100 * mean(inc) / base)
  }))
}

# first-order accumulated local effects (Apley and Zhu, 2020), in kg
ale_1d <- function(predict_fun, X, feature, K = ale_bins) {
  x <- X[[feature]]
  z <- unique(quantile(x, probs = seq(0, 1, length.out = K + 1), type = 1, na.rm = TRUE))
  if (length(z) < 3) return(NULL)
  nb <- length(z) - 1
  bin <- pmin(pmax(findInterval(x, z, rightmost.closed = TRUE, left.open = TRUE), 1), nb)
  lo <- X; lo[[feature]] <- z[bin]
  hi <- X; hi[[feature]] <- z[bin + 1]
  diffs <- predict_fun(hi) - predict_fun(lo)
  eff <- numeric(nb)
  m <- tapply(diffs, bin, mean)
  eff[as.integer(names(m))] <- m
  cnt <- tabulate(bin, nbins = nb)
  acc <- c(0, cumsum(eff))
  mid <- (acc[-1] + acc[-length(acc)]) / 2
  acc <- acc - sum(mid * cnt) / sum(cnt)
  tibble(feature = feature, x = z, ale = acc)
}


# ============================================================
# Main loop: one scale, all horizons, folds, predictor sets and models
# ============================================================

run_part <- function(part, scale, fsets, folds) {
  base <- load_scale(scale) %>% trim_batch_end()
  preds <- list(); tuned <- list(); imps <- list()

  for (h in horizons) {
    dh <- build_features(base, h)
    needed <- unique(unlist(lapply(fsets, scale_features, scale = scale)))

    for (fo in folds) {
      # all models are compared on exactly the same rows
      d <- dh %>%
        filter(batch %in% c(fo$train, fo$test)) %>%
        filter(!is.na(y), if_all(all_of(needed), ~ !is.na(.x)))
      train <- d %>% filter(batch %in% fo$train)
      test  <- d %>% filter(batch == fo$test)
      if (nrow(test) == 0 || nrow(train) == 0) next
      message(sprintf("[%s] %s | h=%d | test=%s | n_train=%d n_test=%d",
                      part, scale, h, fo$test, nrow(train), nrow(test)))

      keys <- test %>%
        select(batch, entity, pen, day, date, y, any_of("single_pen"),
               f_mean5, e_temp, eo_temp) %>%
        mutate(part = part, scale = scale, horizon = h, test_batch = fo$test)

      for (m in baseline_models) {
        preds[[length(preds) + 1]] <- keys %>%
          mutate(fset = "none", model = m, pred = predict_baseline(m, test))
      }

      for (fs in names(fsets)) {
        vars <- scale_features(fsets[[fs]], scale)
        for (m in fitted_models) {
          key <- paste(scale, fs, m, fo$test, sep = "|")
          if (is.null(tuned[[key]])) tuned[[key]] <- tune_model(m, train, vars, scale)
          with_q <- m == "rf_rel" && scale != "animal"
          fit <- fit_model(m, train, vars, tuned[[key]], quantiles = with_q)
          out <- keys %>% mutate(fset = fs, model = m, pred = fit$predict(test))
          if (with_q) {
            q <- fit$predict_q(test, interval_probs)
            out <- out %>% mutate(q025 = q[, 1], q05 = q[, 2], q50 = q[, 3],
                                  q95 = q[, 4], q975 = q[, 5])
          }
          preds[[length(preds) + 1]] <- out

          if (m == "rf_rel" && h == h1 && fs == importance_set[[part]]) {
            imps[[length(imps) + 1]] <- perm_importance(fit, test, vars, perm_reps) %>%
              mutate(part = part, scale = scale, test_batch = fo$test, fset = fs)
          }
        }
      }
    }
  }

  tuned_tbl <- bind_rows(lapply(names(tuned), function(k) {
    s <- strsplit(k, "|", fixed = TRUE)[[1]]
    if (length(tuned[[k]]) == 0) return(NULL)
    as_tibble(tuned[[k]]) %>%
      mutate(part = part, scale = s[1], fset = s[2], model = s[3], test_batch = s[4])
  }))
  list(preds = bind_rows(preds), tuned = tuned_tbl, imp = bind_rows(imps))
}


# ============================================================
# Figures: labels, colours and theme
# ============================================================

if (!exists("fig_dpi")) fig_dpi <- 600
if (!exists("fig_pdf")) fig_pdf <- FALSE

save_fig <- function(p, name, w, h) {
  dir.create(out_figures, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_figures, name), p, width = w, height = h, dpi = fig_dpi, bg = "white")
  if (fig_pdf) {
    ggsave(file.path(out_figures, sub("\\.png$", ".pdf", name)), p,
           width = w, height = h, bg = "white",
           device = if (capabilities("cairo")) cairo_pdf else pdf)
  }
  invisible(file.path(out_figures, name))
}

scale_names <- c(animal = "Animal", pen = "Pen", batch = "Batch")

fset_labels <- c(
  "F"      = "Feed history",
  "F+G"    = "+ Growth",
  "F+B"    = "+ Behaviour",
  "F+E"    = "+ Environment",
  "F+E+EO" = "+ Environment (incl. target day)",
  "F+G+B"  = "+ Growth + Behaviour",
  "Full"   = "Full (G + B + E)"
)
fset_label_none <- "No predictors"     # naive references

fsetA_order  <- c("F", "F+G", "F+B", "F+E", "F+E+EO", "Full")
fsetB_order  <- c("F", "F+G", "F+B", "F+G+B")
fsetA_levels <- rev(unname(fset_labels[fsetA_order]))
fsetB_levels <- rev(c(fset_label_none, unname(fset_labels[fsetB_order])))

feature_labels <- c(
  f_lag_a = "Intake, t-h", f_lag_b = "Intake, t-h-1", f_lag_c = "Intake, t-h-2",
  f_mean5 = "5-day mean intake", g_day = "Fattening day", g_weight = "Body weight",
  g_nanimals = "Animals in pen", b_visits = "Feeder visits", b_visits3 = "Visits, 3-day mean",
  b_duration = "Feeder occupation time", b_rate = "Feeding rate",
  e_temp = "Temperature", e_hum = "Humidity", e_co2 = "CO2", e_nh3 = "NH3",
  eo_temp = "Temperature, target day", eo_hum = "Humidity, target day"
)

model_labels <- c(
  persistence = "Persistence",
  mean5       = "5-day mean",
  rf_level    = "RF (kg target)",
  rf_rel      = "RF (relative target)",
  lm          = "Linear regression"
)

# Okabe-Ito colours (colour-blind safe). Models, predictor groups and batches
# never share a colour; batches are told apart by shape only.
model_colours <- c(
  "Persistence"          = "#999999",
  "5-day mean"           = "#555555",
  "RF (kg target)"       = "#E69F00",
  "RF (relative target)" = "#009E73",
  "Linear regression"    = "#CC79A7"
)
model_shapes <- c("Persistence" = 4, "5-day mean" = 3, "RF (kg target)" = 17,
                  "RF (relative target)" = 15, "Linear regression" = 16)

group_colours <- c("Feed history" = "#0072B2", "Growth" = "#D55E00",
                   "Behaviour" = "#56B4E9", "Environment" = "#8856A7")

batch_shapes <- c(batch1 = 16, batch2 = 17, batch3 = 15, batch4 = 18)
batch_labels <- c(batch1 = "Batch 1", batch2 = "Batch 2", batch3 = "Batch 3",
                  batch4 = "Batch 4 (external)")

status_colours <- c("Inside interval" = "grey25", "Below 90% interval" = "#B2182B",
                    "Above 90% interval" = "#4575B4")
status_shapes  <- c("Inside interval" = 16, "Below 90% interval" = 25,
                    "Above 90% interval" = 24)

theme_paper <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"),
        strip.text       = element_text(face = "bold"),
        axis.title.x     = element_text(margin = margin(t = 10)),
        axis.title.y     = element_text(margin = margin(r = 10)),
        axis.text        = element_text(size = rel(0.9), colour = "grey20"),
        panel.spacing.x  = unit(1.4, "lines"),     # stops tick labels of facets touching
        panel.spacing.y  = unit(1.0, "lines"),
        legend.position  = "bottom",
        legend.text      = element_text(size = rel(0.95)),
        legend.key.width = unit(1.1, "lines"),
        plot.margin      = margin(6, 10, 6, 6))

prep_scale <- function(d) d %>% mutate(scale_lab = factor(scale_names[scale], levels = scale_names))
prep_model <- function(d) d %>% mutate(model_label = factor(model_labels[model], levels = model_labels))
