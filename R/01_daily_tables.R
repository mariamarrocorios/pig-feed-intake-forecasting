# From visits to daily tables at three scales: pig, pen and batch.
#
# Input (data/):
#   visits file       one row per visit to the feeder (see README for the columns)
#   environment file  daily climate per pen (optional)
# Output (outputs/tables/):
#   animal_day.csv, pen_day.csv, batch_day.csv, environment_pen_day.csv
#
# To use your own data, export it with the same columns as the simulated file
# and change visits_file / environment_file below (or in run_all.R).
# Visits should already be cleaned (implausible weights set to NA, outlier
# visits removed); see Section 2.2 of the paper for how we did it.

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

if (!exists("visits_file"))      visits_file      <- "data/simulated_visits.csv.gz"
if (!exists("environment_file")) environment_file <- "data/simulated_environment.csv"
if (!exists("tables_path"))      tables_path      <- "outputs/tables"
dir.create(tables_path, recursive = TRUE, showWarnings = FALSE)

visits <- read_csv(visits_file, show_col_types = FALSE) %>%
  mutate(date = as.Date(date),
         rate_g_min = 1000 * feed_kg / (duration_s / 60)) %>%
  group_by(batch) %>%
  mutate(batch_day_fattening = as.numeric(date - min(date)) + 1) %>%
  ungroup()

# ---- pig-day -----------------------------------------------------------------
animal_day <- visits %>%
  group_by(batch, animal_number = animal_id, date) %>%
  summarise(batch_day_fattening = first(batch_day_fattening),
            daily_feed_kg      = sum(feed_kg, na.rm = TRUE),
            daily_visits       = n(),
            daily_duration_sec = sum(duration_s, na.rm = TRUE),
            mean_rate_g_min    = mean(rate_g_min[is.finite(rate_g_min)]),
            mean_weight_kg     = mean(weight_kg, na.rm = TRUE),
            n_pens_day         = n_distinct(pen),
            pens_day           = paste(sort(unique(pen)), collapse = ","),
            .groups = "drop")

# ---- pen-day (means per pig) ----------------------------------------------------
pen_day <- visits %>%
  group_by(batch, location = pen, date) %>%
  summarise(batch_day_fattening = first(batch_day_fattening),
            n_animals_pen_day = n_distinct(animal_id),
            total_feed_kg     = sum(feed_kg, na.rm = TRUE),
            mean_feed_kg      = total_feed_kg / n_animals_pen_day,
            mean_visits       = n() / n_animals_pen_day,
            mean_duration_h   = sum(duration_s, na.rm = TRUE) / 3600 / n_animals_pen_day,
            mean_rate_g_min   = mean(rate_g_min[is.finite(rate_g_min)]),
            mean_weight_kg    = mean(weight_kg, na.rm = TRUE),
            .groups = "drop")

# ---- batch-day (means per pig) -------------------------------------------------
batch_day <- animal_day %>%
  group_by(batch, date) %>%
  summarise(batch_day_fattening = first(batch_day_fattening),
            n_animals_day   = n_distinct(animal_number),
            mean_feed_kg    = mean(daily_feed_kg),
            mean_visits     = mean(daily_visits),
            mean_duration_h = mean(daily_duration_sec) / 3600,
            mean_rate_g_min = mean(mean_rate_g_min, na.rm = TRUE),
            mean_weight_kg  = mean(mean_weight_kg, na.rm = TRUE),
            .groups = "drop")

write_csv(animal_day, file.path(tables_path, "animal_day.csv"))
write_csv(pen_day,    file.path(tables_path, "pen_day.csv"))
write_csv(batch_day,  file.path(tables_path, "batch_day.csv"))

# ---- climate ------------------------------------------------------------------
# If there are no climate data, write an empty table and 02 leaves the
# environment group out.
if (!is.null(environment_file) && !is.na(environment_file) && file.exists(environment_file)) {
  env <- read_csv(environment_file, show_col_types = FALSE) %>%
    transmute(batch, location = pen, date = as.Date(date),
              inside_temp_c_day_mean       = temp_c,
              inside_humidity_pct_day_mean = humidity_pct,
              co2_ppm_day_mean             = co2_ppm,
              nh3_ppm_day_mean             = nh3_ppm)
} else {
  env <- tibble(batch = character(), location = integer(), date = as.Date(character()),
                inside_temp_c_day_mean = numeric(), inside_humidity_pct_day_mean = numeric(),
                co2_ppm_day_mean = numeric(), nh3_ppm_day_mean = numeric())
}
write_csv(env, file.path(tables_path, "environment_pen_day.csv"))

message("Daily tables written to ", tables_path)
