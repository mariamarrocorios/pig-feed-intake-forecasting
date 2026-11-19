# Simulated feeder data with the same structure as the real data used in the paper.
# The real data belong to a commercial farm and cannot be shared, so this script
# creates a fake dataset that lets you run the whole analysis.
# These are NOT the data of the paper and will not give the same results.
#
# Output:
#   data/simulated_visits.csv.gz       one row per visit to the feeder
#   data/simulated_environment.csv     daily climate per pen (batches 1-3 only)

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

set.seed(2026)

n_batches  <- 4
n_pens     <- 10
pigs_range <- 11:13          # pigs per pen at entry
n_days     <- 115
start_dates <- as.Date(c("2025-01-13", "2025-06-02", "2025-11-03", "2026-05-18"))
batch_level <- c(1.00, 1.08, 0.97, 0.95)   # batch 2 eats more, as in the real data

# zones of the climate sensors (each one covers 2 or 3 pens)
zone_of_pen <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4)

# ---- daily intake per pig ----------------------------------------------------

# typical intake curve (kg/day) and body weight (kg) along fattening
intake_curve <- function(d) 1.05 + 1.75 * (1 - exp(-d / 45))
weight_curve <- function(d) 22 + 0.95 * d - 0.0012 * d^2

simulate_pig <- function(batch, pen, pig, start) {
  # each pig has its own level, and daily noise with some memory (AR1)
  pig_effect <- rnorm(1, 1, 0.10)
  noise <- as.numeric(arima.sim(list(ar = 0.3), n = n_days, sd = 0.15))
  d <- seq_len(n_days)
  intake <- intake_curve(d) * batch_level[batch] * pig_effect * (1 + noise)

  # a few pigs get sick: intake drops for some days and may not recover
  if (runif(1) < 0.08) {
    onset <- sample(15:(n_days - 10), 1)
    len <- sample(4:10, 1)
    drop <- seq(0.3, 0.8, length.out = len)
    idx <- onset:min(n_days, onset + len - 1)
    intake[idx] <- intake[idx] * (1 - drop[seq_along(idx)])
    if (runif(1) < 0.4) intake[d > max(idx)] <- NA   # died or removed
  }

  id <- sprintf("B%dP%02dA%02d", batch, pen, pig)
  tibble(batch = paste0("batch", batch), pen = pen, animal_id = id,
         date = start + d - 1, day = d,
         intake_kg = pmax(intake, 0),
         weight_kg = weight_curve(d) * pig_effect^0.5)
}

pig_days <- list()
for (b in seq_len(n_batches)) {
  for (p in seq_len(n_pens)) {
    for (a in seq_len(sample(pigs_range, 1))) {
      pig_days[[length(pig_days) + 1]] <- simulate_pig(b, p, a, start_dates[b])
    }
  }
}
pig_days <- bind_rows(pig_days) %>% filter(!is.na(intake_kg))

# ---- split each day into visits ------------------------------------------------

# number of visits per day goes down as pigs grow; eating rate goes up
visits_mean <- function(d) 11 - 4 * d / n_days
rate_mean   <- function(d) 25 + 30 * d / n_days          # g/min

# visits are more frequent during the day than at night
hour_weights <- c(2, 1, 1, 1, 1, 2, 4, 6, 7, 7, 6, 6, 6, 6, 6, 6, 7, 7, 6, 5, 4, 3, 3, 2)

n_vis <- pmax(1, rpois(nrow(pig_days), visits_mean(pig_days$day)))
low <- pig_days$intake_kg < 0.05
n_vis[low] <- rbinom(sum(low), 1, 0.5)            # some days with no visits at all
row <- rep(seq_len(nrow(pig_days)), n_vis)

share <- rgamma(length(row), 2)
share <- share / ave(share, row, FUN = sum)
feed  <- pig_days$intake_kg[row] * share
rate  <- pmax(rnorm(length(row), rate_mean(pig_days$day[row]), 6), 8)
secs  <- sample(0:86399, length(row), replace = TRUE, prob = rep(hour_weights, each = 3600))

visits <- tibble(batch = pig_days$batch[row], animal_id = pig_days$animal_id[row],
                 pen = pig_days$pen[row],
                 datetime = as.POSIXct(pig_days$date[row], tz = "UTC") + secs,
                 duration_s = round(60 * 1000 * feed / rate),
                 feed_kg = round(feed, 3),
                 weight_kg = round(pig_days$weight_kg[row] + rnorm(length(row), 0, 1.5), 1)) %>%
  filter(duration_s > 0) %>%
  mutate(date = as.Date(datetime), time = format(datetime, "%H:%M:%S", tz = "UTC")) %>%
  select(batch, animal_id, pen, date, time, duration_s, feed_kg, weight_kg) %>%
  arrange(batch, date, time)

# the scale sometimes fails: some weights are missing
visits$weight_kg[runif(nrow(visits)) < 0.05] <- NA

dir.create("data", showWarnings = FALSE)
write_csv(visits, "data/simulated_visits.csv.gz")

# ---- climate inside the barn (batches 1 to 3) ----------------------------------

season_temp <- c(17, 25, 19)      # winter-spring, summer, autumn-winter
env <- list()
for (b in 1:3) {
  for (z in 1:4) {
    d <- seq_len(n_days)
    temp <- season_temp[b] + 3 * sin(2 * pi * d / 60) + rnorm(n_days, 0, 1.2) + (z - 2.5) * 0.4
    hum  <- 65 - 0.8 * (temp - 20) + rnorm(n_days, 0, 4)
    co2  <- 1500 + 8 * d - 40 * (temp - 20) + rnorm(n_days, 0, 150)
    nh3  <- 6 + 0.06 * d + rnorm(n_days, 0, 1.5)
    env[[length(env) + 1]] <- tibble(batch = paste0("batch", b), zone = z,
                                     date = start_dates[b] + d - 1,
                                     temp_c = round(temp, 1),
                                     humidity_pct = round(pmin(pmax(hum, 30), 95), 1),
                                     co2_ppm = round(pmax(co2, 400)),
                                     nh3_ppm = round(pmax(nh3, 0), 1))
  }
}
environment <- bind_rows(env) %>%
  inner_join(tibble(pen = seq_len(n_pens), zone = zone_of_pen), by = "zone",
             relationship = "many-to-many") %>%
  select(batch, pen, date, temp_c, humidity_pct, co2_ppm, nh3_ppm) %>%
  arrange(batch, pen, date)

write_csv(environment, "data/simulated_environment.csv")

message("Simulated ", nrow(visits), " visits of ", n_distinct(visits$animal_id), " pigs.")
