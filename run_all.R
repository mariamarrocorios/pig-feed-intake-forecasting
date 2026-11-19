# Runs the whole analysis, from the visit records to the tables and figures.
#
# Open R in this folder (or open the .Rproj / set the working directory here)
# and run:
#   source("run_all.R")
#
# With the simulated data a full run takes about 1 hour on 4 cores.
# quick <- TRUE runs a reduced version in a few minutes, only to check that
# everything works.

quick <- FALSE

# Simulated data (included) or your own files with the same columns
visits_file      <- "data/simulated_visits.csv.gz"
environment_file <- "data/simulated_environment.csv"   # NA if you have no climate data

tables_path <- "outputs/tables"
out_path    <- "outputs/results"

# ------------------------------------------------------------------------------

pkgs <- c("readr", "dplyr", "tidyr", "ranger", "ggplot2")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install these packages first: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
if (!file.exists("R/functions.R")) stop("Run this from the folder that contains run_all.R")

t0 <- Sys.time()

if (!file.exists(visits_file) && grepl("simulated", visits_file)) {
  message("Simulating data...")
  source("R/00_simulate_data.R")
}

message("Daily tables...")
source("R/01_daily_tables.R")

message("Forecasting study", if (quick) " (quick mode)", "...")
source("R/02_forecasting_study.R")

message("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min. ",
        "Results in ", out_path, "/ (start with study_summary.txt)")
