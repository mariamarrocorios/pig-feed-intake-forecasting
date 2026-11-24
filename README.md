# Feed intake forecasting in growing-finishing pigs

[![DOI](https://zenodo.org/badge/1464046293.svg)](https://doi.org/10.5281/zenodo.24367704)

R code for the paper:

> Marroco M, Segalés J, Aguilar-Vega C, Guevara RD, Llonch P, Garcia-Morante B.
> *Effects of aggregation scale, predictor groups and forecast horizon on feed intake
> forecasting in growing-finishing pigs.*
> Smart Agricultural Technology (submitted).

The code forecasts daily feed intake 1, 3 and 7 days ahead with Random Forests (RF)
at three scales (pig, pen and batch) and compares predictor groups: feed intake
history, growth, feeding behaviour and barn environment. RF forecasts are compared
with simple references (yesterday's intake and the mean of the last 5 days) and,
in the supplementary material, with linear regression.

## Data

The data of the paper come from a commercial farm and cannot be shared.
To let anyone run the code, `data/` contains a **simulated** dataset with the same
structure (4 batches, 10 pens, about 480 pigs, 115 days). The simulated data are
not the data of the paper, so the numbers you get will not match the paper.

To use your own data, prepare two files with these columns:

**Feeder visits** (one row per visit)

| column | description |
|---|---|
| batch | batch name, e.g. `batch1` |
| animal_id | pig identifier |
| pen | pen number |
| date | date of the visit (YYYY-MM-DD) |
| time | time of the visit (HH:MM:SS) |
| duration_s | duration of the visit (s) |
| feed_kg | feed eaten in the visit (kg) |
| weight_kg | body weight recorded in the visit (kg), NA if missing |

**Barn environment** (optional, one row per pen and day)

| column | description |
|---|---|
| batch, pen, date | as above |
| temp_c | mean temperature (°C) |
| humidity_pct | mean relative humidity (%) |
| co2_ppm | mean CO2 (ppm) |
| nh3_ppm | mean NH3 (ppm) |

If there are no climate data, set `environment_file <- NA` in `run_all.R` and the
models with the environment group are left out. Visits should already be cleaned
(implausible weights set to NA, outlier visits removed).

The batch names `batch1` to `batch4` are used in the validation: batches 1-3 for
leave-one-batch-out and batch 4 as external test (see the settings at the top of
`R/02_forecasting_study.R`).

## How to run it

R (≥ 4.1) and these packages are needed:

```r
install.packages(c("readr", "dplyr", "tidyr", "ranger", "ggplot2"))
```

Open R in the repository folder and run:

```r
source("run_all.R")
```

A full run takes about 1 hour on 4 cores. With `quick <- TRUE` (top of `run_all.R`)
only the 1-day horizon is run, with fewer trees and smaller tuning grids, in about
5 minutes. The quick mode is only to check that everything works.

## Files

| file | what it does |
|---|---|
| `run_all.R` | runs everything |
| `R/00_simulate_data.R` | creates the simulated dataset in `data/` |
| `R/01_daily_tables.R` | visits → daily tables per pig, pen and batch |
| `R/02_forecasting_study.R` | the study: models, validation, metrics, tests, importance, ALE, prediction intervals, learning curve |
| `R/03_figures.R` | all figures (can be run alone after 02) |
| `R/functions.R` | functions used by 02 and 03 |

Results are written to `outputs/results/`: `tables/` (csv), `figures/` (png),
`predictions/` and `study_summary.txt`, a short text summary of the main results.

## Citation

If you use this code, please cite the paper above and this repository:

> Marroco M. (2026). Feed intake forecasting in growing-finishing pigs (R code).
> Zenodo. https://doi.org/10.5281/zenodo.24367704

## License

MIT, see `LICENSE`.

## Contact

Maria Marroco, IRTA-CReSA
