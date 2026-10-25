# Usage: Rscript Run_forecasts.R [dataset|all] [workers]
# Optional: FORECAST_MAX_WINDOWS=12 for a resumable pilot.
local_library <- file.path(getwd(), ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
source("Ecomp-Rfns.R")
source("Forecasts.R")
arguments <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(arguments)) arguments[1] else "all"
workers <- if (length(arguments) >= 2) as.integer(arguments[2]) else 4L
datasets <- if (dataset == "all") c("stationary", "nasdaq", "structural_volatility", "structural_shape") else dataset
max_windows <- as.numeric(Sys.getenv("FORECAST_MAX_WINDOWS", "Inf"))
for (name in datasets) {
  prepare_data(name)
  distributions <- if (name == "stationary") c("norm", "std", "sstd") else c("norm", "sstd")
  for (distribution in distributions)
    prepare_forecasts(name, distribution, workers = workers, max_windows = max_windows)
}
