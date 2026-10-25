# Forecasts for the original stationary AR(1)-GARCH(1,1) experiment.
# Run from the repository root after installing the packages in README.md.
# Forecasts are cached in data/; interrupted runs resume from checkpoints.
if (dir.exists('.Rlib')) .libPaths(c(normalizePath('.Rlib'), .libPaths()))
source('Ecomp-Rfns.R')
source('Forecasts.R')
workers <- as.integer(Sys.getenv('FORECAST_WORKERS', '4'))
prepare_data('stationary')
for (distribution in c('norm', 'std', 'sstd'))
  prepare_forecasts('stationary', distribution, workers = workers)
