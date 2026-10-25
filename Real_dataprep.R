# Risk forecasts for the unchanged NASDAQ sample, January 2003--May 2025.
# The first 500 returns fit the models; the next 500 forecasts initialize
# calibration and betting parameters before the backtesting period.
if (dir.exists('.Rlib')) .libPaths(c(normalizePath('.Rlib'), .libPaths()))
source('Ecomp-Rfns.R')
source('Forecasts.R')
workers <- as.integer(Sys.getenv('FORECAST_WORKERS', '4'))
prepare_data('nasdaq')
for (distribution in c('norm', 'sstd'))
  prepare_forecasts('nasdaq', distribution, workers = workers)
