# Regenerate figures from completed results without refitting forecasts.
if(dir.exists(".Rlib")) .libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("Analysis.R")
source("Publication_figures.R")
for(dataset in c("stationary","nasdaq","structural_volatility","structural_shape")) {
  overall <- read.csv(file.path("results",dataset,"overall.csv"))
  models <- readRDS(file.path("results",dataset,"tests.rds"))$models
  save_overall_heatmaps(overall,models,dataset)
}
make_publication_figures()
