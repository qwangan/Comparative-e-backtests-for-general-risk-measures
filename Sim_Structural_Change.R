if(dir.exists(".Rlib")) .libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("Analysis.R")
for(dataset in c("structural_volatility","structural_shape")) {
  prepare_data(dataset)
  for(distribution in c("norm","sstd")) prepare_forecasts(dataset,distribution,workers=4L)
  run_dataset(dataset)
  run_structural_assessments(dataset)
}
