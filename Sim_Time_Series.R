if(dir.exists(".Rlib")) .libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("Analysis.R")
run_dataset("stationary")
