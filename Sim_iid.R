if(dir.exists(".Rlib")) .libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("Iid_experiments.R")
run_traditional_iid()
