# Install once, explicitly. Experiment scripts never install software.
local_library <- file.path(getwd(),".Rlib")
dir.create(local_library,showWarnings=FALSE)
.libPaths(c(local_library,.libPaths()))
packages <- c("rugarch","fGarch","expectreg","ismev","skewt","sgt","ggplot2","patchwork")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org",lib=local_library)
cat("Numerical-study dependencies are ready.\n")
