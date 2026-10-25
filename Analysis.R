# End-to-end analyses using the original losses and forecast models.
source("Forecasts.R")
source("refined-tests.R")
source("Plotting.R")

numeric_config <- list(window=500L,epsilon=.1,sensitivity=c(0,.1,.25,.5),
                       lambda_grid=2^seq(-6,12),alpha=.05,homogeneous=1L,
                       thresholds=c(2,5,10,40))

collect_forecasts <- function(dataset,cache_dir="data") {
  load_forecast_packages()
  data <- readRDS(file.path(cache_dir,paste0(dataset,"_data.rds")))
  distributions <- if(dataset=="stationary") c("norm","std","sstd") else c("norm","sstd")
  labels <- c(norm="n",std="t",sstd="st")
  models <- list()
  for(distribution in distributions) {
    out <- readRDS(file.path(cache_dir,paste0(dataset,"_",distribution,".rds")))
    for(method in c("FP","FHS","EVT")) {
      suffix <- c(FP="",FHS="fhs",EVT="evt")[[method]]
      models[[paste(labels[distribution],method,sep="-")]] <- list(
        VaR=out[[paste0("VaR",suffix)]],EXP=out[[paste0("EXP",suffix)]],
        ES=out[[paste0("ES",suffix)]])
    }
  }
  if(dataset=="stationary") {
    nu <- 5; xi <- 1.5
    mst <- mean.st(nu,xi); sst <- sqrt(var.st(nu,xi))
    q <- fGarch::qsstd(c(data$levels$VaR,data$levels$ES),nu=nu,xi=xi)
    e <- (est(data$levels$EXP,shape=nu,skew=xi)-mst)/sst
    es <- vapply(data$levels$ES,function(p) {
      q0 <- fGarch::qsstd(p,nu=nu,xi=xi)
      integrate(function(x) x*fGarch::dsstd(x,nu=nu,xi=xi),q0,Inf)$value/(1-p)
    },numeric(1))
    transform <- function(a) outer(data$true_sigma,a,"*")+data$true_mu
    models$opt <- list(VaR=transform(q),EXP=transform(e),ES=transform(es))
  }
  list(data=data,models=models)
}

training_grid <- function(parts,start,grid=numeric_config$lambda_grid) {
  before <- seq_len(start-1L)
  scale <- sqrt(mean((parts$u[before]-parts$ell[before])^2))
  # Identical forecasts yield zero capital increments for every lambda.
  grid/ifelse(is.finite(scale)&&scale>1e-12,scale,1)
}

baseline_tolerance <- function(data,risk,multiplier,p=NULL) {
  if(risk=="EXP" && (is.null(p)||length(p)!=1L)) stop("Supply the expectile level p.")
  scale <- if(risk=="EXP") (1-p)*sd(data$loss[seq_len(data$backtest_start-1L)]) else 1
  rep(multiplier*scale,if(risk=="VaRES") 2L else 1L)
}

parts_for <- function(bundle,risk,k,i,j) {
  reports <- function(model) {
    a <- bundle$models[[model]]
    if(risk=="VaRES") cbind(ES=a$ES[,k],VaR=a$VaR[,length(bundle$data$levels$VaR)+k])
    else a[[risk]][,k]
  }
  level <- if(risk=="VaRES") bundle$data$levels$ES[k] else bundle$data$levels[[risk]][k]
  refined_score(bundle$data$loss,reports(i),reports(j),level,risk,
                  homogeneous=numeric_config$homogeneous)
}

save_overall_heatmaps <- function(overall_rows,models,dataset) {
  for(risk in c("VaR","EXP","VaRES")) for(threshold in numeric_config$thresholds) {
    # Base R has no tidy-evaluation operators: explicit indexing avoids capture.
    cells <- overall_rows[overall_rows$risk==risk & overall_rows$epsilon==numeric_config$epsilon,,drop=FALSE]
    shown_levels <- sort(unique(cells$level))
    if(length(shown_levels)>2) shown_levels <- shown_levels[c(1,length(shown_levels))]
    cells <- cells[cells$level %in% shown_levels,,drop=FALSE]
    cells$color <- ifelse(cells$logMminus>=log(threshold),"red",
                          ifelse(cells$logMplus>=log(threshold),"green","yellow"))
    diagonals <- expand.grid(internal=models,level=shown_levels,stringsAsFactors=FALSE)
    diagonals$standard <- diagonals$internal; diagonals$color <- "diagonal"
    cells <- rbind(cells[,c("internal","standard","level","color")],diagonals[,c("internal","standard","level","color")])
    cells$internal <- factor(cells$internal,levels=models)
    cells$standard <- factor(cells$standard,levels=models)
    fig <- plot_heatmap(cells,title=sprintf("%s | terminal threshold %g",risk,threshold))
    save_figure(fig,file.path("figures",dataset,paste0(risk,"_overall_",threshold,".pdf")),
                  if(length(shown_levels)>1) 8 else 5.2,4.7)
  }
}

run_dataset <- function(dataset,output="results",cache_dir="data") {
  bundle <- collect_forecasts(dataset,cache_dir)
  data <- bundle$data; models <- names(bundle$models); start <- data$backtest_start
  dir.create(file.path(output,dataset),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path("figures",dataset),recursive=TRUE,showWarnings=FALSE)
  overall_rows <- reset_rows <- score_rows <- invalid <- list()
  comparisons <- list(); index <- 0L; invalid_index <- 0L
  for(risk in c("VaR","EXP","VaRES")) {
    levels <- if(risk=="VaRES") data$levels$ES else data$levels[[risk]]
    for(k in seq_along(levels)) for(i in models) for(j in setdiff(models,i)) {
      # Both orientations are retained: the standard model determines calibration.
      parts <- tryCatch(parts_for(bundle,risk,k,i,j),error=function(e)e)
      if(inherits(parts,"error")) {
        invalid_index <- invalid_index+1L
        invalid[[invalid_index]] <- data.frame(dataset,risk,level=levels[k],internal=i,
                                              standard=j,reason=conditionMessage(parts))
        next
      }
      grid <- training_grid(parts,start)
      for(epsilon in numeric_config$sensitivity) {
        baseline <- baseline_tolerance(data,risk,epsilon,levels[k])
        overall <- comparison_test(parts,start,numeric_config$window,baseline,grid,
                                   alpha=numeric_config$alpha,scheme="overall")
        index <- index+1L
        terminal <- tail(overall$path,1)
        overall_rows[[index]] <- data.frame(dataset,risk,level=levels[k],internal=i,standard=j,
          epsilon,epsilon_components=paste(signif(baseline,8),collapse=";"),
          logMminus=terminal$logMminus,logMplus=terminal$logMplus,
          threshold=overall$config$threshold,color=terminal$color,
          mean_raw=mean(parts$raw[start:length(parts$raw)]),
          mean_D=mean(parts$D[start:length(parts$D)]),
          mean_cminus=mean(terminal$sum_cminus/nrow(overall$path)),
          mean_cplus=mean(terminal$sum_cplus/nrow(overall$path)),
          lambda_min=min(grid),lambda_max=max(grid),lambda_count=length(grid))
        # Match the original fifteen real-data comparisons and all structural
        # pairs. Stationary heatmaps are overall terminal comparisons.
        primary <- match(i,models)<match(j,models)
        if(primary && dataset!="stationary") {
          reset <- comparison_test(parts,start,numeric_config$window,baseline,grid,
                                     alpha=numeric_config$alpha,scheme="resetting")
          reset_rows[[length(reset_rows)+1L]] <- data.frame(dataset,risk,level=levels[k],
            internal=i,standard=j,epsilon,threshold=reset$config$threshold,
            reports=nrow(reset$reports),red=sum(reset$reports$color=="red"),
            green=sum(reset$reports$color=="green"))
          if(epsilon==numeric_config$epsilon) {
            tag <- paste(risk,format(levels[k],trim=TRUE),i,j,sep="_")
            comparisons[[tag]] <- list(overall=overall,resetting=reset)
            export_path <- reset$path
            export_path$date <- as.character(data$forecast_dates[export_path$time])
            write.csv(export_path,file.path(output,dataset,paste0(tag,"_path.csv")),row.names=FALSE)
            reports <- reset$reports
            if(nrow(reports)) {
              reports$start_date <- as.character(data$forecast_dates[reports$start])
              reports$report_date <- as.character(data$forecast_dates[reports$end])
            }
            write.csv(reports,file.path(output,dataset,paste0(tag,"_reports.csv")),row.names=FALSE)
            plot <- plot_comparison(reset,x=data$forecast_dates[reset$path$time],
                      title=paste(i,"vs",j),
                      subtitle=sprintf("Resetting | %s %.5g | baseline %.2g",risk,levels[k],epsilon))
            if(dataset=="nasdaq") plot <- plot+
              geom_vline(xintercept=as.Date(c("2008-09-15","2020-03-11")),
                          colour="#B3B3B3",linetype="dotted",linewidth=.35)
            save_figure(plot,file.path("figures",dataset,paste0(tag,"_resetting.pdf")),6.6,3.3)
          }
        }
        if(epsilon==numeric_config$epsilon && dataset=="stationary" && primary && k==length(levels)) {
          tag <- paste(risk,format(levels[k],trim=TRUE),i,j,sep="_")
          # Retain full trajectories for all original stationary pairwise tests.
          comparisons[[tag]] <- list(overall=overall)
        }
      }
      if(match(j,models)==length(models)) message(dataset,": ",risk," ",levels[k]," / ",i)
    }
  }
  overall_rows <- do.call(rbind,overall_rows)
  write.csv(overall_rows,file.path(output,dataset,"overall.csv"),row.names=FALSE)
  if(length(reset_rows)) write.csv(do.call(rbind,reset_rows),file.path(output,dataset,"resetting_summary.csv"),row.names=FALSE)
  if(length(invalid)) write.csv(do.call(rbind,invalid),file.path(output,dataset,"invalid_forecasts.csv"),row.names=FALSE)
  saveRDS(list(comparisons=comparisons,config=numeric_config,models=models,
               data_setup=data$setup),file.path(output,dataset,"tests.rds"))
  save_overall_heatmaps(overall_rows,models,dataset)
  invisible(list(overall=overall_rows,comparisons=comparisons))
}

run_structural_assessments <- function(dataset,output="results") {
  bundle <- collect_forecasts(dataset)
  start <- bundle$data$backtest_start
  change <- start+2000L
  chosen <- if(dataset=="structural_volatility") c("n-FHS","n-EVT") else c("st-FP","st-FHS")
  parts <- parts_for(bundle,"VaR",1L,chosen[1],chosen[2])
  grid <- training_grid(parts,start)
  baseline <- baseline_tolerance(bundle$data,"VaR",numeric_config$epsilon)
  # Two prespecified assessments, not an extra reset of a discovery stream.
  # No threshold-triggered restarts are used in the displayed two assessments.
  # Construct their uninterrupted overall paths separately.
  first_parts <- lapply(parts,function(x) x)
  ids1 <- seq_len(change-1L)
  for(name in c("raw","D","ell","u")) first_parts[[name]] <- parts[[name]][ids1]
  for(name in c("h","I")) first_parts[[name]] <- parts[[name]][ids1,,drop=FALSE]
  first <- comparison_test(first_parts,start,500L,baseline,grid,alpha=.05,scheme="overall")
  second <- comparison_test(parts,change,500L,baseline,grid,alpha=.05,scheme="overall")
  second$path$episode <- 2L
  display <- first; display$path <- rbind(first$path,second$path)
  display$reports <- display$reports[0,,drop=FALSE]
  path <- display$path; path$time <- path$time-start+1L
  display$path <- path
  saveRDS(list(first=first,second=second,display=display),file.path(output,dataset,"prespecified_assessments.rds"))
  write.csv(display$path,file.path(output,dataset,"prespecified_assessments.csv"),row.names=FALSE)
  fig <- plot_comparison(display,title=paste(chosen,collapse=" vs "),
                          subtitle="Two prespecified assessments | VaR 0.90",change=2000.5)
  save_figure(fig,file.path("figures",dataset,"prespecified_assessments.pdf"),6.6,3.5)
  invisible(display)
}
