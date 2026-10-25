# Known-DGP conditional score comparisons for the simulated experiments.
# These are realized-path diagnostics, not estimates of an expected FDR.
source("Calibration_diagnostics.R")

student_t_tail_second <- function(q,shape) {
  # Integral_q^Inf x^2 f_shape(x) dx, shape > 2.
  (shape*pt(q,df=shape,lower.tail=FALSE)+
     q*(shape+q^2)*dt(q,df=shape))/(shape-2)
}

fs_t_squared_hinge <- function(q,shape,skew) {
  k <- dt(0,df=shape)
  m <- 2*k*shape/(shape-1)*(skew-1/skew)
  second <- shape/(shape-2)*(skew^2+skew^(-2)-1)
  s <- sqrt(second-m^2)
  y <- m+s*q
  c <- 2/(skew+1/skew)
  right_survival <- c*skew*pt(y/skew,df=shape,lower.tail=FALSE)
  left_survival <- 1-c/skew*pt(skew*y,df=shape)
  survival <- ifelse(y<0,left_survival,right_survival)
  first <- ifelse(y<0,
    m+c/skew^2*(shape+(skew*y)^2)/(shape-1)*dt(skew*y,df=shape),
    c*skew^2*(shape+(y/skew)^2)/(shape-1)*dt(y/skew,df=shape))
  tail_second <- ifelse(y<0,
    second-c/skew^3*(shape/(shape-2)-student_t_tail_second(skew*y,shape)),
    c*skew^3*student_t_tail_second(y/skew,shape))
  pmax(0,(tail_second-2*y*first+y^2*survival)/s^2)
}

true_loss_squared_hinge <- function(data,threshold) {
  q <- (threshold-data$true_mu)/data$true_sigma
  if(data$dataset=="structural_volatility") {
    survival <- pnorm(q,lower.tail=FALSE)
    hinge <- pmax(0,(1+q^2)*survival-q*dnorm(q))
  } else {
    shape <- if(data$dataset=="stationary") 5 else
      ifelse(data$forecast_index<=3000L,6,3)
    skew <- if(data$dataset=="stationary") 1.5 else 1
    hinge <- fs_t_squared_hinge(q,shape,skew)
  }
  data$true_sigma^2*hinge
}

conditional_score <- function(data,forecast,p,risk) {
  if(risk=="VaRES") {
    forecast <- as.matrix(forecast)
    r <- forecast[,1L]; z <- forecast[,2L]
    return((true_loss_tail(data,z)$hinge+(1-p)*(r+z))/(2*sqrt(r)))
  }
  forecast <- as.numeric(forecast)
  if(risk=="VaR") return((1-p)*forecast+true_loss_tail(data,forecast)$hinge)
  (2*p-1)*true_loss_squared_hinge(data,forecast)+
    (1-p)*(forecast^2-2*forecast*data$true_mu)
}

.truth_gauss_legendre <- function(n) {
  k <- seq_len(n-1L); off <- k/sqrt(4*k^2-1)
  jacobi <- matrix(0,n,n)
  jacobi[cbind(k,k+1L)] <- off; jacobi[cbind(k+1L,k)] <- off
  eig <- eigen(jacobi,symmetric=TRUE)
  list(nodes=eig$values,weights=2*eig$vectors[1L,]^2)
}

conditional_refined_expectile <- function(data,r,rs,p,quadrature=32L) {
  rule <- .truth_gauss_legendre(quadrature)
  lo <- pmin(r,rs); length <- abs(r-rs)
  expected <- (1-p)*(r-rs)^2
  for(k in seq_along(rule$nodes)) {
    x <- lo+length*(1+rule$nodes[k])/2
    derivative <- 2*(2*p-1)*ifelse(r>=rs,r-x,x-r)
    expected <- expected+rule$weights[k]*length/2*derivative*true_loss_tail(data,x)$survival
  }
  expected
}

check_simulation_truth <- function() {
  load_forecast_packages()
  rows <- list(); count <- 0L
  for(shape in c(3,5,6,30)) for(skew in c(.6,1,1.5)) for(q in c(-3,-.1,0,1,3,8)) {
    formula <- fs_t_squared_hinge(q,shape,skew)
    integrated <- integrate(function(z) (z-q)^2*fGarch::dsstd(z,nu=shape,xi=skew),
                            q,Inf,rel.tol=1e-9,subdivisions=1000L)$value
    count <- count+1L
    rows[[count]] <- data.frame(check="FS-t squared hinge",shape,skew,q,
                               formula,reference=integrated,absolute_error=abs(formula-integrated))
  }
  set.seed(831)
  for(dataset in c("stationary","structural_shape","structural_volatility")) {
    n <- 60L
    data <- list(dataset=dataset,true_mu=runif(n,-.3,.3),true_sigma=runif(n,.5,2),
                 forecast_index=rep(c(2000L,4000L),length.out=n))
    r <- data$true_mu+data$true_sigma*runif(n,-1,3)
    rs <- data$true_mu+data$true_sigma*runif(n,-1,3)
    p <- .98761
    raw <- conditional_score(data,r,p,"EXP")-conditional_score(data,rs,p,"EXP")
    correction <- -2*(r-rs)*as.numeric(conditional_identification(data,rs,p,"EXP"))
    expected32 <- conditional_refined_expectile(data,r,rs,p,32L)+correction
    expected64 <- conditional_refined_expectile(data,r,rs,p,64L)+correction
    for(name in c("32-point refined identity","64-point refined identity","32/64-point agreement")) {
      error <- switch(name,"32-point refined identity"=abs(raw-expected32),
        "64-point refined identity"=abs(raw-expected64),"32/64-point agreement"=abs(expected32-expected64))
      count <- count+1L
      rows[[count]] <- data.frame(check=paste(dataset,name),shape=NA,skew=NA,q=NA,
        formula=NA,reference=NA,absolute_error=max(error))
    }
  }
  result <- do.call(rbind,rows)
  stopifnot(max(result$absolute_error)<1e-6)
  result
}

.truth_reports <- function(bundle,risk,k,model) {
  a <- bundle$models[[model]]
  if(risk=="VaRES") cbind(ES=a$ES[,k],VaR=a$VaR[,length(bundle$data$levels$VaR)+k])
  else a[[risk]][,k]
}

run_simulation_truth <- function(datasets=c("stationary","structural_volatility","structural_shape"),
                                 output="results",cache_dir="data",verify=TRUE) {
  if(any(!datasets %in% c("stationary","structural_volatility","structural_shape")))
    stop("Known-DGP truth is available only for the simulated time series.")
  if(numeric_config$homogeneous!=1L) stop("This truth calculation uses the original homogeneous scores.")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  if(verify) write.csv(check_simulation_truth(),file.path(output,"truth_formula_check.csv"),row.names=FALSE)
  combined_summary <- checks <- list()
  for(dataset in datasets) {
    bundle <- collect_forecasts(dataset,cache_dir); data <- bundle$data
    saved <- readRDS(file.path(output,dataset,"tests.rds"))
    overall <- read.csv(file.path(output,dataset,"overall.csv"),stringsAsFactors=FALSE)
    indices <- seq.int(data$backtest_start,length(data$loss))
    scores <- expected_I <- list()
    for(risk in c("VaR","EXP","VaRES")) {
      levels <- if(risk=="VaRES") data$levels$ES else data$levels[[risk]]
      for(k in seq_along(levels)) for(model in names(bundle$models)) {
        key <- paste(risk,levels[k],model,sep="|")
        forecast <- .truth_reports(bundle,risk,k,model)
        scores[[key]] <- conditional_score(data,forecast,levels[k],risk)
        expected_I[[key]] <- conditional_identification(data,forecast,levels[k],risk)
      }
    }
    if("opt" %in% names(bundle$models)) for(risk in c("VaR","EXP","VaRES")) {
      levels <- if(risk=="VaRES") data$levels$ES else data$levels[[risk]]
      for(p in levels) {
        opt <- scores[[paste(risk,p,"opt",sep="|")]]
        smallest <- min(vapply(setdiff(names(bundle$models),"opt"),function(model)
          min(scores[[paste(risk,p,model,sep="|")]]-opt),numeric(1)))
        checks[[length(checks)+1L]] <- data.frame(dataset,check=paste(risk,p,"optimal score minimum"),
          value=smallest,passed=smallest>=-1e-8)
        if(smallest< -1e-8) stop("An expected score is below the known conditional optimum.")
      }
    }
    overall$expected_gap <- NA_real_; overall$sign_tolerance <- NA_real_
    overall$calibration_violations <- NA_integer_; overall$max_calibration_excess <- NA_real_
    calibration_cache <- list(); report_rows <- list(); reset_summaries <- list()
    for(row in seq_len(nrow(overall))) {
      risk <- overall$risk[row]; p <- overall$level[row]
      internal <- overall$internal[row]; standard <- overall$standard[row]
      ki <- paste(risk,p,internal,sep="|"); ks <- paste(risk,p,standard,sep="|")
      gap <- scores[[ki]]-scores[[ks]]
      overall$expected_gap[row] <- mean(gap[indices])
      overall$sign_tolerance[row] <- 1e-9*max(1,mean(abs(scores[[ki]][indices])),
                                                        mean(abs(scores[[ks]][indices])))
      cache_key <- paste(ks,overall$epsilon[row],sep="|")
      if(is.null(calibration_cache[[cache_key]])) {
        k <- match(p,if(risk=="VaRES") data$levels$ES else data$levels[[risk]])
        forecasts <- .truth_reports(bundle,risk,k,standard)
        ii <- refined_score(data$loss,forecasts,forecasts,p,risk)$I
        baseline <- baseline_tolerance(data,risk,overall$epsilon[row],p)
        tolerance <- refined_tolerances(ii,numeric_config$window,baseline)
        excess <- pmax(expected_I[[ks]]-tolerance$plus,-expected_I[[ks]]-tolerance$minus,0)
        calibration_cache[[cache_key]] <- list(
          violations=rowSums(excess>1e-8)>0,max_excess=apply(excess,1L,max))
      }
      calibration <- calibration_cache[[cache_key]]
      overall$calibration_violations[row] <- sum(calibration$violations[indices])
      overall$max_calibration_excess[row] <- max(calibration$max_excess[indices])
      # The original resetting comparisons use one orientation for every pair.
      if(dataset!="stationary" && match(internal,saved$models)<match(standard,saved$models)) {
        k <- match(p,if(risk=="VaRES") data$levels$ES else data$levels[[risk]])
        if(overall$epsilon[row]==numeric_config$epsilon) {
          tag <- paste(risk,format(p,trim=TRUE),internal,standard,sep="_")
          test <- saved$comparisons[[tag]]$resetting
        } else {
          parts <- parts_for(bundle,risk,k,internal,standard)
          test <- comparison_test(parts,data$backtest_start,numeric_config$window,
            baseline_tolerance(data,risk,overall$epsilon[row],p),training_grid(parts,data$backtest_start),
            alpha=numeric_config$alpha,scheme="resetting")
        }
        r <- test$reports
        if(nrow(r)) {
          r$expected_gap <- vapply(seq_len(nrow(r)),function(z) mean(gap[r$start[z]:r$end[z]]),numeric(1))
          r$sign_tolerance <- overall$sign_tolerance[row]
          r$false_report <- (r$color=="red" & r$expected_gap<=0) |
            (r$color=="green" & r$expected_gap>=0)
          r$numerically_uncertain_sign <- abs(r$expected_gap)<=r$sign_tolerance
          r$calibration_violations <- vapply(seq_len(nrow(r)),function(z)
            sum(calibration$violations[r$start[z]:r$end[z]]),integer(1))
          r$max_calibration_excess <- vapply(seq_len(nrow(r)),function(z)
            max(calibration$max_excess[r$start[z]:r$end[z]]),numeric(1))
          r$realized_fdp <- cumsum(r$false_report)/seq_len(nrow(r))
          r$dataset <- dataset; r$risk <- risk; r$level <- p
          r$internal <- internal; r$standard <- standard; r$epsilon <- overall$epsilon[row]
          r$conditioning <- "true DGP states and past observations"
          report_rows[[length(report_rows)+1L]] <- r
        }
        reset_summaries[[length(reset_summaries)+1L]] <- data.frame(dataset,risk,level=p,
          internal,standard,epsilon=overall$epsilon[row],reports=nrow(r),
          false_reports=if(nrow(r)) sum(r$false_report) else 0L,
          numerically_uncertain_reports=if(nrow(r)) sum(r$numerically_uncertain_sign) else 0L,
          realized_fdp=if(nrow(r)) tail(r$realized_fdp,1L) else 0,
          realized_uniform_fdp=if(nrow(r)) max(r$realized_fdp) else 0,
          calibration_violations=sum(calibration$violations[indices]),
          note="Realized proportions on one path; not expected FDR estimates")
      }
    }
    overall$false_report <- (overall$color=="red" & overall$expected_gap<=0) |
      (overall$color=="green" & overall$expected_gap>=0)
    overall$numerically_uncertain_sign <- abs(overall$expected_gap)<=overall$sign_tolerance
    overall$true_direction <- ifelse(overall$expected_gap>overall$sign_tolerance,"standard better",
      ifelse(overall$expected_gap< -overall$sign_tolerance,"internal better","numerically indistinguishable"))
    overall$conditioning <- "true DGP states and past observations"
    write.csv(overall,file.path(output,dataset,"overall_truth.csv"),row.names=FALSE)
    if(length(report_rows)) write.csv(do.call(rbind,report_rows),
      file.path(output,dataset,"resetting_report_truth.csv"),row.names=FALSE)
    if(length(reset_summaries)) {
      resets <- do.call(rbind,reset_summaries)
      original <- read.csv(file.path(output,dataset,"resetting_summary.csv"),stringsAsFactors=FALSE)
      key <- function(a) paste(a$risk,a$level,a$internal,a$standard,a$epsilon,sep="|")
      matched <- match(key(resets),key(original))
      identical_counts <- all(!is.na(matched)) && all(resets$reports==original$reports[matched])
      checks[[length(checks)+1L]] <- data.frame(dataset,check="resetting report counts reproduced",
        value=sum(resets$reports),passed=identical_counts)
      if(!identical_counts) stop("Reconstructed resetting counts differ from the saved analysis.")
      write.csv(resets,file.path(output,dataset,"resetting_truth_summary.csv"),row.names=FALSE)
    }
    saveRDS(scores,file.path(output,dataset,"conditional_scores.rds"))
    groups <- split(overall,interaction(overall$risk,overall$level,overall$epsilon,drop=TRUE))
    summaries <- lapply(groups,function(a) data.frame(dataset,risk=a$risk[1],level=a$level[1],
      epsilon=a$epsilon[1],comparisons=nrow(a),reports=sum(a$color!="yellow"),
      false_reports=sum(a$false_report),uncertain_signs=sum(a$numerically_uncertain_sign),
      comparisons_with_calibration_violations=sum(a$calibration_violations>0),
      note="Realized counts across model pairs; not FDR estimates"))
    summary <- do.call(rbind,summaries)
    write.csv(summary,file.path(output,dataset,"overall_truth_summary.csv"),row.names=FALSE)
    combined_summary[[dataset]] <- summary
    message("Simulation truth: ",dataset," complete")
  }
  result <- do.call(rbind,combined_summary)
  write.csv(result,file.path(output,"simulation_truth_summary.csv"),row.names=FALSE)
  if(length(checks)) write.csv(do.call(rbind,checks),file.path(output,"truth_checks.csv"),row.names=FALSE)
  invisible(result)
}

run_betting_diagnostics <- function(datasets=c("stationary","nasdaq","structural_volatility","structural_shape"),
                                    output="results",cache_dir="data",epsilon=.1) {
  best_exponent <- function(A,V) {
    value <- numeric(length(A))
    positive <- A>0 & V>0
    value[positive] <- 2*A[positive]^2/V[positive]
    value[A>0 & V==0] <- Inf
    value
  }
  combined <- list()
  for(dataset in datasets) {
    bundle <- collect_forecasts(dataset,cache_dir)
    saved <- read.csv(file.path(output,dataset,"overall.csv"),stringsAsFactors=FALSE)
    saved <- saved[abs(saved$epsilon-epsilon)<1e-12,,drop=FALSE]
    rows <- lapply(seq_len(nrow(saved)),function(row) {
      risk <- saved$risk[row]; p <- saved$level[row]
      k <- match(p,if(risk=="VaRES") bundle$data$levels$ES else bundle$data$levels[[risk]])
      parts <- parts_for(bundle,risk,k,saved$internal[row],saved$standard[row])
      test <- comparison_test(parts,bundle$data$backtest_start,numeric_config$window,
        baseline_tolerance(bundle$data,risk,epsilon,p),training_grid(parts,bundle$data$backtest_start),
        alpha=numeric_config$alpha,scheme="overall")
      path <- test$path
      minus <- best_exponent(path$sumD-path$sum_cminus,path$variance)
      plus <- best_exponent(-path$sumD-path$sum_cplus,path$variance)
      data.frame(dataset,risk,level=p,internal=saved$internal[row],standard=saved$standard[row],
        epsilon,primary=match(saved$internal[row],names(bundle$models))<match(saved$standard[row],names(bundle$models)),
        terminal_best_fixed_log_e_minus=tail(minus,1L),
        terminal_best_fixed_log_e_plus=tail(plus,1L),
        prefix_best_fixed_log_e_minus=max(minus),prefix_best_fixed_log_e_plus=max(plus),
        terminal_mixture_log_e_minus=tail(path$logMminus,1L),
        terminal_mixture_log_e_plus=tail(path$logMplus,1L),
        prefix_mixture_log_e_minus=max(path$logMminus),prefix_mixture_log_e_plus=max(path$logMplus),
        log_threshold=log(test$config$threshold),
        note="Supremum over fixed lambda is a descriptive upper bound, not a valid selected test")
    })
    result <- do.call(rbind,rows)
    stopifnot(all(result$terminal_mixture_log_e_minus<=result$terminal_best_fixed_log_e_minus+1e-8),
              all(result$terminal_mixture_log_e_plus<=result$terminal_best_fixed_log_e_plus+1e-8),
              all(result$prefix_mixture_log_e_minus<=result$prefix_best_fixed_log_e_minus+1e-8),
              all(result$prefix_mixture_log_e_plus<=result$prefix_best_fixed_log_e_plus+1e-8))
    write.csv(result,file.path(output,dataset,"betting_diagnostics.csv"),row.names=FALSE)
    combined[[dataset]] <- result
    message("Betting diagnostics: ",dataset," complete")
  }
  invisible(do.call(rbind,combined))
}

if(sys.nframe()==0L) {
  arguments <- commandArgs(trailingOnly=TRUE)
  run_simulation_truth(if(length(arguments)) arguments else
    c("stationary","structural_volatility","structural_shape"))
}
