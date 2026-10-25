# Audit completed result files without refitting any forecast model.
# Run after Run_all.R: Rscript Validate_results.R
source("refined-tests.R")

validate_results <- function(output="results",cache_dir="data",figure_dir="figures",
                             stop_on_failure=TRUE) {
  checks <- list()
  record <- function(test,passed,detail="") {
    checks[[length(checks)+1L]] <<- data.frame(test=test,passed=isTRUE(passed),detail=detail)
    invisible(isTRUE(passed))
  }
  close <- function(a,b,tolerance=1e-8) length(a)==length(b) &&
    all(is.finite(a)) && all(is.finite(b)) && all(abs(a-b)<=tolerance*pmax(1,abs(a),abs(b)))
  read_csv <- function(name) read.csv(file.path(output,name),check.names=FALSE,
                                      stringsAsFactors=FALSE)
  audit_path <- function(test,label) {
    p <- test$path; cfg <- test$config
    required <- c("time","episode","start","D","ell","u","cminus","cplus",
                  "sumD","sum_cminus","sum_cplus","variance","logMminus","logMplus")
    record(paste0(label,": finite paths"),
      all(required %in% names(p)) && all(is.finite(as.matrix(p[,required]))))
    record(paste0(label,": bounds"),all(p$ell<=p$D+1e-10 & p$D<=p$u+1e-10))
    lm <- p$logMminus>=cfg$log_threshold; lp <- p$logMplus>=cfg$log_threshold
    expected_color <- ifelse(lm,"red",ifelse(lp,"green","yellow"))
    record(paste0(label,": exclusive colors"),!any(lm & lp) && identical(p$color,expected_color))
    sums <- function(a) as.numeric(ave(a,p$episode,FUN=cumsum))
    record(paste0(label,": episode sums"),
      close(p$sumD,sums(p$D)) && close(p$sum_cminus,sums(p$cminus)) &&
      close(p$sum_cplus,sums(p$cplus)) && close(p$variance,sums((p$u-p$ell)^2)))
    # Reconstruct fixed-lambda mixtures independently of their sequential update.
    sample_rows <- unique(c(1L,as.integer(nrow(p)/2L),nrow(p),which(p$report)))
    sample_rows <- sample_rows[sample_rows>=1L]
    direct <- function(k,side) {
      a <- if(side=="minus") p$sumD[k]-p$sum_cminus[k] else -p$sumD[k]-p$sum_cplus[k]
      values <- cfg$lambda_grid*a-cfg$lambda_grid^2*p$variance[k]/8
      largest <- max(values)
      largest+log(sum(exp(values-largest))/length(values))
    }
    record(paste0(label,": mixture values"),
      close(p$logMminus[sample_rows],vapply(sample_rows,direct,numeric(1),side="minus")) &&
      close(p$logMplus[sample_rows],vapply(sample_rows,direct,numeric(1),side="plus")))
    training <- seq_len(cfg$start-1L)
    training_scale <- sqrt(mean((test$parts$u[training]-test$parts$ell[training])^2))
    if(!is.finite(training_scale) || training_scale<=1e-12) training_scale <- 1
    record(paste0(label,": training-only grid"),
      close(cfg$lambda_grid,2^seq(-6,12)/training_scale))
    if (test$parts$risk=="EXP") record(paste0(label,": expectile restriction"),
      all(p$D>=-1e-11) && all(p$logMplus<=1e-10) && !any(p$color=="green"))
    if (cfg$scheme=="resetting") {
      r <- test$reports; where <- match(r$end,p$time)
      emitted <- p$time[p$report]
      record(paste0(label,": report times"),identical(as.integer(r$end),as.integer(emitted)) &&
        all(!is.na(where)) && all(r$color==p$color[where]))
      previous <- p$episode==c(NA_integer_,head(p$episode,-1L))
      previous[is.na(previous)] <- FALSE
      prior_crossing <- c(FALSE,head(lm|lp,-1L)) & previous
      record(paste0(label,": first crossings"),!any(prior_crossing) &&
        all(p$report==(lm|lp)) && all(p$start[where]==r$start))
      if(nrow(r)>1L) record(paste0(label,": report intervals"),
        all(tail(r$start,-1L)==head(r$end,-1L)+1L) && all(diff(r$episode)>0))
      record(paste0(label,": FDR threshold"),
        isTRUE(cfg$uniform_fdr_claim) && close(cfg$nominal*(1-log(cfg$nominal)),cfg$alpha) &&
        close(cfg$threshold,2/cfg$nominal))
    } else record(paste0(label,": terminal result"),
      close(tail(p$logMminus,1L),test$overall$logMminus) &&
      close(tail(p$logMplus,1L),test$overall$logMplus) &&
      identical(tail(p$color,1L),test$overall$color))
    # At selected dates directly average the *preceding* errors. Reports do not
    # shorten this window; initialization follows the available training data.
    selected_times <- unique(c(head(p$time,1L),tail(p$time,1L),test$reports$end+1L))
    selected_times <- selected_times[selected_times<=nrow(test$parts$I)]
    baseline <- cfg$epsilon
    if(length(baseline)==1L) baseline <- rep(baseline,ncol(test$parts$I))
    predictable <- vapply(selected_times,function(t) {
      bias <- if(t>cfg$w) colMeans(test$parts$I[(t-cfg$w):(t-1L),,drop=FALSE]) else rep(0,length(baseline))
      close(test$tolerance$minus[t,],baseline+pmax(-bias,0)) &&
        close(test$tolerance$plus[t,],baseline+pmax(bias,0))
    },logical(1))
    record(paste0(label,": past-only calibration"),all(predictable))
  }
  for(dataset in c("stationary","nasdaq","structural_volatility","structural_shape")) {
    tryCatch({
      metadata <- readRDS(file.path(cache_dir,paste0(dataset,"_data.rds")))
      saved <- readRDS(file.path(output,dataset,"tests.rds"))
      overall <- read_csv(file.path(dataset,"overall.csv"))
      expected_models <- if(dataset=="stationary") 10L else 6L
      risk_levels <- if(dataset=="stationary") 9L else 3L
      expected_rows <- expected_models*(expected_models-1L)*risk_levels*4L
      record(paste0(dataset,": all overall comparisons"),nrow(overall)==expected_rows &&
        length(saved$models)==expected_models && !any(overall$internal==overall$standard),
        paste(nrow(overall),"rows; expected",expected_rows))
      expected_paths <- if(dataset=="stationary") 135L else 45L
      record(paste0(dataset,": retained comparison paths"),length(saved$comparisons)==expected_paths)
      expected_levels <- if(dataset=="stationary")
        list(VaR=c(.90,.95,.99),EXP=c(.96561,.98761,.99855),ES=c(.754,.875,.975)) else
        list(VaR=if(dataset=="nasdaq") .99 else .9,EXP=.99855,ES=.975)
      record(paste0(dataset,": original risk levels"),
        all(vapply(names(expected_levels),function(risk)
          close(metadata$levels[[risk]],expected_levels[[risk]]),logical(1))))
      expected_keys <- expand.grid(internal=saved$models,standard=saved$models,
                                  epsilon=saved$config$sensitivity,stringsAsFactors=FALSE)
      expected_keys <- expected_keys[expected_keys$internal!=expected_keys$standard,,drop=FALSE]
      key <- function(a) paste(a$internal,a$standard,a$epsilon,sep="|")
      complete <- TRUE
      for(risk in c("VaR","EXP","VaRES")) {
        levels <- if(risk=="VaRES") metadata$levels$ES else metadata$levels[[risk]]
        for(level in levels) {
          rows <- overall[overall$risk==risk & abs(overall$level-level)<1e-10,,drop=FALSE]
          complete <- complete && identical(sort(key(rows)),sort(key(expected_keys)))
        }
      }
      record(paste0(dataset,": exact model/risk/level coverage"),complete)
      numeric_columns <- vapply(overall,is.numeric,logical(1))
      record(paste0(dataset,": finite overall results"),all(is.finite(as.matrix(overall[,numeric_columns]))))
      record(paste0(dataset,": overall threshold/colors"),all(overall$threshold==40) &&
        !any(overall$logMminus>=log(40) & overall$logMplus>=log(40)) &&
        identical(overall$color,ifelse(overall$logMminus>=log(40),"red",
                              ifelse(overall$logMplus>=log(40),"green","yellow"))))
      grouped <- split(overall,interaction(overall$risk,overall$level,overall$internal,overall$standard,drop=TRUE))
      monotone <- vapply(grouped,function(a) {
        a <- a[order(a$epsilon),]
        all(diff(a$logMminus)<=1e-8) && all(diff(a$logMplus)<=1e-8)
      },logical(1))
      record(paste0(dataset,": fixed-period sensitivity monotone"),all(monotone))
      var_rows <- overall[overall$risk=="VaR",]
      grouped_var <- split(var_rows,interaction(var_rows$level,var_rows$internal,var_rows$standard,drop=TRUE))
      record(paste0(dataset,": VaR independent of tolerance"),all(vapply(grouped_var,function(a)
        diff(range(a$logMminus))<1e-9 && diff(range(a$logMplus))<1e-9,logical(1))))
      exp_rows <- overall[overall$risk=="EXP",]
      record(paste0(dataset,": no expectile green reports"),all(exp_rows$logMplus<=1e-10) &&
        !any(exp_rows$color=="green"))
      invalid_path <- file.path(output,dataset,"invalid_forecasts.csv")
      invalid <- if(file.exists(invalid_path)) read.csv(invalid_path) else data.frame()
      record(paste0(dataset,": no omitted invalid forecasts"),nrow(invalid)==0L,
             if(nrow(invalid)) paste(nrow(invalid),"invalid forecast comparisons") else "")
      record(paste0(dataset,": unchanged sample size"),
        nrow(saved$comparisons[[1L]]$overall$path)==metadata$setup$assessed)
      distributions <- if(dataset=="stationary") c("norm","std","sstd") else c("norm","sstd")
      caches_ok <- vapply(distributions,function(distribution) {
        forecasts <- readRDS(file.path(cache_dir,paste0(dataset,"_",distribution,".rds")))
        fields <- c("VaR","VaRfhs","VaRevt","EXP","EXPfhs","EXPevt","ES","ESfhs","ESevt")
        all(vapply(forecasts[fields],function(a) is.matrix(a) && nrow(a)==length(metadata$loss) &&
                       all(is.finite(a)),logical(1)))
      },logical(1))
      record(paste0(dataset,": complete finite forecast caches"),all(caches_ok))
      if(dataset=="nasdaq") record("nasdaq: original input checksum",
        identical(unname(tools::md5sum("nasdaq0325.csv")),metadata$file_md5))
      for(tag in names(saved$comparisons)) for(scheme in names(saved$comparisons[[tag]]))
        audit_path(saved$comparisons[[tag]][[scheme]],paste(dataset,tag,scheme,sep=" / "))
      if(dataset!="stationary") {
        resetting <- read_csv(file.path(dataset,"resetting_summary.csv"))
        record(paste0(dataset,": all resetting comparisons"),nrow(resetting)==180L &&
          all(resetting$reports==resetting$red+resetting$green) &&
          all(abs(resetting$threshold-2/refined_fdr_level(.05))<1e-7))
        path_exports <- list.files(file.path(output,dataset),pattern="_path\\.csv$")
        report_exports <- list.files(file.path(output,dataset),pattern="_reports\\.csv$")
        record(paste0(dataset,": exported paths and reports"),
          length(path_exports)==45L && length(report_exports)==45L)
      }
      if(grepl("^structural",dataset)) {
        assessments <- readRDS(file.path(output,dataset,"prespecified_assessments.rds"))
        record(paste0(dataset,": two prespecified periods"),
          nrow(assessments$first$path)==2000L && nrow(assessments$second$path)==2000L &&
          tail(assessments$first$path$time,1L)+1L==head(assessments$second$path$time,1L) &&
          close(assessments$second$path$sumD[1L],assessments$second$path$D[1L]) &&
          nrow(assessments$display$reports)==0L)
      }
      if(dataset!="nasdaq") {
        truth <- read_csv(file.path(dataset,"overall_truth.csv"))
        record(paste0(dataset,": complete conditional truth"),nrow(truth)==nrow(overall) &&
          all(is.finite(truth$expected_gap)) && all(is.finite(truth$sign_tolerance)) &&
          identical(truth$color,overall$color))
        record(paste0(dataset,": exact false-report criterion"),
          identical(truth$false_report,(truth$color=="red" & truth$expected_gap<=0) |
            (truth$color=="green" & truth$expected_gap>=0)))
        truth_key <- function(a,reverse=FALSE) paste(a$risk,a$level,
          if(reverse) a$standard else a$internal,if(reverse) a$internal else a$standard,a$epsilon,sep="|")
        matched <- match(truth_key(truth,TRUE),truth_key(truth))
        record(paste0(dataset,": antisymmetric conditional gaps"),all(!is.na(matched)) &&
          close(truth$expected_gap,-truth$expected_gap[matched]))
        record(paste0(dataset,": calibration violation counts"),
          all(truth$calibration_violations>=0 & truth$calibration_violations<=metadata$setup$assessed))
        if(grepl("^structural",dataset)) {
          rt <- read_csv(file.path(dataset,"resetting_truth_summary.csv"))
          record(paste0(dataset,": resetting truth counts"),nrow(rt)==180L &&
            all(rt$false_reports>=0 & rt$false_reports<=rt$reports) &&
            close(rt$realized_fdp,rt$false_reports/pmax(rt$reports,1L)) &&
            all(rt$realized_uniform_fdp>=rt$realized_fdp & rt$realized_uniform_fdp<=1))
          reports_path <- file.path(output,dataset,"resetting_report_truth.csv")
          if(sum(rt$reports)>0) {
            reported <- read.csv(reports_path,stringsAsFactors=FALSE)
            record(paste0(dataset,": every reported interval has truth"),nrow(reported)==sum(rt$reports) &&
              sum(reported$false_report)==sum(rt$false_reports))
          }
        }
      }
      number_figures <- length(list.files(file.path(figure_dir,dataset),pattern="\\.pdf$"))
      minimum_figures <- if(dataset=="stationary") 12L else if(dataset=="nasdaq") 57L else 58L
      record(paste0(dataset,": exported figures"),number_figures>=minimum_figures,
             paste(number_figures,"PDF figures"))
    },error=function(e) record(paste0(dataset,": audit completed"),FALSE,conditionMessage(e)))
  }
  tryCatch({
    traditional <- read_csv("traditional_iid_rejections.csv")
    comparative <- read_csv("comparative_iid_rejections.csv")
    expected <- read_csv("comparative_iid_expected_gaps.csv")
    record("iid: complete replications",nrow(traditional)==24L && nrow(comparative)==120L &&
      nrow(expected)==30L && all(traditional$B==1000L) && all(comparative$B==1000L))
    record("iid: rejection probabilities",all(traditional$rejection>=0 & traditional$rejection<=1) &&
      all(comparative$red>=0 & comparative$green>=0 & comparative$yellow>=0) &&
      close(comparative$red+comparative$green+comparative$yellow,rep(1,nrow(comparative))))
    record("iid: false-report Monte Carlo errors",close(comparative$false_mc_se,
      sqrt(comparative$false_report*(1-comparative$false_report)/comparative$B)))
    record("iid: nominal errors match thresholds",close(comparative$nominal_error,
      pmin(1,2/comparative$threshold)))
    for(name in c("traditional_iid_paths.pdf","traditional_iid_loss.pdf","comparative_iid_red_rates.pdf"))
      record(paste0("iid: figure ",name),file.exists(file.path(figure_dir,name)))
  },error=function(e) record("iid: audit completed",FALSE,conditionMessage(e)))
  tryCatch({
    truth_formula <- read_csv("truth_formula_check.csv")
    truth_checks <- read_csv("truth_checks.csv")
    record("truth: independent analytic and quadrature checks",max(truth_formula$absolute_error)<1e-6)
    record("truth: conditional optima and report reproduction",all(truth_checks$passed))
  },error=function(e) record("truth: audit completed",FALSE,conditionMessage(e)))
  checks <- do.call(rbind,checks)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  write.csv(checks,file.path(output,"validation.csv"),row.names=FALSE)
  failures <- checks[!checks$passed,,drop=FALSE]
  if(nrow(failures)) {
    print(failures,row.names=FALSE)
    if(stop_on_failure) stop(nrow(failures)," result validation checks failed.")
  } else message("All ",nrow(checks)," result validation checks passed.")
  invisible(checks)
}

if(sys.nframe()==0L) validate_results()
