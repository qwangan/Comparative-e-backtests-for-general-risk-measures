source("Plotting.R")
library(patchwork)

comparison_heading <- function(key) {
  parts <- strsplit(key,"_",fixed=TRUE)[[1L]]
  risk <- c(VaR="VaR",EXP="Expectile",VaRES="(ES, VaR)")[[parts[1L]]]
  paste0(risk," ",parts[2L]," | ",parts[3L]," vs ",parts[4L])
}

make_publication_figures <- function() {
  dir.create("figures/paper",recursive=TRUE,showWarnings=FALSE)
  stationary <- readRDS("results/stationary/tests.rds")
  chosen <- c("VaR_0.99_n-FHS_t-FHS","EXP_0.99855_t-EVT_st-EVT","VaRES_0.975_n-FP_n-EVT")
  plots <- lapply(chosen,function(key) {
    test <- stationary$comparisons[[key]]$overall
    if(is.null(test)) stop("Missing stationary plot ",key)
    plot_comparison(test,title=comparison_heading(key),show_zones=TRUE)
  })
  combined <- wrap_plots(plots,ncol=1,guides="collect") & theme(legend.position="bottom")
  save_figure(combined,"figures/paper/stationary_paths.pdf",7.2,7.0)
  real <- readRDS("results/nasdaq/tests.rds")
  data <- readRDS("data/nasdaq_data.rds")
  keys <- names(real$comparisons)
  chosen <- c("VaR_0.99_n-EVT_st-EVT","VaR_0.99_st-FP_st-EVT",
               "EXP_0.99855_n-EVT_st-EVT","EXP_0.99855_st-FP_st-EVT",
               "VaRES_0.975_n-EVT_st-EVT","VaRES_0.975_st-FP_st-EVT")
  plots <- lapply(chosen,function(key) {
    test <- real$comparisons[[key]]$resetting
    if(is.null(test)) stop("Missing NASDAQ plot ",key)
    plot_comparison(test,x=data$forecast_dates[test$path$time],
                      title=comparison_heading(key),show_zones=FALSE) +
      geom_vline(xintercept=as.Date(c("2008-09-15","2020-03-11")),
                  colour="#A8A8A8",linetype="dotted",linewidth=.3) +
      labs(subtitle=NULL) + theme(plot.title=element_text(size=9))
  })
  combined <- wrap_plots(plots,ncol=2,guides="collect") & theme(legend.position="bottom")
  save_figure(combined,"figures/paper/nasdaq_resetting.pdf",7.2,7.4)
  # One atlas per risk, with overall and resetting curves for all fifteen pairs.
  for(risk in c("VaR","EXP","VaRES")) {
    relevant <- keys[startsWith(keys,paste0(risk,"_"))]
    for(page in seq_len(ceiling(length(relevant)/3))) {
      selected <- relevant[seq.int((page-1)*3+1,min(page*3,length(relevant)))]
      panels <- list()
      for(key in selected) for(scheme in c("overall","resetting")) {
        test <- real$comparisons[[key]][[scheme]]
        panels[[length(panels)+1L]] <- plot_comparison(test,
          x=data$forecast_dates[test$path$time],title=comparison_heading(key),
          subtitle=if(scheme=="overall") "Overall" else "Resetting",show_zones=FALSE) +
          theme(plot.title=element_text(size=9))
      }
      atlas <- wrap_plots(panels,ncol=2,guides="collect") & theme(legend.position="bottom")
      save_figure(atlas,sprintf("figures/paper/nasdaq_%s_atlas_%02d.pdf",risk,page),7.4,8.1)
    }
  }
  structural <- lapply(c("structural_volatility","structural_shape"),function(dataset) {
    test <- readRDS(file.path("results",dataset,"prespecified_assessments.rds"))$display
    plot_comparison(test,title=if(dataset=="structural_volatility") "n-FHS vs n-EVT" else "st-FP vs st-FHS",
                      subtitle=if(dataset=="structural_volatility") "Volatility coefficient: 0.1 to 0.9" else "Degrees of freedom: 6 to 3",
                      change=2000.5,show_zones=FALSE)
  })
  combined <- wrap_plots(structural,ncol=1,guides="collect") & theme(legend.position="bottom")
  save_figure(combined,"figures/paper/structural_changes.pdf",7.2,5.9)
  invisible(TRUE)
}
