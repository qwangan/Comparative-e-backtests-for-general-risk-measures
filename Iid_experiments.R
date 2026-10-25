# The original iid data-generating settings; terminal comparisons are explicit.
source("refined-tests.R")
source("Plotting.R")

run_traditional_iid <- function(B=1000L,n=1000L,l=10L,seed=20251007L,
                                 output="results") {
  set.seed(seed)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  scenarios <- data.frame(scenario=c("-5% VaR","-10% VaR","-5% ES","-10% ES","-5% Both","Baseline"),
                          k=c(.95,.9,1,1,.95,1),es=c(1,1,.95,.9,.95,1))
  values <- c(0,.1*(1:5),-.1*(1:5)); Tlen <- n+l; p <- .95
  thresholds <- c(2,5,10,20); rates <- paths <- list()
  for (s in seq_len(nrow(scenarios))) {
    logM <- matrix(0,Tlen,B)
    for (b in seq_len(B)) {
      eps <- sample(values,Tlen,replace=TRUE)
      z <- 1.64*scenarios$k[s]+eps
      r <- 2.06*scenarios$es[s]+eps
      loss <- rnorm(Tlen)
      # There are only eleven possible reports. Cumulative sums evaluate the
      # original GREL rule at today's report without rescoring the whole past.
      zvalues <- 1.64*scenarios$k[s]+values
      gap <- 2.06*scenarios$es[s]-1.64*scenarios$k[s]
      if (gap<=0) stop("ES must exceed VaR in the original iid design.")
      errors <- vapply(zvalues,function(z0) pmax(loss-z0,0)/((1-p)*gap)-1,
                        numeric(Tlen))
      sums <- apply(errors,2,cumsum); squares <- apply(errors^2,2,cumsum)
      times <- (l+1L):Tlen
      locations <- cbind(times-1L,match(eps[ times ],values))
      den <- squares[locations]
      lambda <- pmin(.5,pmax(0,ifelse(den>0,sums[locations]/den,0)))
      current <- errors[cbind(times,match(eps[times],values))]
      logM[times,b] <- cumsum(log1p(lambda*current))
      if (s==1L && b==1L) example <- data.frame(time=seq_len(Tlen),loss=loss,z=z,r=r)
    }
    rates[[s]] <- data.frame(scenario=scenarios$scenario[s],threshold=thresholds,
      rejection=vapply(thresholds,function(a) mean(logM[Tlen,]>=log(a)),numeric(1)),B=B)
    paths[[s]] <- data.frame(time=seq_len(Tlen),scenario=scenarios$scenario[s],
                              mean_log=rowMeans(logM),
                              log_mean=apply(logM,1,function(a) .refined_logmeanexp(a)))
    message("Traditional iid: ",scenarios$scenario[s]," complete")
  }
  rates <- do.call(rbind,rates); paths <- do.call(rbind,paths)
  rates$mc_se <- sqrt(rates$rejection*(1-rates$rejection)/B)
  write.csv(rates,file.path(output,"traditional_iid_rejections.csv"),row.names=FALSE)
  write.csv(paths,file.path(output,"traditional_iid_paths.csv"),row.names=FALSE)
  write.csv(example,file.path(output,"traditional_iid_example.csv"),row.names=FALSE)
  paths$scenario <- factor(paths$scenario,levels=c("Baseline",scenarios$scenario[1:5]))
  palette <- c("Baseline"="#444444","-5% VaR"="#B55F46","-10% VaR"="#D9A02E",
               "-5% ES"="#437FAD","-10% ES"="#6C5898","-5% Both"="#2B8E79")
  fig <- ggplot(paths,aes(time,mean_log/log(10),colour=scenario)) +
    geom_line(linewidth=.55) + scale_colour_manual(values=palette) +
    labs(x="Observation",y=expression("Mean "*log[10]~"e-value")) +
    paper_theme() + theme(legend.position="bottom") +
    guides(colour=guide_legend(nrow=2))
  save_figure(fig,"figures/traditional_iid_paths.pdf",6.8,3.7)
  ex <- rbind(data.frame(time=example$time,value=example$loss,series="Loss"),
              data.frame(time=example$time,value=example$z,series="VaR forecast"),
              data.frame(time=example$time,value=example$r,series="ES forecast"))
  fig <- ggplot(ex,aes(time,value,colour=series)) + geom_line(linewidth=.35) +
    scale_colour_manual(values=c("Loss"="#9E9E9E","VaR forecast"="#237B65","ES forecast"="#B34440")) +
    labs(x="Observation",y="Loss and forecast") + paper_theme() +
    theme(legend.position="bottom")
  save_figure(fig,"figures/traditional_iid_loss.pdf",6.8,3.7)
  invisible(list(rejections=rates,paths=paths,example=example))
}

normal_pinball_mean <- function(r,p=.95) (1-p)*r+dnorm(r)-r*pnorm(r,lower.tail=FALSE)

run_comparative_iid <- function(B=1000L,n=1000L,l=10L,output="results") {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  # The original noise was sampled before set.seed(233); seed it separately.
  set.seed(20251007)
  eps <- sample(-5:5,n+l,replace=TRUE)/20
  truth <- qnorm(.95)
  names <- c("-50%","-20%","Noisy","+20%","+50%","Exact")
  forecasts <- outer(rep(truth,n+l),c(.5,.8,1,1.2,1.5),"*")+eps
  forecasts <- cbind(forecasts,rep(truth,n+l)); colnames(forecasts) <- names
  set.seed(233)
  y <- matrix(rnorm((n+l)*B),nrow=B,ncol=n+l)
  grid <- 2^seq(-6,12)
  thresholds <- c(2,5,10,40)
  rows <- list(); rates <- list(); z <- 0L
  for (i in seq_along(names)) for (j in setdiff(seq_along(names),i)) {
    D <- vapply(seq_len(B),function(b)
      refined_score(y[b,],forecasts[,i],forecasts[,j],.95,"VaR")$D,numeric(n+l))
    width <- abs(forecasts[,i]-forecasts[,j])
    scale <- sqrt(mean(width[seq_len(l)]^2))
    lambdas <- grid/max(scale,1e-12)
    index <- (l+1L):(n+l)
    sumD <- colSums(D[index,,drop=FALSE]); variance <- sum(width[index]^2)/8
    lm <- vapply(sumD,function(a) .refined_logmeanexp(lambdas*a-lambdas^2*variance),numeric(1))
    lp <- vapply(sumD,function(a) .refined_logmeanexp(-lambdas*a-lambdas^2*variance),numeric(1))
    delta <- mean(normal_pinball_mean(forecasts[index,i])-normal_pinball_mean(forecasts[index,j]))
    z <- z+1L
    rows[[z]] <- data.frame(internal=names[i],standard=names[j],expected_gap=delta,
                            mean_logMminus=mean(lm),mean_logMplus=mean(lp))
    rates[[z]] <- do.call(rbind,lapply(thresholds,function(a) {
      red <- lm>=log(a); green <- lp>=log(a)
      data.frame(internal=names[i],standard=names[j],threshold=a,expected_gap=delta,
                  red=mean(red),green=mean(green),yellow=mean(!red & !green),
                  false_report=mean((red & delta<=0)|(green & delta>=0)),
                  nominal_error=min(1,2/a),B=B)
    }))
    message("Comparative iid: ",names[i]," versus ",names[j])
  }
  rates <- do.call(rbind,rates); rows <- do.call(rbind,rows)
  rates$false_mc_se <- sqrt(rates$false_report*(1-rates$false_report)/B)
  write.csv(rates,file.path(output,"comparative_iid_rejections.csv"),row.names=FALSE)
  write.csv(rows,file.path(output,"comparative_iid_expected_gaps.csv"),row.names=FALSE)
  write.csv(data.frame(time=seq_len(n+l),forecasts),file.path(output,"comparative_iid_forecasts.csv"),row.names=FALSE)
  cells <- subset(rates,threshold==40)
  cells$color <- ifelse(cells$red>cells$green,"red","green")
  fig <- ggplot(cells,aes(internal,standard,fill=red)) +
    geom_tile(colour="white",linewidth=.5) +
    geom_text(aes(label=sprintf("%.2f",red)),size=3) +
    scale_fill_gradient(low="#FFF8F1",high="#C35A54",limits=c(0,1),name="Red rate") +
    coord_equal() + labs(x="Internal model",y="Standard model") + paper_theme() +
    theme(panel.grid=element_blank(),legend.position="right")
  save_figure(fig,"figures/comparative_iid_red_rates.pdf",5.7,4.5)
  invisible(list(rejections=rates,expected_gaps=rows))
}
