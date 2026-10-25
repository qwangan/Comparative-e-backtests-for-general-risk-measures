# Run with Rscript tests_refined.R. Only base R is required.
source("refined-tests.R")
source("Ecomp-Rfns.R")

near <- function(x,y,tol=1e-10) {
  if (!isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol)))
    stop("Numerical identity failed: maximum difference ",max(abs(x-y)))
}
set.seed(719)
n <- 1200L
x <- rnorm(n,0,7)
r <- runif(n,.4,8); rs <- runif(n,.4,8)
z <- r-runif(n,.1,.4); zs <- rs-runif(n,.1,.4)
for (homogeneous in c(0L,1L)) {
  a <- refined_score(x,r,rs,.95,"VaR",homogeneous)
  near(a$D,sfVaR(r,x,.95,homogeneous)-sfVaR(rs,x,.95,homogeneous))
  b <- refined_score(x,r,rs,.975,"EXP",homogeneous)
  near(b$D,sfEXP(r,x,.975,homogeneous)-sfEXP(rs,x,.975,homogeneous)-
         rowSums(b$h*b$I))
  c <- refined_score(x,cbind(r,z),cbind(rs,zs),.975,"VaRES",homogeneous)
  near(c$D,sfVaRES(z,r,x,.975,homogeneous)-sfVaRES(zs,rs,x,.975,homogeneous)-
         rowSums(c$h*c$I))
  for (parts in list(a,b,c)) {
    stopifnot(all(parts$D >= parts$ell-1e-11),all(parts$D <= parts$u+1e-11))
    # Both tails attain the stated bounds, without a finite loss support.
    if (parts$risk == "VaRES") {
      lo <- refined_score(rep(-1e150,n),cbind(r,z),cbind(rs,zs),.975,"VaRES",homogeneous)
      hi <- refined_score(rep(1e150,n),cbind(r,z),cbind(rs,zs),.975,"VaRES",homogeneous)
    } else {
      lo <- refined_score(rep(-1e150,n),r,rs,parts$p,parts$risk,homogeneous)
      hi <- refined_score(rep(1e150,n),r,rs,parts$p,parts$risk,homogeneous)
    }
    near(pmin(lo$D,hi$D),parts$ell)
    near(pmax(lo$D,hi$D),parts$u)
  }
  # Exact cancellation leaves a nonnegative expectile remainder.
  bb <- comparison_test(b,start=11,w=10,epsilon=.1,lambda_grid=c(.5,2,8))
  stopifnot(all(bb$path$logMplus <= 1e-13),all(b$D >= 0))
}

# Predictability: at t=4 the window uses errors 1,2,3, never error 4.
cal <- refined_tolerances(cbind(a=1:8,b=-(1:8)),w=3,epsilon=c(.1,.2))
near(cal$bias[1:3,],matrix(0,3,2))
near(cal$bias[4,],c(2,-2))
near(cal$minus[4,],c(.1,2.2))
near(cal$plus[4,],c(2.1,.2))
changed <- cbind(a=c(1:3,999,5:8),b=-(1:8))
near(refined_tolerances(changed,3,c(.1,.2))$plus[4,],cal$plus[4,])

# The finite mixture is the mean of whole fixed-parameter processes.
a <- refined_score(x,r,rs,.95,"VaR",1)
grid <- c(.03,.1,.3,1)
aa <- comparison_test(a,start=11,w=10,lambda_grid=grid,alpha=.05)
for (k in c(1,20,nrow(aa$path))) {
  q <- aa$path[k,]
  vals <- grid*q$sumD-grid^2*q$variance/8
  near(q$logMminus,.refined_logmeanexp(vals))
}
stopifnot(!any(aa$path$logMminus >= log(40) & aa$path$logMplus >= log(40)))

# Resetting uses only episode observations, while calibration keeps its history.
constant <- list(D=rep(1,30),ell=rep(1,30),u=rep(1,30),
                 h=matrix(0,30,1),I=matrix(1:30,30,1),raw=rep(1,30))
rr <- comparison_test(constant,w=3,lambda_grid=1,scheme="resetting",alpha=.05)
delay <- ceiling(rr$config$log_threshold)
near(rr$reports$end,seq(delay,30,by=delay))
near(rr$path$sumD[delay+1],1)
near(rr$tolerance$plus[delay+1,1],mean((delay-2):delay))
near(rr$overall$logMminus,30)
stopifnot(rr$config$uniform_fdr_claim)
near(rr$config$nominal*(1-log(rr$config$nominal)),.05)
planned <- comparison_test(constant,w=3,lambda_grid=1,scheme="resetting",alpha=.05,breaks=5)
stopifnot(!planned$config$uniform_fdr_claim)
near(planned$path$sumD[5],1)
near(planned$reports$end[1],5+delay-1)

# Overall tests run to their terminal date after earlier crossings.
reversal <- constant
reversal$D <- c(rep(1,15),rep(-1,15))
reversal$ell <- rep(-1,30); reversal$u <- rep(1,30)
tt <- comparison_test(reversal,w=3,lambda_grid=1,alpha=.2)
stopifnot(any(tt$path$color == "red"),tt$overall$color == "yellow")

# Baseline sensitivity is componentwise monotone, with forecasts and period fixed.
v <- refined_score(x,cbind(r,z),cbind(rs,zs),.975,"VaRES",0)
v0 <- comparison_test(v,start=11,w=10,epsilon=0,lambda_grid=grid)
v1 <- comparison_test(v,start=11,w=10,epsilon=c(.1,.2),lambda_grid=grid)
stopifnot(all(v1$path$logMminus <= v0$path$logMminus+1e-12),
          all(v1$path$logMplus <= v0$path$logMplus+1e-12))
near(v1$parts$D,v0$parts$D)

# Asymmetric calibration corrections account for the sign of h.
manual <- list(D=c(.3,.2),ell=c(-1,-1),u=c(1,1),
               h=matrix(c(2,-3),2,1),I=matrix(0,2,1),raw=c(.3,.2))
tol <- list(minus=matrix(c(.1,.2),2,1),plus=matrix(c(.4,.5),2,1))
mt <- comparison_test(manual,w=1,lambda_grid=1,tolerance=tol)
near(mt$path$cminus,c(.2,1.5)); near(mt$path$cplus,c(.8,.6))

# Close forecasts retain a nonnegative quadratic logarithmic remainder.
nearparts <- refined_score(c(-1,1,3),rep(1+1e-9,3),rep(1,3),.99,"EXP",0)
stopifnot(all(nearparts$D >= 0),all(nearparts$u > 0))
cat("All refined comparative-test checks passed.\n")
