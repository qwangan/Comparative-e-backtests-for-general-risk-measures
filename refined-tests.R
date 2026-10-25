# Refined comparative e-backtests, Sections 3.1--3.4 of Revised.tex.
# No loss range is estimated from the verification sample.
# VaRES forecasts have columns (ES, VaR), in that order.

.refined_log_remainder <- function(u) {
  # log(1+u)-u, with its quadratic behaviour retained near zero.
  ans <- log1p(u) - u
  near <- abs(u) < 1e-4
  if (any(near)) {
    z <- u[near]
    ans[near] <- -z^2/2 + z^3/3 - z^4/4 + z^5/5 - z^6/6
  }
  ans
}

.refined_log_divergence <- function(a, b) {
  # log(a/b)+b/a-1 >= 0.
  u <- (a-b)/b
  pmax(0, .refined_log_remainder(u) + u^2/(1+u))
}

.refined_vector <- function(a, n, name) {
  if (length(a) == 1L) a <- rep(a, n)
  if (length(a) != n || any(!is.finite(a)))
    stop(name, " must be finite, with length one or length(x).")
  as.numeric(a)
}

refined_score <- function(x, internal, standard, p,
                          risk = c("VaR", "EXP", "VaRES"),
                          homogeneous = 1L) {
  risk <- match.arg(risk)
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
    stop("x must be a nonempty finite numeric vector.")
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1)
    stop("p must be in (0,1).")
  if (!(homogeneous %in% c(0L, 1L))) stop("homogeneous must be 0 or 1.")
  if (risk == "EXP" && p < .5) stop("The expectile formulas require p >= 1/2.")
  n <- length(x)
  if (risk == "VaRES") {
    internal <- as.matrix(internal)
    standard <- as.matrix(standard)
    if (nrow(internal) == 1L) internal <- internal[rep(1L, n), , drop=FALSE]
    if (nrow(standard) == 1L) standard <- standard[rep(1L, n), , drop=FALSE]
    if (!identical(dim(internal), c(n, 2L)) ||
        !identical(dim(standard), c(n, 2L)) ||
        any(!is.finite(internal)) || any(!is.finite(standard)))
      stop("VaRES forecasts must be finite n-by-2 matrices (ES, VaR).")
    r <- internal[,1L]; z <- internal[,2L]
    rs <- standard[,1L]; zs <- standard[,2L]
    if (any(r <= 0 | rs <= 0 | z > r | zs >= rs))
      stop("VaRES requires positive ES, internal VaR <= ES and standard VaR < ES.")
    I <- cbind(VaR = as.numeric(x > zs)/(1-p)-1,
               ES = pmax(x-zs, 0)/((1-p)*(rs-zs))-1)
    if (homogeneous == 1L) {
      h <- cbind(VaR = rep(0, n),
                 ES = (1-p)*(rs-zs)*(1/(2*sqrt(r))-1/(2*sqrt(rs))))
      constant <- (1-p)*(sqrt(r)-sqrt(rs))^2
      denominator <- 2*sqrt(r)
      raw <- (pmax(x-z,0)+(1-p)*(r+z))/(2*sqrt(r)) -
        (pmax(x-zs,0)+(1-p)*(rs+zs))/(2*sqrt(rs))
    } else {
      h <- cbind(VaR = rep(0, n), ES = (1-p)*(rs-zs)*(1/r-1/rs))
      constant <- (1-p)*r*.refined_log_divergence(r, rs)
      denominator <- r
      raw <- pmax(x-z,0)/r + (1-p)*(z/r-1+log(r)) -
        pmax(x-zs,0)/rs - (1-p)*(zs/rs-1+log(rs))
    }
    # Difference of hinges, evaluated without subtracting large tail levels.
    dz <- z-zs
    hinge <- ifelse(dz >= 0, -pmin(pmax(x-zs, 0), dz),
                    pmin(pmax(x-z, 0), -dz))
    D <- (constant+hinge+(1-p)*dz)/denominator
    ell <- (constant+pmin((1-p)*dz, -p*dz))/denominator
    u <- (constant+pmax((1-p)*dz, -p*dz))/denominator
  } else {
    r <- .refined_vector(internal, n, "internal")
    rs <- .refined_vector(standard, n, "standard")
    if (homogeneous == 0L && any(r <= 0 | rs <= 0))
      stop("Scores with zero-homogeneous differences require positive forecasts.")
    if (risk == "VaR") {
      I <- cbind(VaR = as.numeric(x > rs)/(1-p)-1)
      h <- matrix(0, n, 1L, dimnames=list(NULL, "VaR"))
      if (homogeneous == 1L) {
        d <- r-rs
        hinge <- ifelse(d >= 0, -pmin(pmax(x-rs,0),d),
                        pmin(pmax(x-r,0),-d))
        D <- (1-p)*d+hinge
        ell <- pmin((1-p)*d, -p*d); u <- pmax((1-p)*d, -p*d)
      } else {
        d <- log(r/rs)
        low <- pmin(r,rs); high <- pmax(r,rs)
        between <- pmin(pmax(x,low),high)
        D <- (1-p)*d + ifelse(r >= rs, -log(between/rs), log(between/r))
        ell <- pmin((1-p)*d, -p*d); u <- pmax((1-p)*d, -p*d)
      }
      raw <- D
    } else {
      weight <- ifelse(x > rs, p, 1-p)
      I <- cbind(EXP = weight*(x-rs))
      if (homogeneous == 1L) {
        d <- r-rs
        h <- cbind(EXP = -2*d)
        # The remainder is constant in both tails and nonnegative everywhere.
        y <- pmin(pmax(x-pmin(r,rs),0),abs(d))
        D <- (1-p)*d^2 + (2*p-1)*ifelse(d >= 0, 2*d*y-y^2, y^2)
        ell <- (1-p)*d^2; u <- p*d^2
        raw <- D+as.numeric(h)*as.numeric(I)
      } else {
        h <- cbind(EXP = 1/r-1/rs)
        constant <- .refined_log_divergence(r,rs)
        between <- pmin(pmax(x,pmin(r,rs)),pmax(r,rs))
        v <- (between-rs)/rs
        extra_up <- .refined_log_remainder(v) + v*(r-rs)/r
        extra_down <- .refined_log_divergence(r,between)
        D <- (1-p)*constant + (2*p-1)*ifelse(r >= rs, extra_up, extra_down)
        ell <- (1-p)*constant; u <- p*constant
        # Removing roundoff at the degenerate equal-report case preserves D=0.
        D <- pmin(pmax(D,ell),u)
        raw <- D+as.numeric(h)*as.numeric(I)
      }
    }
  }
  structure(list(raw=raw, D=D, h=h, I=I, ell=ell, u=u,
                 risk=risk, p=p, homogeneous=homogeneous), class="refined_score")
}

refined_tolerances <- function(I, w=500L, epsilon=0) {
  I <- as.matrix(I)
  n <- nrow(I); d <- ncol(I)
  if (any(!is.finite(I))) stop("Identification errors must be finite.")
  if (length(w) != 1L || w < 1 || w != as.integer(w)) stop("w must be a positive integer.")
  if (length(epsilon) == 1L) epsilon <- rep(epsilon, d)
  if (length(epsilon) != d || any(!is.finite(epsilon) | epsilon < 0))
    stop("epsilon must be a finite nonnegative baseline vector.")
  bias <- matrix(0,n,d,dimnames=dimnames(I))
  if (n > w) {
    cs <- rbind(rep(0,d),apply(I,2L,cumsum))
    if (d == 1L) cs <- matrix(cs,ncol=1L)
    times <- (w+1L):n
    bias[times,] <- (cs[times,,drop=FALSE]-cs[times-w,,drop=FALSE])/w
  }
  baseline <- matrix(rep(epsilon,each=n),n,d,dimnames=dimnames(I))
  list(minus=baseline+pmax(-bias,0),plus=baseline+pmax(bias,0),bias=bias)
}

refined_fdr_level <- function(alpha=.05) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0,1).")
  uniroot(function(a) a*(1-log(a))-alpha,
          c(.Machine$double.xmin,1),tol=.Machine$double.eps)$root
}

.refined_logmeanexp <- function(a) {
  m <- max(a)
  if (is.infinite(m)) return(m)
  m+log(mean(exp(a-m)))
}

.refined_logmixture_path <- function(total, variance, lambda_grid) {
  values <- outer(total,lambda_grid)-outer(variance,lambda_grid^2/8)
  largest <- values[cbind(seq_len(nrow(values)),max.col(values,ties.method="first"))]
  largest+log(rowMeans(exp(values-largest)))
}

comparison_test <- function(parts, start=1L, w=500L, epsilon=0,
                            lambda_grid=2^seq(-6,12), alpha=.05,
                            scheme=c("overall","resetting"),
                            fdr=identical(scheme[1L],"resetting"),
                            breaks=integer(), tolerance=NULL) {
  scheme <- match.arg(scheme)
  n <- length(parts$D); d <- ncol(parts$I)
  if (length(start) != 1L || start < 1 || start > n || start != as.integer(start))
    stop("start must index an observation in parts.")
  if (!length(lambda_grid) || any(!is.finite(lambda_grid) | lambda_grid <= 0))
    stop("lambda_grid must contain positive finite fixed parameters.")
  lambda_grid <- sort(unique(lambda_grid))
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0,1).")
  if (any(!is.finite(c(parts$D,parts$ell,parts$u,parts$h))) ||
      any(parts$ell > parts$u)) stop("Invalid refined scores or bounds.")
  if (any(parts$D < parts$ell-1e-10 | parts$D > parts$u+1e-10))
    stop("A refined score falls outside its bounds.")
  if (is.null(tolerance)) tolerance <- refined_tolerances(parts$I,w,epsilon)
  for (side in c("minus","plus")) {
    tolerance[[side]] <- as.matrix(tolerance[[side]])
    if (!identical(dim(tolerance[[side]]),c(n,d)) ||
        any(!is.finite(tolerance[[side]]) | tolerance[[side]] < 0))
      stop("tolerance must have finite nonnegative n-by-d minus and plus matrices.")
  }
  if (any(!is.finite(breaks)) || any(breaks != as.integer(breaks)))
    stop("breaks must contain observation indices.")
  breaks <- sort(unique(breaks[breaks > start & breaks <= n]))
  if (scheme == "overall" && length(breaks)) stop("Prespecified breaks require a resetting test.")
  hpos <- pmax(parts$h,0); hneg <- pmax(-parts$h,0)
  cminus <- rowSums(hpos*tolerance$minus+hneg*tolerance$plus)
  cplus <- rowSums(hpos*tolerance$plus+hneg*tolerance$minus)
  nominal <- if (scheme == "resetting" && fdr) refined_fdr_level(alpha) else alpha
  log_threshold <- log(2/nominal)
  times <- seq.int(start,n); nt <- length(times)
  path_episode <- rep(1L,nt); path_start <- rep(start,nt)
  path_D <- path_cminus <- path_cplus <- path_variance <- numeric(nt)
  path_lm <- path_lp <- numeric(nt); path_color <- rep("yellow",nt)
  path_report <- rep(FALSE,nt)
  reports <- list(); episode <- 1L; beginning <- start
  sumD <- sumcminus <- sumcplus <- variance <- 0
  allD <- sum(parts$D[times]); allcminus <- sum(cminus[times])
  allcplus <- sum(cplus[times]); allvariance <- sum((parts$u[times]-parts$ell[times])^2)
  if (scheme == "overall") {
    path_D <- cumsum(parts$D[times]); path_cminus <- cumsum(cminus[times])
    path_cplus <- cumsum(cplus[times]); path_variance <- cumsum((parts$u[times]-parts$ell[times])^2)
    path_lm <- .refined_logmixture_path(path_D-path_cminus,path_variance,lambda_grid)
    path_lp <- .refined_logmixture_path(-path_D-path_cplus,path_variance,lambda_grid)
    path_color <- ifelse(path_lm >= log_threshold,"red",
                          ifelse(path_lp >= log_threshold,"green","yellow"))
  } else for (k in seq_along(times)) {
    t <- times[k]
    if (t %in% breaks) {
      # A planned break does not count as a discovery or a false report.
      episode <- episode+1L; beginning <- t
      sumD <- sumcminus <- sumcplus <- variance <- 0
    }
    sumD <- sumD+parts$D[t]; sumcminus <- sumcminus+cminus[t]
    sumcplus <- sumcplus+cplus[t]; variance <- variance+(parts$u[t]-parts$ell[t])^2
    lm <- .refined_logmeanexp(lambda_grid*(sumD-sumcminus)-lambda_grid^2*variance/8)
    lp <- .refined_logmeanexp(lambda_grid*(-sumD-sumcplus)-lambda_grid^2*variance/8)
    color <- if (lm >= log_threshold) "red" else if (lp >= log_threshold) "green" else "yellow"
    path_episode[k] <- episode; path_start[k] <- beginning
    path_D[k] <- sumD; path_cminus[k] <- sumcminus; path_cplus[k] <- sumcplus
    path_variance[k] <- variance; path_lm[k] <- lm; path_lp[k] <- lp
    path_color[k] <- color
    if (scheme == "resetting" && color != "yellow") {
      path_report[k] <- TRUE
      reports[[length(reports)+1L]] <- data.frame(episode=episode,start=beginning,end=t,
        color=color,logMminus=lm,logMplus=lp,meanD=sumD/(t-beginning+1L))
      episode <- episode+1L; beginning <- t+1L
      sumD <- sumcminus <- sumcplus <- variance <- 0
    }
  }
  path <- data.frame(time=times,episode=path_episode,start=path_start,
                     D=parts$D[times],ell=parts$ell[times],u=parts$u[times],
                     cminus=cminus[times],cplus=cplus[times],
                     sumD=path_D,sum_cminus=path_cminus,sum_cplus=path_cplus,
                     variance=path_variance,logMminus=path_lm,logMplus=path_lp,
                     color=path_color,report=path_report)
  # The overall result always uses the complete prescribed assessment period.
  overall_minus <- .refined_logmeanexp(lambda_grid*(allD-allcminus)-lambda_grid^2*allvariance/8)
  overall_plus <- .refined_logmeanexp(lambda_grid*(-allD-allcplus)-lambda_grid^2*allvariance/8)
  overall_color <- if (overall_minus >= log(2/alpha)) "red" else if (overall_plus >= log(2/alpha)) "green" else "yellow"
  overall <- data.frame(start=start,end=n,logMminus=overall_minus,logMplus=overall_plus,
                        color=overall_color,meanD=allD/nt)
  if (!length(reports)) reports <- data.frame(episode=integer(),start=integer(),end=integer(),
      color=character(),logMminus=numeric(),logMplus=numeric(),meanD=numeric())
  else reports <- do.call(rbind,reports)
  list(path=path,reports=reports,overall=overall,parts=parts,tolerance=tolerance,
       config=list(scheme=scheme,start=start,w=w,epsilon=epsilon,lambda_grid=lambda_grid,
                   alpha=alpha,nominal=nominal,log_threshold=log_threshold,
                   threshold=2/nominal,breaks=breaks,
                   uniform_fdr_claim=scheme == "resetting" && fdr && !length(breaks)))
}
