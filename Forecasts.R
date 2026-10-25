# Reproducible forecast preparation for all numerical experiments.
# Run from this repository; no packages are installed by these functions.

forecast_packages <- c("rugarch", "fGarch", "expectreg", "ismev", "skewt", "sgt")

load_forecast_packages <- function() {
  missing <- forecast_packages[!vapply(forecast_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  for (package in forecast_packages)
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  if (!exists("RM.forecasts3", mode = "function")) source("Ecomp-Rfns.R")
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not save ", path)
}

prepare_data <- function(dataset, cache_dir = "data", force = FALSE) {
  dataset <- match.arg(dataset, c("stationary", "nasdaq", "structural_volatility",
                                  "structural_shape"))
  path <- file.path(cache_dir, paste0(dataset, "_data.rds"))
  if (file.exists(path) && !force) return(readRDS(path))
  load_forecast_packages()
  w <- 500L
  if (dataset == "stationary") {
    set.seed(2345)
    spec <- fGarch::garchSpec(model = list(mu = -.05, ar = .3, omega = .01,
                                          alpha = .1, beta = .85, skew = 1.5,
                                          shape = 5), cond.dist = "sstd")
    simulation <- fGarch::garchSim(spec, n = 6000, n.start = 1000, extended = TRUE)
    x <- as.numeric(simulation$garch)
    true_sigma <- as.numeric(simulation$sigma)
    true_mu <- x - true_sigma * as.numeric(simulation$eps)
    levels <- list(VaR = c(.90, .95, .99), EXP = c(.96561, .98761, .99855),
                   ES = c(.754, .875, .975))
    setup <- list(seed = 2345L, burn_in = 1000L, assessed = 5000L,
                  mu = -.05, ar = .3, omega = .01, alpha = .1, beta = .85,
                  skew = 1.5, shape = 5)
    dates <- seq_along(x)
  } else if (dataset == "nasdaq") {
    raw <- read.csv("nasdaq0325.csv", check.names = FALSE)
    price <- rev(raw$Price)
    x <- -100 * log(price[-1] / price[-length(price)])
    dates <- as.Date(rev(raw$Date), format = "%m/%d/%Y")[-1]
    stopifnot(!anyNA(dates), all(diff(dates) > 0), all(is.finite(x)))
    true_mu <- true_sigma <- rep(NA_real_, length(x))
    levels <- list(VaR = .99, EXP = .99855, ES = .975)
    setup <- list(source = "nasdaq0325.csv", first_price = min(as.Date(raw$Date,
                   format = "%m/%d/%Y")), last_price = max(as.Date(raw$Date,
                   format = "%m/%d/%Y")), assessed = length(x) - 1000L)
  } else {
    # Preserve the original concatenation of separately burned-in regimes.
    # The original volatility experiment uses beta=.9 after the break.
    shape_change <- dataset == "structural_shape"
    beta_after <- .9
    set.seed(233)
    model <- list(mu = -.05, ar = .1, omega = .3,
                  alpha = if (shape_change) .1 else .01,
                  beta = if (shape_change) .5 else .1)
    if (shape_change) model <- c(model, list(skew = 1, shape = 6))
    distribution <- if (shape_change) "sstd" else "norm"
    first <- fGarch::garchSim(fGarch::garchSpec(model = model, cond.dist = distribution),
                              n = 3000, n.start = 1000, extended = TRUE)
    if (shape_change) model$shape <- 3 else model$beta <- beta_after
    second <- fGarch::garchSim(fGarch::garchSpec(model = model, cond.dist = distribution),
                               n = 2000, n.start = 1000, extended = TRUE)
    x <- c(as.numeric(first$garch), as.numeric(second$garch))
    true_sigma <- c(as.numeric(first$sigma), as.numeric(second$sigma))
    true_mu <- x - true_sigma * c(as.numeric(first$eps), as.numeric(second$eps))
    dates <- seq_along(x)
    levels <- list(VaR = .9, EXP = .99855, ES = .975)
    setup <- list(seed = 233L, burn_in_each_regime = 1000L, assessed = 4000L,
                  change_at = 2000L, mu = -.05, ar = .1, omega = .3,
                  alpha = if (shape_change) .1 else .01,
                  beta = if (shape_change) .5 else c(.1, beta_after), skew = 1,
                  shape = if (shape_change) c(6, 3) else NA_real_)
  }
  forecast_index <- seq.int(w + 1L, length(x))
  result <- list(dataset = dataset, x = x, dates = dates, window = w,
                 levels = levels, setup = setup, forecast_index = forecast_index,
                 forecast_dates = dates[forecast_index], loss = x[forecast_index],
                 backtest_start = w + 1L, true_mu = true_mu[forecast_index],
                 true_sigma = true_sigma[forecast_index],
                 file_md5 = if (dataset == "nasdaq") unname(tools::md5sum("nasdaq0325.csv")) else NULL)
  save_rds_atomic(result, path)
  result
}

model_risk_values <- function(distribution, levels, parameters = NULL) {
  quantile_levels <- c(levels$VaR, levels$ES)
  if (distribution == "norm") return(list(
    qmodel = qnorm(quantile_levels),
    emodel = expectreg::enorm(asy = levels$EXP),
    esmodel = dnorm(qnorm(levels$ES)) / (1 - levels$ES)))
  shape <- unname(parameters["shape"])
  if (distribution == "std") return(list(
    qmodel = rugarch::qdist("std", p = quantile_levels, mu = 0, sigma = 1, shape = shape),
    emodel = sqrt((shape - 2) / shape) * expectreg::et(asy = levels$EXP, df = shape),
    esmodel = sqrt((shape - 2) / shape) * dt(qt(levels$ES, df = shape), df = shape) /
      (1 - levels$ES) * (shape + qt(levels$ES, df = shape)^2) / (shape - 1)))
  skew <- unname(parameters["skew"])
  qmodel <- rugarch::qdist("sstd", p = quantile_levels, mu = 0, sigma = 1,
                           skew = skew, shape = shape)
  # Integrate the same standardized Fernandez--Steel distribution. This
  # avoids overflow in the original gamma-ratio formula at large fitted df.
  esmodel <- vapply(levels$ES, function(p) {
    q <- rugarch::qdist("sstd", p = p, mu = 0, sigma = 1, skew = skew, shape = shape)
    integrate(function(z) z * rugarch::ddist("sstd", z, mu = 0, sigma = 1,
                 skew = skew, shape = shape), lower = q, upper = Inf,
              rel.tol = 1e-7)$value / (1 - p)
  }, numeric(1))
  # The existing analytic expectile routine uses the unstandardized law.
  emodel <- (est(asy = levels$EXP, shape = shape, skew = skew) -
               mean.st(shape = shape, skew = skew)) / sqrt(var.st(shape = shape, skew = skew))
  list(qmodel = qmodel, emodel = emodel, esmodel = esmodel)
}

prepare_forecasts <- function(dataset, distribution, workers = 4L, cache_dir = "data",
                              chunk_size = 120L, force = FALSE, max_windows = Inf) {
  distribution <- match.arg(distribution, c("norm", "std", "sstd"))
  load_forecast_packages()
  data <- prepare_data(dataset, cache_dir)
  output_path <- file.path(cache_dir, paste0(dataset, "_", distribution, ".rds"))
  checkpoint_path <- sub("[.]rds$", "_checkpoint.rds", output_path)
  if (file.exists(output_path) && !force) return(invisible(readRDS(output_path)))
  count <- length(data$loss)
  results <- vector("list", count)
  if (file.exists(checkpoint_path) && !force) results <- readRDS(checkpoint_path)$results
  spec <- rugarch::ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
                              distribution.model = distribution)
  risk_norm <- if (distribution == "norm") model_risk_values(distribution, data$levels) else NULL
  compute <- function(i) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers({
      # Each fit is seeded separately so parallel scheduling does not affect it.
      set.seed(i)
      fit <- rugarch::ugarchfit(spec, data$x[i:(i + data$window - 1L)], solver = "hybrid")
      if (fit@fit$convergence != 0L) stop("GARCH optimizer did not converge")
      parameters <- coef(fit)
      risks <- if (distribution == "norm") risk_norm else
        model_risk_values(distribution, data$levels, parameters)
      out <- RM.forecasts3(fit, alpha = data$levels$VaR, tau = data$levels$EXP,
                    nu = data$levels$ES, qmodel = risks$qmodel, emodel = risks$emodel,
                    esmodel = risks$esmodel, n = 100000L, seed = i)
      out$parameters <- parameters
      out
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(error = conditionMessage(e)))
    list(value = value, warnings = unique(warnings))
  }
  todo <- which(vapply(results, function(x) is.null(x) || !is.null(x$value$error), logical(1)))
  if (is.finite(max_windows)) todo <- head(todo, max_windows)
  started <- Sys.time()
  chunks <- split(todo, ceiling(seq_along(todo) / chunk_size))
  for (indices in chunks) {
    computed <- if (workers > 1L && .Platform$OS.type != "windows")
      parallel::mclapply(indices, compute, mc.cores = workers, mc.set.seed = FALSE) else
      lapply(indices, compute)
    results[indices] <- computed
    save_rds_atomic(list(results = results, completed = sum(!vapply(results, is.null, logical(1))),
                         dataset = dataset, distribution = distribution), checkpoint_path)
    message(dataset, " / ", distribution, ": ", sum(!vapply(results, is.null, logical(1))),
            "/", count, " windows; elapsed ", round(as.numeric(difftime(Sys.time(), started,
                units = "mins")), 1), " min")
    errors <- indices[vapply(computed, function(x) !is.null(x$value$error), logical(1))]
    if (length(errors)) {
      failure <- paste(vapply(results[errors], function(x) x$value$error, character(1)), collapse = "; ")
      stop("Forecast failures at ", paste(errors, collapse = ","), ": ", failure)
    }
  }
  if (any(vapply(results, is.null, logical(1)))) return(invisible(results))
  fields <- c("VaRmodel", "VaRfhs", "VaRevt", "EXPmodel", "EXPfhs", "EXPevt",
               "ESmodel", "ESfhs", "ESevt")
  output_names <- c("VaR", "VaRfhs", "VaRevt", "EXP", "EXPfhs", "EXPevt",
                    "ES", "ESfhs", "ESevt")
  out <- setNames(lapply(fields, function(field) do.call(rbind,
                     lapply(results, function(x) as.numeric(x$value[[field]])))), output_names)
  for (pair in list(c("xi", "evt.shape"), c("beta", "evt.scale"), c("mut", "mut"), c("sigt", "sigt")))
    out[[pair[1]]] <- vapply(results, function(x) as.numeric(x$value[[pair[2]]]), numeric(1))
  out$par <- do.call(rbind, lapply(results, function(x) x$value$parameters))
  out$metadata <- list(dataset = dataset, distribution = distribution, window = data$window,
                        bootstrap_size = 100000L, forecast_seed = "window index",
                        fit_seed = "window index", warnings = lapply(results, `[[`, "warnings"),
                        session = capture.output(sessionInfo()))
  stopifnot(all(vapply(out[output_names], function(x) all(is.finite(x)), logical(1))))
  save_rds_atomic(out, output_path)
  # Legacy caches retain compatibility with the original helper functions.
  legacy_name <- switch(dataset,
    stationary = switch(distribution, norm = "Sim3norm3.RDATA", std = "Sim3std3.RDATA", sstd = "Sim3sstd3.RDATA"),
    nasdaq = switch(distribution, norm = "NASDAQnorm0325.RDATA", sstd = "NASDAQsstd0325.RDATA", std = "NASDAQstd0325.RDATA"),
    NULL)
  if (!is.null(legacy_name)) save(out, file = file.path(cache_dir, legacy_name))
  invisible(out)
}
