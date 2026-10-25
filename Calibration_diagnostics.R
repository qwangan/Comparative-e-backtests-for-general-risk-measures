# Simulation diagnostics under the known conditional innovation law.
# This checks the tolerance assumption; it does not estimate or certify it
# from the observations. NASDAQ has no known conditional law and is omitted.

if (dir.exists(".Rlib")) .libPaths(c(normalizePath(".Rlib"), .libPaths()))
source("Analysis.R")

fs_t_tail <- function(q, shape, skew) {
  # Survival probability and positive-part expectation of a standardized
  # Fernandez--Steel t random variable, with the fGarch parameterization.
  k <- dt(0, df = shape)
  m <- 2*k*shape/(shape-1)*(skew-1/skew)
  s <- sqrt(shape/(shape-2)*(skew^2+skew^(-2)-1)-m^2)
  y <- m+s*q
  constant <- 2/(skew+1/skew)
  survival <- ifelse(y < 0,
    1-constant/skew*pt(skew*y, df = shape),
    constant*skew*pt(y/skew, df = shape, lower.tail = FALSE))
  tail_first <- ifelse(y < 0,
    m+constant/skew^2*(shape+(skew*y)^2)/(shape-1)*dt(skew*y, df = shape),
    constant*skew^2*(shape+(y/skew)^2)/(shape-1)*dt(y/skew, df = shape))
  list(survival = survival, positive_part = pmax(0, (tail_first-y*survival)/s))
}

true_loss_tail <- function(data, threshold) {
  q <- (threshold-data$true_mu)/data$true_sigma
  if (data$dataset == "structural_volatility") {
    survival <- pnorm(q, lower.tail = FALSE)
    hinge <- data$true_sigma*pmax(0, dnorm(q)-q*survival)
  } else {
    shape <- if (data$dataset == "stationary") 5 else
      ifelse(data$forecast_index <= 3000L, 6, 3)
    skew <- if (data$dataset == "stationary") 1.5 else 1
    tail <- fs_t_tail(q, shape, skew)
    survival <- tail$survival
    hinge <- data$true_sigma*tail$positive_part
  }
  list(survival = survival, hinge = hinge)
}

conditional_identification <- function(data, standard, p, risk) {
  if (risk == "VaRES") {
    standard <- as.matrix(standard)
    tail <- true_loss_tail(data, standard[,2])
    return(cbind(VaR = tail$survival/(1-p)-1,
      ES = tail$hinge/((1-p)*(standard[,1]-standard[,2]))-1))
  }
  tail <- true_loss_tail(data, as.numeric(standard))
  if (risk == "VaR") return(cbind(VaR = tail$survival/(1-p)-1))
  cbind(EXP = (2*p-1)*tail$hinge+(1-p)*(data$true_mu-as.numeric(standard)))
}

check_tail_formula <- function() {
  load_forecast_packages()
  error <- numeric()
  for (shape in c(3, 5, 6, 30)) for (skew in c(.6, 1, 1.5))
    for (q in c(-3, -.1, 0, 1, 3)) {
      calculated <- fs_t_tail(q, shape, skew)
      survival <- 1-fGarch::psstd(q, nu = shape, xi = skew)
      hinge <- integrate(function(z) (z-q)*fGarch::dsstd(z, nu = shape, xi = skew),
                          lower = q, upper = Inf, rel.tol = 1e-8)$value
      error <- c(error, abs(calculated$survival-survival),
                 abs(calculated$positive_part-hinge))
    }
  stopifnot(max(error) < 1e-7)
  max(error)
}

run_calibration_diagnostics <- function(datasets = c("stationary", "structural_volatility",
                                                     "structural_shape"), output = "results") {
  max_formula_error <- check_tail_formula()
  rows <- list()
  for (dataset in datasets) {
    bundle <- collect_forecasts(dataset)
    data <- bundle$data
    indices <- seq.int(data$backtest_start, length(data$loss))
    periods <- list(full = indices)
    if (startsWith(dataset, "structural_")) {
      periods$before <- indices[data$forecast_index[indices] <= 3000L]
      periods$after <- indices[data$forecast_index[indices] > 3000L]
    }
    for (model in names(bundle$models)) for (risk in c("VaR", "EXP", "VaRES")) {
      levels <- if (risk == "VaRES") data$levels$ES else data$levels[[risk]]
      for (level_index in seq_along(levels)) {
        p <- levels[level_index]
        standard <- if (risk == "VaRES") cbind(ES = bundle$models[[model]]$ES[,level_index],
          VaR = bundle$models[[model]]$VaR[,length(data$levels$VaR)+level_index]) else
          bundle$models[[model]][[risk]][,level_index]
        parts <- tryCatch(refined_score(data$loss, standard, standard, p, risk),
                           error = function(e) e)
        if (inherits(parts, "error")) stop(dataset, "/", model, "/", risk, ": ", conditionMessage(parts))
        expected <- conditional_identification(data, standard, p, risk)
        for (epsilon in c(0, .1, .25, .5)) {
          baseline <- baseline_tolerance(data, risk, epsilon, p = p)
          tolerance <- refined_tolerances(parts$I, w = 500L, epsilon = baseline)
          numerical_tolerance <- 1e-8
          below <- expected < -tolerance$minus-numerical_tolerance
          above <- expected > tolerance$plus+numerical_tolerance
          outside <- below | above
          required <- pmax(expected-pmax(tolerance$bias,0),
                            -expected-pmax(-tolerance$bias,0),0)
          required <- sweep(required, 2L, baseline_tolerance(data, risk, 1, p = p), "/")
          for (period in names(periods)) {
            ids <- periods[[period]]
            for (component in c(colnames(expected), "joint")) {
              which_components <- if (component == "joint") seq_len(ncol(expected)) else
                match(component, colnames(expected))
              probability <- function(a) mean(rowSums(a[ids,which_components,drop = FALSE]) > 0)
              rows[[length(rows)+1L]] <- data.frame(dataset, standard = model, risk,
                level = p, component, period, epsilon,
                baseline = paste(signif(baseline, 8), collapse = ";"), observations = length(ids),
                violation_fraction = probability(outside),
                lower_violation_fraction = probability(below),
                upper_violation_fraction = probability(above),
                mean_absolute_conditional_error = mean(abs(expected[ids,which_components,drop = FALSE])),
                max_absolute_conditional_error = max(abs(expected[ids,which_components,drop = FALSE])),
                minimum_common_baseline_for_realized_path = max(required[ids,which_components,drop = FALSE]),
                numerical_tolerance,
                refinement_uses_component = risk == "EXP" || (risk == "VaRES" && component != "VaR"),
                conditioning = "true DGP states and past observations")
            }
          }
        }
      }
    }
    message("Calibration diagnostics: ", dataset, " complete")
  }
  result <- do.call(rbind, rows)
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  write.csv(result, file.path(output, "calibration_diagnostics.csv"), row.names = FALSE)
  summary <- aggregate(violation_fraction ~ dataset+risk+epsilon,
    subset(result, period == "full" & component == "joint" & standard != "opt"), range)
  summary <- data.frame(summary[c("dataset", "risk", "epsilon")],
    minimum_violation_fraction = summary$violation_fraction[,1],
    maximum_violation_fraction = summary$violation_fraction[,2])
  write.csv(summary, file.path(output, "calibration_summary.csv"), row.names = FALSE)
  write.csv(data.frame(check = "FS-t tail formulas versus numerical integration",
                       max_absolute_error = max_formula_error),
              file.path(output, "calibration_formula_check.csv"), row.names = FALSE)
  invisible(result)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  run_calibration_diagnostics(if (length(arguments)) arguments else
    c("stationary", "structural_volatility", "structural_shape"))
}
