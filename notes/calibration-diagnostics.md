The simulation diagnostic uses the known normal or standardized Fernandez--Steel
t innovation law to calculate conditional expectations of the identification
functions. The same preceding 500 observations and baseline tolerances as in
the tests are used. Expectile baselines are multiplied by `(1-p)` times the
standard deviation of the first 500 verification losses, which precede the
assessed observations. This accounts for the scale of the canonical
identification function at extreme expectile levels.

`calibration_diagnostics.csv` reports the fraction of assessed observations at
which each conditional mean falls outside its adaptive interval, for each
standard model, risk measure, risk level and baseline. Joint rows check all
components simultaneously. A numerical tolerance of `1e-8` avoids counting
root-finding roundoff as a calibration failure for optimal forecasts. The
minimum common baseline column is calculated for the realized path, and does
not establish a bound for future observations. These are diagnostics under the known simulation
law; historical calibration errors do not certify the conditional assumption.

The conditioning includes the initial states of the simulated regimes. The
original structural experiments concatenate two separately burned-in series;
their new initial state at the change point is not observed from the preceding
losses. The state-conditional diagnostic is therefore stronger than conditioning
on the observed losses alone. A violation here does not by itself prove failure
of the restriction under the smaller filtration.

The VaR score uses no identification adjustment. Its e-process validity does
not need the reported VaR calibration restriction. For VaR--ES, the chosen
adjustment uses only the ES identification component; the `refinement_uses_component`
column records this distinction. No corresponding diagnostic is possible for
NASDAQ without an additional conditional-distribution assumption.
