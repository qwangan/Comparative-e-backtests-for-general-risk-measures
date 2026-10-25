# Comparative e-backtests for general risk measures

R code for the numerical study in *Comparative e-backtests for general risk measures* by Z. Jiao, Q. Wang and Y. Zhao. The revised implementation follows Sections 3.1–3.4. `Readme.rtf` retains the historical description.

## Run the experiments

Run from this repository with a recent R installation:

```sh
Rscript Setup.R
Rscript tests_refined.R
Rscript Run_all.R
```

`Setup.R` installs dependencies in `.Rlib`; entry scripts add it to R's library paths. Packages include `rugarch`, `fGarch`, `expectreg`, `ismev`, `skewt`, `sgt`, `ggplot2`, and `patchwork`. Source installations may require R compilation tools.

`Run_all.R` prepares forecasts, runs the experiments, and exports results and figures. An optional argument sets the worker count, for example `Rscript Run_all.R 6`. Rolling model estimation is the expensive step. Caches and checkpoints in `data/` are reused, so interrupted forecasting can resume.

| Script | Experiment |
|---|---|
| `Sim_iid.R` | Traditional ES tests with iid normal losses |
| `Sim_Type1.R` | Comparative VaR type-I error and power |
| `Sim_ts_dataprep.R` / `Sim_Time_Series.R` | Stationary forecasts / comparative tests |
| `Sim_Structural_Change.R` | Both structural-change experiments |
| `Real_dataprep.R` / `Real_testing.R` | NASDAQ forecasts / comparative tests |

To prepare forecasts independently:

```sh
Rscript Run_forecasts.R all 6
Rscript Run_forecasts.R stationary 6
Rscript Run_forecasts.R nasdaq 6
```

Other dataset names are `structural_volatility` and `structural_shape`. Set `FORECAST_MAX_WINDOWS` for a resumable pilot; pilot outputs are not complete results.

## Data setup

The original distributions, forecast models, sample sizes, and attached real-data sample are retained, with the implementation corrections listed below.

| Experiment | Setup |
|---|---|
| Traditional iid | 10 training and 1000 assessment observations; 1000 replications; N(0,1); noisy VaR/ES .95 forecasts; seed `20251007` |
| Comparative iid | 10 training and 1000 assessment observations; 1000 replications; N(0,1); exact, noisy, ±20%, and ±50% VaR .95 forecasts; loss seed `233` |
| Stationary | AR(1)–GARCH(1,1): mean parameters −.05/.3, variance parameters .01/.1/.85; standardized skewed-t shape 5/skewness 1.5; seed `2345`; 1000 burn-in, 500 forecast-training, 500 test-training, and 5000 assessment observations |
| Volatility change | Mean parameters −.05/.1; variance intercept/ARCH coefficient .3/.01; normal innovations; GARCH coefficient .1→.9 at assessment observation 2000; 4000 assessment observations; seed `233` |
| Shape change | Mean parameters −.05/.1; variance parameters .3/.1/.5; skewed-t with fixed skewness 1 and shape 6→3 at assessment observation 2000; 4000 assessment observations; seed `233` |
| NASDAQ | Attached `nasdaq0325.csv`, January 3, 2003–May 30, 2025; negative percentage log returns; 500 observations for estimation and 500 initial forecasts before assessment |

Structural experiments preserve the original concatenation of two separately burned-in regimes. The author confirmed preserving the original code's post-change GARCH coefficient .9.

Stationary models combine normal/t/skewed-t innovations with fully parametric (FP), filtered historical simulation (FHS), and extreme-value (EVT) methods, plus exact forecasts (`opt`): ten models. Real-data and structural experiments use six models, normal/skewed-t × FP/FHS/EVT.

Stationary levels are VaR .90/.95/.99, expectile .96561/.98761/.99855, and VaR–ES .754/.875/.975. All generated levels remain in result files; low/high levels reproduce the original heatmap layout. NASDAQ levels are .99/.99855/.975, respectively. Structural VaR uses .90.

Each fitting window is seeded separately, making random numbers independent of parallel scheduling. Cache metadata include the NASDAQ checksum. Different package versions or optimizers can still change fitted forecasts.

## Revised tests

`refined-tests.R` computes the refined score difference

\[
D_t=S(L_t,R_t)-S(L_t,R_t^*)-h_t^\top I(L_t,R_t^*)
\]

and the explicit bounds in the paper. They hold for every real loss size, without estimating a support bound from the assessment sample.

Each direction uses an equal-weight mixture of fixed betting parameters: `2^seq(-6, 12)` divided by the root mean squared bound width in the initial training observations. The scale is fixed before assessment. Component allocations adapt with accumulated capital; parameters within individual components remain fixed.

Calibration allowances use the standard model's preceding 500 identification errors. The baseline multiplier is **.1**, with sensitivity runs **0/.1/.25/.5**. VaR and VaR–ES use normalized identification components. For expectiles, the allowance in the unnormalized identification function's units is the multiplier times **(1−p) times the standard deviation of the first 500 losses**. This is equivalent to normalizing the identification function by `(1−p) × training SD`. The scale is fixed before assessment and was chosen before examining any time-series comparison results.

Rolling allowances are predictable tolerance rules, not guarantees about the next conditional identification error. Expectile and VaR–ES validity requires the stated conditional calibration assumption for the standard model. VaR uses no correction, so its capital does not depend on this allowance. Sensitivity results do not verify calibration.

- **Overall:** use terminal e-values, two-sided level .05, threshold **40**. Earlier crossings do not fix the final color.
- **Resetting:** report at the first crossing in either direction, then restart both capitals at one. Each report concerns data since its previous restart; forecast and calibration windows continue. For uniform FDR .05, solve `alpha0*(1+log(1/alpha0))=.05` and use threshold `2/alpha0 ≈ 229.7546`.
- **Prespecified structural assessments:** also assess the two periods separated at the experimental change point. These differ from rejection-triggered resets.

Red indicates internal underperformance, green internal outperformance, and yellow inconclusive evidence. Heatmaps place internal models in columns and standards in rows. Reversing roles need not give complementary colors because calibration depends on the standard. Guarantees apply separately to each comparison, not jointly across pairs or sensitivity choices.

Supplementary thresholds 2/5/10 retain comparison with earlier experiments; their two-sided nominal bounds are 1/.4/.2, not .05. Uniform FDR claims use the corrected resetting threshold.

The squared expectile refinement is nonnegative for every loss, so its positive-direction e-value cannot exceed one or produce a green report. Yellow therefore does not establish equal accuracy. Traditional ES tests retain the original predictable GREL rule; their noisy baseline is not an exactly calibrated null and its rejection rate is not a type-I error rate.

## Outputs

The pipeline exports machine-readable results before plotting. This manifest describes outputs, not a claim that a particular run has completed.

| Location | Contents |
|---|---|
| `data/` | Losses, metadata, forecasts, resumable checkpoints |
| `results/traditional_iid_*.csv` | Rejection rates, Monte Carlo standard errors, capital summaries, example path |
| `results/comparative_iid_*.csv` | Directed rejection rates, average-score truth, forecasts |
| `results/<dataset>/overall.csv` | Terminal capitals, colors, score summaries, tolerances, grid |
| `results/<dataset>/resetting_summary.csv` | Red/green report counts by comparison and baseline |
| `results/<dataset>/*_path.csv`, `*_reports.csv` | Main-baseline trajectories and report intervals/dates |
| `results/<dataset>/invalid_forecasts.csv` | Domain violations, when present; forecasts are not silently clipped |
| `results/<dataset>/prespecified_assessments.*` | Two-period structural assessments |
| `results/<dataset>/tests.rds` | Saved trajectories and configuration |
| `results/session-info.txt`, `run-config.rds` | Completed pipeline's R environment and configuration |
| `figures/` | Vector PDF figures and PNG previews |

`Analysis.R` coordinates time-series analyses, `Iid_experiments.R` handles iid experiments, `Plotting.R` supplies the style, and `Publication_figures.R` assembles paper figures. Report intervals here are complete comparative assessments, not block-minimum refinements from the sequential resetting paper.

Reproducibility corrections:

- Replace the missing `Rfns.R` path with `Ecomp-Rfns.R`.
- Seed comparative iid forecast noise with `20251007`; retain its distribution and common sequence across replications.
- Use the stated 1000 held-out observations after ten training observations; the old loop used 1009.
- Correct structural descriptions: GARCH .1→.9, VaR .90, degrees of freedom 6→3 with fixed skewness.
- Distinguish mean log capital from log mean capital.
- Remove orange cells, capital-difference rankings, and unmatched claims of superiority over p-value tests.

Retain run logs and `sessionInfo()` when reproducing results across installations.

## Completed rerun

The full pipeline completed with 44,776 rolling fits and all 3,729 result checks passing. At the main baseline, the stationary experiment gives eight red conclusions, all agreeing with the known average conditional score gaps. Both structural experiments and all NASDAQ comparisons are inconclusive at the stated error levels. The iid test detects 50% VaR under-reporting against the exact benchmark in 90.1% of runs at threshold 40, while smaller differences are often inconclusive.

The diagnostics also show the limits of these results. Estimated simulation benchmarks frequently violate the proposed calibration allowances; baseline zero produces four false stationary conclusions on such paths. NASDAQ's lack of reports persists even under the descriptive supremum over every fixed betting parameter. These results support neither a general power advantage nor unconditional error control for uncalibrated benchmarks.

`Rscript Replot.R` redraws the heatmaps and paper figures from saved results. `Rscript Validate_results.R` audits completed outputs without fitting models again. The pipeline also writes conditional-score truth and calibration and betting diagnostics to `results/`.
