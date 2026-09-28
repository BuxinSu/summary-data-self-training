# Self-training of summary data — code and data

Each folder owns one part of the analysis; every paper figure is produced by one
entry in the table below.

## What each part does

| Folder | Responsible for |
|---|---|
| `simulation/generators/` | **Step 1** — generate the synthetic data: `X`, `y`, `β`, reference panel `W`, and the ridge-shrinkage matrices. |
| `simulation/ridge/` | **Steps 2–4** — ridge self-training (standard design): `02_` train (`R²` over a `λ` grid) → `03_` aggregate over iterations → `04_` out-of-sample test. |
| `simulation/stress_test/heavy_tail/`, `…/long_range/` | The same Step 1→4 chain under a heavy-tailed (`t₅`) design and a long-range-dependence design. |
| `simulation/covariance_recovery/` | Covariance-recovery experiment: `‖Σ̂ − Σ‖` vs reference-panel size `n_w`. |
| `realdata_ukbb/` | UK Biobank DXA pipeline (R): `step1` build LD + standardize GWAS sumstats → `step2` LDpred2 / lassosum2 fit + tuning → validation scripts. |
| [`DXA/`](DXA/) | **Figure 2 raw data** — 100 iteration-level prediction $R^2$ values for each of 71 traits, comparing LDpred2-pseudo with individual-level LDpred2 at tuning sample sizes of 100, 500, and 1,000. |
| `plots/` | Turn the results into the paper PDFs. |
| `figures/` | The figures themselves — `main_text/` and `supplementary/`. |

## Figure 2 raw data

[`DXA/`](DXA/) contains the raw data for all three panels of Figure 2. The analysis uses GWAS summary statistics from 45,622 unrelated White British UK Biobank participants, approximately one million HapMap 3 variants, and a European-ancestry LD reference panel from the 1000 Genomes Project. LDpred2-pseudo uses an 8:2 pseudo-training-to-pseudo-validation ratio. Individual-level LDpred2 uses 100, 500, or 1,000 tuning samples from 2,618 unrelated White non-British participants, with the remaining samples used for testing both methods. Results are averaged across 100 random splits.

See the [DXA README](DXA/README.md) for the experimental setup, file descriptions, trait identifiers, and how to recover the Figure 2 data points.

## Which code makes which figure

### Main text

| Figure | Made by |
|---|---|
| `R_squared_ref_ridge_dense`, `testing_R_squared_ref_ridge` | `simulation/generators/` + `simulation/ridge/` → `plots/plot_test_ref_ridge_vary_ref.py` |
| `DXA_R_squared_self_training_val_{100,500,1000}` | `realdata_ukbb/val-single-ancestry-analysis-UKBB-varying-tuning-sample-size.R` |

### Supplementary

| Figure | Made by |
|---|---|
| `covariance_recovery_{operator,frobenius}_norm` | `simulation/covariance_recovery/covariance_recovery.py` → `plots/plot_covariance_recovery.py` |
| `testing_R_squared_ref_ridge_vary_ref_nw_*` | `simulation/ridge/` (`*_vary_ref`) → `plots/plot_test_ref_ridge_vary_ref.py` |
| `test_results_ref_ridge_vary_ref_heavy_tail` | `simulation/stress_test/heavy_tail/` → `plots/plot_test_ref_ridge_vary_ref.py` |
| `test_results_ref_ridge_vary_ref_long_range` | `simulation/stress_test/long_range/` → `plots/plot_test_ref_ridge_vary_ref.py` |
| `testing_R_squared_ref_ridge_mismatch` | `simulation/ridge/` (`*_vary_ref`, `W` from a mismatched `Σ`) → `plots/plot_test_ref_ridge_vary_ref.py` |
| `DXA_R_squared_lassosum2_versus_LDpred2`, `DXA_R_squared_ensemble_learning{,_lassosum}` | `realdata_ukbb/single-ancestry-step2-DXA.R` |
| `LDpred2_ind_vs_sum_top10_sd_paired_boxplot` | `plots/plot_ldpred2_sd_paired_boxplots.py` |
| `LDpred2_ind_vs_sum_mean_sample{5,10}_scatter` | `plots/plot_ldpred2_sample10_5_mean_scatter.py` |
| `LDpred2_ind_vs_sum_mean_sample20_scatter` | `plots/plot_ldpred2_sample20_mean_scatter.py` |
| `LDpred2_ind_vs_sum_mean_scatter` | `plots/plot_ldpred2_mean_scatter.py` |
| `R_squared_marginal_cs_sparse`, `testing_R_squared_marginal_cs` | *figure kept; marginal V1 code not included* |
| `pennprs_oct_tuning_overfit.jpeg` | external screenshot — no code |

## Notes

- Generate simulation result CSVs with the simulation chain before running the corresponding scripts in `plots/`. Figure 2 raw data are provided in [`DXA/`](DXA/).
- Scripts use placeholder cluster paths (`/path/to/summary_training/…`); edit those or the `argparse`/`optparse` flags for your environment.
