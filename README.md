# Organized code — *Self-training of summary data*

Cleaned, de-duplicated subset of the reproducibility repository. Only the code that
backs a figure in the paper is kept, grouped by experiment. The `results/` CSVs are
**not** included here — only the final `figures/`.

> Filenames keep the original `summery_…` spelling on purpose: the SLURM `.sh`
> wrappers call their `.py` by that exact name, so renaming would break the links.

---

## Layout

```
organized/
├── simulation/                      # synthetic-data studies
│   ├── generators/                  #   STEP 1  base-setting data generators (01_*)
│   │                                #            X, y, β, W, ridge-shrinkage matrices
│   ├── ridge/                       #   STEP 2-4 standard-setting ridge pipeline
│   │   ├── 02_summery_training_*.py #     R²_sum / R²_ind over a 100-point λ grid
│   │   ├── 03_aggregate_*.py        #     average over Monte-Carlo iterations
│   │   └── 04_test_*.py             #     out-of-sample test R² at the argmax λ
│   ├── stress_test/                 #   robustness variants (full 01→04 chains)
│   │   ├── heavy_tail/              #     rows of X from multivariate t_5
│   │   └── long_range/             #     Σ + low-rank cross-block term
│   ├── covariance_recovery/         #   ‖Σ̂ − Σ‖ vs n_w experiment
│   └── data_generation_description.tex
├── realdata_ukbb/                   # UK Biobank DXA real-data pipeline (R)
│   ├── single-ancestry-step1-DXA.R  #   ingest + standardize GWAS sumstats, build LD
│   ├── single-ancestry-step2-DXA.R  #   LDpred2 / lassosum2 fit + tuning (ind vs sum)
│   ├── Validation-single-ancestry-DXA.R
│   └── val-single-ancestry-analysis-UKBB-varying-tuning-sample-size.R
├── plots/                           # plotting scripts that make the paper PDFs
└── figures/
    ├── main_text/                   # the 5 figures in main_revision.tex
    └── supplementary/               # the 18 figures/panels in supp_revision.tex
```

**Run order (simulation, standard setting):** `generators/` (STEP 1) → `ridge/`
02 → 03 → 04 → `plots/`. Each `stress_test/{heavy_tail,long_range}/` folder is a
self-contained 01→04 chain for that design.

---

## Figure → code map

### Main text (`main_revision.tex`)

| Figure (`figures/main_text/…`) | Code in this package |
|---|---|
| `R_squared_ref_ridge_dense.pdf`, `testing_R_squared_ref_ridge.pdf` | `simulation/generators/` → `simulation/ridge/` (02→04) → `plots/plot_test_ref_ridge_vary_ref.py`. *The exact paper PDF was produced by the V1 ridge scripts, which were dropped; this pipeline reproduces the same `sum`-vs-`ind` comparison.* |
| `DXA_R_squared_self_training_val_{100,500,1000}.pdf` | `realdata_ukbb/val-single-ancestry-analysis-UKBB-varying-tuning-sample-size.R` (+ `single-ancestry-step1/step2-DXA.R`) |

### Supplementary (`supp_revision.tex`)

| Figure (`figures/supplementary/…`) | Code in this package |
|---|---|
| `covariance_recovery_operator_norm_…`, `…_frobenius_norm_…` | `simulation/covariance_recovery/covariance_recovery.py` → `plots/plot_covariance_recovery.py` |
| `testing_R_squared_ref_ridge_vary_ref_nw_200_num_iter_20.pdf`, `…_nw_20_num_iter_100.pdf` | `simulation/generators/` → `simulation/ridge/` (`*_vary_ref`, 02→04) → `plots/plot_test_ref_ridge_vary_ref.py` |
| `test_results_ref_ridge_vary_ref_heavy_tail.pdf` | `simulation/stress_test/heavy_tail/` (01→04) → `plots/plot_test_ref_ridge_vary_ref.py` |
| `test_results_ref_ridge_vary_ref_long_range.pdf` | `simulation/stress_test/long_range/` (01→04) → `plots/plot_test_ref_ridge_vary_ref.py` |
| `testing_R_squared_ref_ridge_mismatch.pdf` | `simulation/ridge/` (`*_vary_ref`), run with `W` drawn from a covariance ≠ that of `X` → `plots/plot_test_ref_ridge_vary_ref.py` |
| `R_squared_marginal_cs_sparse.pdf`, `testing_R_squared_marginal_cs.pdf` | **no code in package** — the V1 marginal-thresholding scripts were intentionally dropped; the figure PDFs are kept for completeness |
| `DXA_R_squared_lassosum2_versus_LDpred2.pdf`, `DXA_R_squared_ensemble_learning{,_lassosum}.pdf` | `realdata_ukbb/single-ancestry-step2-DXA.R` (+ step1, Validation) |
| `LDpred2_ind_vs_sum_top10_sd_paired_boxplot.pdf` | `plots/plot_ldpred2_sd_paired_boxplots.py` |
| `LDpred2_ind_vs_sum_mean_sample{5,10}_scatter.pdf` | `plots/plot_ldpred2_sample10_5_mean_scatter.py` |
| `LDpred2_ind_vs_sum_mean_sample20_scatter.pdf` | `plots/plot_ldpred2_sample20_mean_scatter.py` |
| `LDpred2_ind_vs_sum_mean_scatter.pdf` | `plots/plot_ldpred2_mean_scatter.py` |
| `pennprs_oct_tuning_overfit.jpeg` | external screenshot — no code |

---

## Running notes

* **Simulation plots** read a results CSV and write the PDF. The `results/` CSVs were
  excluded, so regenerate them with the `generators → ridge` (or stress_test) chain,
  or `covariance_recovery.py`, before plotting.
* **LDpred2 plots** read real-data result CSVs from `Results/Summery_Validation/`
  (override via `LDPRED2_RESULTS_DIR`). UK Biobank data is **not redistributable**; these
  scripts require approved UKBB access.
* **Hardcoded paths.** The `.sh` submitters and some `.py` defaults point at cluster
  paths like `/path/to/summary_training/…`; edit them (or the `argparse` flags).

## What was dropped (legacy / unused)

* **V1 simulation scripts** (root of the original repo): the marginal experiment
  (`summery_training_marginal_*`, `Test_R_Plot_marginal.py`), the V1 ridge scripts
  (`summery_train.py`, `summery_training_ref_ridge_01_09/12_22.py`), and the V1 data
  generators (`Individual_Data_Generation.py`, `Reference_Panel_Generation.py`,
  `External_Reference_Panel_Generation.py`).
* Scratch experiments: `summery_train_test.py`, `summery_training_test_2.py`,
  `summery_training_test_3.py`.
* `03_aggregate_ref_ridge_vary_ref_all_iter.py` (alternate, no `.sh`, superseded).
* Plot variants producing figures **not** in the paper: `plot_bar_std_ref_ridge_vary_ref.py`,
  `plot_ldpred2_all_sd_top10.py`, `plot_ldpred2_boxplot_23.py`,
  `plot_ldpred2_mean_paired_bars.py`, `plot_ldpred2_sample20_mean_paired_bars.py`,
  `plot_ldpred2_sample20_sd_paired_boxplots.py`, `plot_ldpred2_sd_paired_bars.py`.
* The 1207-line root `single-ancestry-step1-DXA.R` (V1); the cleaned 537-line version is kept.
* The `results/` CSVs and the auxiliary figures `bar_std_ref_ridge_nw*.pdf`, `variance_cost.png`.
