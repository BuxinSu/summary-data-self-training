# Reproducibility code: Self-training with summary data

This repository contains the numerical code for the paper *"Self-training of summary
data"* (JRSSB submission JRSSB-May-2025-0335). It reproduces every figure in the main
text and the Supplementary Material, covering both the **simulation** studies and the
**UK Biobank** real-data analysis.

The methods compared throughout are:

* `ind` — individual-level (oracle) train/validate tuning (Algorithm "individual-level
  data-based model training");
* `sum` — summary-data self-training via resampled pseudo-training / pseudo-validation
  statistics (Algorithm "summary data-based model training").

Each figure plots the out-of-sample `R^2` (or a covariance-recovery error) of `sum`
against `ind`.

---

## 1. Repository layout

```
.
├── revised_code/                 # CURRENT pipeline used for the revision (start here)
│   ├── 01_*_Generation*.py/.sh   # STEP 1 — data generation (X, y, beta, W, shrinkage)
│   ├── 02_summery_training_*.py  # STEP 2 — summary-data training, R^2 over a lambda grid
│   ├── 03_aggregate_*.py         # STEP 3 — average results over Monte-Carlo iterations
│   ├── 04_test_*.py              # STEP 4 — out-of-sample test R^2 at the selected lambda
│   ├── covariance_recovery.py    # covariance-recovery experiment (compute)
│   ├── plot_*.py                 # plotting scripts (one figure family each)
│   ├── *-DXA.R                   # UK Biobank real-data pipeline (LDpred2 / lassosum2)
│   └── data_generation_description.tex   # formal description of the data-generating models
├── results/                      # pre-computed simulation outputs (small; enables plot-only reruns)
│   ├── covariance_recovery_avg_p5000_rho0.6.csv
│   ├── ref_ridge/                # standard-setting aggregated + test results
│   ├── ref_ridge_heavy_tail/
│   └── ref_ridge_long_range/
├── figures/                      # reference copies of the generated figures
└── <top-level legacy scripts>    # ORIGINAL (V1) scripts used for the first submission;
                                  # kept for provenance and referenced as "closest match"
                                  # where the revised pipeline does not regenerate a V1 figure
```

`revised_code/` is the authoritative, current pipeline. The loose scripts in the
repository root are the **original V1 scripts** from the first submission; they are
retained because several main-text and supplementary figures were first produced by
them, and they are listed as the *closest matching code* in Section 5 below.

---

## 2. Software requirements

* **Python ≥ 3.9** with `numpy`, `pandas`, `matplotlib` (a working LaTeX installation is
  needed for the plotting scripts that set `text.usetex = True`).
* **R ≥ 4.1** with `bigsnpr`, `bigreadr`, `data.table`, `dplyr`, `optparse`, `stringr`,
  `scales` (for the UK Biobank `*-DXA.R` scripts only).

> **Hardcoded paths.** The compute scripts (`02_*`, `04_*`, etc.) default to absolute
> cluster paths such as `/path/to/scratch/...` via `argparse`. Override them with the
> `--data_dir`, `--output_dir`, `--aggregated_csv` flags (see each `*.sh`), or edit the
> defaults. The **plotting** scripts in `revised_code/` already use repository-relative
> paths (`../results`, `../figures`) and run as-is against the bundled `results/`.

---

## 3. Reproduction tracks

### Track A — Simulations (`revised_code/`, four-stage pipeline)

```
STEP 1  01_Individual_Data_Generation*.py     # generate X, y, beta            (DATA GENERATION)
        01_Reference_Panel_Generation*.py     # generate W and (W'W + n_w λ I)^{-1}  (DATA GENERATION)
STEP 2  02_summery_training_ref_ridge*.py     # R^2_sum / R^2_ind over a 100-point λ grid
STEP 3  03_aggregate_ref_ridge_vary_ref*.py   # average over Monte-Carlo iterations
STEP 4  04_test_ref_ridge_vary_ref*.py        # out-of-sample test R^2 at the argmax λ
PLOT    plot_*.py                             # produce the PDF in figures/
```

Each experimental setting has a `_heavy_tail` and `_long_range` variant of steps 1–4
(see Section 4 for the generating models). Because the bundled `results/` already
contains the STEP 3/4 outputs, you can regenerate any figure by running only the
matching `plot_*.py`.

### Track B — UK Biobank real data (`*-DXA.R`)

```
single-ancestry-step1-DXA.R   # ingest + standardize GWAS summary statistics, build LD
single-ancestry-step2-DXA.R   # LDpred2 / lassosum2 fitting + tuning (ind vs sum)
Validation-single-ancestry-DXA.R
val-single-ancestry-analysis-UKBB-varying-tuning-sample-size.R   # vary tuning-set size (100/500/1000)
plot_ldpred2_*.py             # produce the LDpred2 figures
```

> **Real data is not redistributable.** UK Biobank individual-level genotypes/phenotypes
> are governed by a data-use agreement and are **not** included here. The `*-DXA.R`
> scripts are provided for transparency; they require approved UK Biobank access to run.

---

## 4. Data-generation scripts (explicit)

These are the only scripts that *synthesize* data; everything else consumes their
output. The exact generating models are written up in
[`revised_code/data_generation_description.tex`](revised_code/data_generation_description.tex).

| Script | What it generates |
|---|---|
| `revised_code/01_Individual_Data_Generation.py` | `X` (Gaussian, block-AR(1) `Σ`, ρ=0.6, 20 blocks), `β` (sparse Gaussian, var `4/p`), `y = Xβ + ε` |
| `revised_code/01_Reference_Panel_Generation.py` | reference panel `W` (same law as `X`) and the ridge shrinkage matrices `(W'W + n_w λ I)^{-1}` over a 100-point λ grid |
| `revised_code/01_Individual_Data_Generation_heavy_tail.py` | as above but rows of `X` from a multivariate `t_ν` (ν=5), `Cov = Σ` |
| `revised_code/01_Reference_Panel_Generation_heavy_tail.py` | heavy-tailed reference panel + shrinkage matrices |
| `revised_code/01_Individual_Data_Generation_long_range.py` | `Σ = Σ_block + α U Uᵀ` (low-rank cross-block term, rank 5, α=0.1) |
| `revised_code/01_Reference_Panel_Generation_long_range.py` | long-range reference panel + shrinkage matrices |
| `revised_code/covariance_recovery.py` | repeated draws of `W` to measure ‖Σ̂ − Σ‖ vs `n_w` |
| `Individual_Data_Generation.py`, `Reference_Panel_Generation.py`, `External_Reference_Panel_Generation.py` (root) | **V1** generators (identity/`Σ` design) used by the original main-text ridge and marginal simulations |

For the **real-data** track, "data generation" is the GWAS / LD construction in
`single-ancestry-step1-DXA.R`; no synthetic data is involved.

---

## 5. Figure → code map

`Figures_revision/*` and `Figures_V1/*` below are the paths as they appear in the LaTeX
sources. ✅ = the listed code regenerates that exact figure from the bundled `results/`;
↳ = *closest match* (the figure was produced by this code but not regenerated verbatim by
the current pipeline — see notes).

### 5.1 Main text (`main_revision.tex`)

| Figure (file in paper) | Description | Code | Status |
|---|---|---|---|
| `Figures_V1/R_squared_ref_ridge_dense.pdf` + `Figures_V1/testing_R_squared_ref_ridge.pdf` (Fig. `ridge_R2`) | Reference-panel ridge: `R^2` vs `θ`, and test `R^2_sum` vs `R^2_ind` | revised pipeline `revised_code/01–04_*_ref_ridge*.py` + `plot_test_ref_ridge_vary_ref.py`; original: `summery_train.py`, `summery_training_ref_ridge_01_09.py`, `summery_training_ref_ridge_12_22.py` | ↳ closest match (V1 figure; revised pipeline reproduces the same comparison) |
| `Figures_V1/DXA_R_squared_self_training_val_{100,500,1000}.pdf` (Fig. `test_R_alg_comparison`) | UK Biobank DXA: `sum` vs `ind` as the tuning-set size varies (100/500/1000) | `val-single-ancestry-analysis-UKBB-varying-tuning-sample-size.R` (+ `single-ancestry-step1/2-DXA.R`), plotted by `revised_code/plot_ldpred2_*_scatter.py` | ↳ closest match (real data; not redistributable) |

### 5.2 Supplementary Material (`supp_revision.tex`)

| Figure (file in paper) | Description | Code | Status |
|---|---|---|---|
| `Figures_revision/covariance_recovery_operator_norm_p5000_rho0.6.pdf` + `..._frobenius_norm_...pdf` | ‖Σ̂ − Σ‖ (operator / Frobenius) vs `n_w` | `revised_code/covariance_recovery.py` → `revised_code/plot_covariance_recovery.py` | ✅ |
| `Figures_revision/testing_R_squared_ref_ridge_vary_ref_nw_200_num_iter_20.pdf`, `..._nw_20_num_iter_100.pdf`, `..._nw_1000_num_iter_100.pdf` | Ridge `sum` vs `ind` for varying reference-panel size `n_w` | `revised_code/02_summery_training_ref_ridge_vary_ref.py` → `03_aggregate_ref_ridge_vary_ref.py` → `04_test_ref_ridge_vary_ref.py` → `plot_test_ref_ridge_vary_ref.py` | ✅ (run per `n_w`/`num_iter`; output renamed with the suffix) |
| `Figures_revision/test_results_ref_ridge_vary_ref_heavy_tail.pdf` | Heavy-tailed design stress test | `..._heavy_tail.py` chain (`01–04`) → `plot_test_ref_ridge_vary_ref.py` (repointed at `results/ref_ridge_heavy_tail/`) | ✅ |
| `Figures_revision/test_results_ref_ridge_vary_ref_long_range.pdf` | Long-range-dependence stress test | `..._long_range.py` chain (`01–04`) → `plot_test_ref_ridge_vary_ref.py` (repointed at `results/ref_ridge_long_range/`) | ✅ |
| `Figures_revision/testing_R_squared_ref_ridge_mismatch.pdf` | LD-reference mismatch between `W` and `X` | no dedicated script in this folder | ↳ closest match: `02_summery_training_ref_ridge_vary_ref.py` + `plot_test_ref_ridge_vary_ref.py`, run with `W` drawn from a mismatched `Σ` |
| `Figures_V1/R_squared_marginal_cs_sparse.pdf` + `Figures_V1/testing_R_squared_marginal_cs.pdf` | Marginal-thresholding simulation (`sum` vs `ind`) | `summery_training_marginal_01_09.py`, `summery_training_marginal_12_22.py` → `Test_R_Plot_marginal.py` | ↳ closest match (V1 marginal pipeline) |
| `Figures_V1/DXA_R_squared_lassosum2_versus_LDpred2.pdf`, `DXA_R_squared_ensemble_learning_lassosum.pdf`, `DXA_R_squared_ensemble_learning.pdf` | UK Biobank: lassosum2 vs LDpred2, and ensemble learning | `single-ancestry-step2-DXA.R` (+ step1, validation) → `revised_code/plot_ldpred2_*.py` | ↳ closest match (real data) |
| `Figures_revision/LDpred2_ind_vs_sum_top10_sd_paired_boxplot.pdf` | LDpred2 `sum` vs `ind`, paired SD boxplot (top-10) | `revised_code/plot_ldpred2_sd_paired_boxplots.py` (and `plot_ldpred2_sample20_sd_paired_boxplots.py`) | ✅ plotting; ↳ inputs are real-data results |
| `Figures_revision/LDpred2_ind_vs_sum_mean_sample{5,10,20}_scatter.pdf` + `LDpred2_ind_vs_sum_mean_scatter.pdf` | LDpred2 `sum` vs `ind`, mean `R^2` scatter at several tuning-set sizes | `revised_code/plot_ldpred2_sample10_5_mean_scatter.py`, `plot_ldpred2_sample20_mean_scatter.py`, `plot_ldpred2_mean_scatter.py` | ✅ plotting; ↳ inputs are real-data results |

**Auxiliary outputs not used as figures in the paper:** `figures/bar_std_ref_ridge_nw*.pdf`
(`revised_code/plot_bar_std_ref_ridge_vary_ref.py`) and `figures/variance_cost.png` are
diagnostic plots kept for reference.

---

## 6. Figures whose exact code is not in this folder (and the closest match)

* `testing_R_squared_ref_ridge_mismatch.pdf` — no dedicated mismatch driver is included.
  Reproduce by running the `vary_ref` pipeline
  (`02_summery_training_ref_ridge_vary_ref.py` → `03` → `04` → `plot_test_ref_ridge_vary_ref.py`)
  with the reference panel `W` sampled from a covariance that differs from the one used
  for `X`.
* The two **main-text** figures and the V1 **marginal**/**ensemble** figures were
  produced by the original (root-level) scripts. The current `revised_code/` pipeline
  reproduces the same `sum`-vs-`ind` comparisons; the V1 scripts are retained as the
  exact provenance and are the listed closest match above.
* The heavy-tail / long-range / mismatch figures share a single plotting script
  (`plot_test_ref_ridge_vary_ref.py`); switch the input CSV under `results/` and rename
  the output to match the paper's filename.
