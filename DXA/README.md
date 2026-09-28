# Figure 2 raw data

This folder contains the raw data underlying Figure 2: out-of-sample prediction $R^2$ from LDpred2-pseudo and individual-level LDpred2 across 71 UK Biobank DXA imaging traits. Each result file contains 100 iteration-level $R^2$ values. Figure 2 plots the mean across these 100 iterations for each trait and method.

## Experimental setup

The analysis uses 71 imaging-derived body composition traits from whole-body dual-energy X-ray absorptiometry (DXA), including bone mass, fat-free mass, tissue mass, and lean mass. GWAS summary statistics are derived from 45,622 unrelated White British individuals in the UK Biobank, adjusting for age at imaging, sex, their interactions, and the top 40 genetic principal components. Approximately one million HapMap 3 variants are used for prediction, with an external LD reference panel of European-ancestry individuals from the 1000 Genomes Project.

LDpred2-pseudo selects hyperparameters using resampling-based self-training with an 8:2 pseudo-training-to-pseudo-validation ratio. The pseudo-validation sample size is fixed at 20% of the GWAS sample size.

Individual-level tuning and testing use 2,618 unrelated White non-British individuals from the UK Biobank. In each random split, $N_{\mathrm{tuning}} = n^{\mathrm{(v)}}$ individuals are used for individual-level LDpred2 parameter tuning, and the remaining individuals are used to calculate prediction $R^2$ for both LDpred2 and LDpred2-pseudo. Results are averaged across 100 random splits.

| Figure 2 panel | Individual-level tuning samples | Testing samples |
|---|---:|---:|
| Left | 100 | 2,518 |
| Middle | 500 | 2,118 |
| Right | 1,000 | 1,618 |

## Files

Here, `i` is the trait index, from 1 to 71. `sum` denotes summary-level resampling-based self-training (LDpred2-pseudo), and `ind` denotes individual-level parameter tuning (LDpred2).

| Filename | Method | Figure 2 coordinate |
|---|---|---|
| `LDpred2_sum_i.csv` | LDpred2-pseudo | Horizontal axis in all three panels |
| `LDpred2_ind_i_val_100.csv` | LDpred2, 100 tuning samples | Vertical axis in the left panel |
| `LDpred2_ind_i_val_500.csv` | LDpred2, 500 tuning samples | Vertical axis in the middle panel |
| `LDpred2_ind_i_val_1000.csv` | LDpred2, 1,000 tuning samples | Vertical axis in the right panel |

The 284 result files each contain two columns:

- `iteration`: iteration number, from 1 to 100.
- `r_squared`: out-of-sample prediction $R^2$ for that iteration.

`traits.csv` maps `trait_index` to `ukb_field_id`, following the trait order in Supplementary Table S1. For example, trait index 26 corresponds to UK Biobank field 23244 (android bone mass).

## Recovering the Figure 2 data points

For each trait, take the arithmetic mean of the 100 `r_squared` values in each result file. Use the mean from `LDpred2_sum_i.csv` as the horizontal coordinate and the mean from the corresponding `LDpred2_ind_i_val_*.csv` file as the vertical coordinate. This gives 71 points in each panel, ordered from left to right by individual-level tuning sample sizes of 100, 500, and 1,000. The diagonal line is $y=x$.
