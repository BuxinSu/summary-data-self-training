import os
import re
import glob
import pandas as pd
import argparse

def parse_filename(filepath):
    """
    Extract n, p, nw, df, kappa, noise from heavy-tail filename.
    Example: results_sum_heavy_tail_n12500_p5000_nw1000_df5.0_kappa0.9_noise2.0_iter4.csv
    Returns: {'n': 12500, 'p': 5000, 'n_w': 1000, 'df': 5.0, 'kappa': 0.9, 'noise_std': 2.0}
    """
    basename = os.path.basename(filepath)
    pattern = r'_n(\d+)_p(\d+)_nw(\d+)_df([\d.]+)_kappa([\d.]+)_noise([\d.]+)_iter(\d+)'
    match = re.search(pattern, basename)
    if match is None:
        raise ValueError(f"Could not parse filename: {basename}")
    return {
        'n': int(match.group(1)),
        'p': int(match.group(2)),
        'n_w': int(match.group(3)),
        'df': float(match.group(4)),
        'kappa': float(match.group(5)),
        'noise_std': float(match.group(6)),
    }

def main():
    parser = argparse.ArgumentParser(description='Aggregate heavy-tail ref_ridge vary_ref results over iterations.')
    parser.add_argument('--results_dir', type=str,
                        default='/path/to/summary_training/results/ref_ridge_heavy_tail',
                        help='Directory containing per-iteration result CSVs')
    parser.add_argument('--output_dir', type=str,
                        default='/path/to/summary_training/results/ref_ridge_heavy_tail',
                        help='Directory to save the aggregated CSV')
    args = parser.parse_args()

    results_dir = args.results_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # ---- Read all results_ind_heavy_tail files ----
    ind_pattern = os.path.join(results_dir, 'results_ind_heavy_tail_n*_p*_nw*_df*_kappa*_noise*_iter*.csv')
    ind_files = sorted(glob.glob(ind_pattern))
    print(f"Found {len(ind_files)} results_ind_heavy_tail files")

    if len(ind_files) == 0:
        print("No results_ind_heavy_tail files found. Check the results_dir path.")
        return

    ind_dfs = []
    for f in ind_files:
        meta = parse_filename(f)
        df = pd.read_csv(f)
        # Assign metadata from filename (overrides any columns in the CSV)
        df['n'] = meta['n']
        df['p'] = meta['p']
        df['n_w'] = meta['n_w']
        df['df'] = meta['df']
        df['kappa'] = meta['kappa']
        df['noise_std'] = meta['noise_std']
        ind_dfs.append(df)
    ind_all = pd.concat(ind_dfs, ignore_index=True)

    # Average R_squared_ind over iterations for each (n, n_w, p, df, kappa, noise_std, lambda)
    group_cols = ['n', 'n_w', 'p', 'df', 'kappa', 'noise_std', 'lambda']
    ind_avg = ind_all.groupby(group_cols, as_index=False)['R_squared_ind'].mean()

    # ---- Read all results_sum_heavy_tail files ----
    sum_pattern = os.path.join(results_dir, 'results_sum_heavy_tail_n*_p*_nw*_df*_kappa*_noise*_iter*.csv')
    sum_files = sorted(glob.glob(sum_pattern))
    print(f"Found {len(sum_files)} results_sum_heavy_tail files")

    if len(sum_files) == 0:
        print("No results_sum_heavy_tail files found. Check the results_dir path.")
        return

    sum_dfs = []
    for f in sum_files:
        meta = parse_filename(f)
        df = pd.read_csv(f)
        # Assign metadata from filename (overrides any columns in the CSV)
        df['n'] = meta['n']
        df['p'] = meta['p']
        df['n_w'] = meta['n_w']
        df['df'] = meta['df']
        df['kappa'] = meta['kappa']
        df['noise_std'] = meta['noise_std']
        sum_dfs.append(df)
    sum_all = pd.concat(sum_dfs, ignore_index=True)

    # Average R_squared_sum over iterations for each (n, n_w, p, df, kappa, noise_std, lambda)
    sum_avg = sum_all.groupby(group_cols, as_index=False)['R_squared_sum'].mean()

    # ---- Merge ind and sum results ----
    merged = pd.merge(ind_avg, sum_avg, on=group_cols, how='outer')

    # Reorder columns
    merged = merged[['n', 'n_w', 'p', 'df', 'kappa', 'noise_std', 'lambda', 'R_squared_ind', 'R_squared_sum']]
    merged.sort_values(by=['n', 'n_w', 'p', 'df', 'kappa', 'noise_std', 'lambda'], inplace=True)

    # ---- Save ----
    output_filepath = os.path.join(output_dir, 'aggregated_ref_ridge_vary_ref_heavy_tail.csv')
    merged.to_csv(output_filepath, index=False)
    print(f"Saved aggregated results ({len(merged)} rows) to {output_filepath}")

if __name__ == '__main__':
    main()
