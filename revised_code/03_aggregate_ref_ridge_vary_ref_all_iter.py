import os
import re
import glob
import pandas as pd
import argparse

def parse_filename(filepath):
    """
    Extract n, p, nw, kappa, noise, iter from filename.
    Example: results_sum_n12500_p5000_nw1000_kappa0.9_noise2.0_iter4.csv
    Returns: {'n': 12500, 'p': 5000, 'n_w': 1000, 'kappa': 0.9, 'noise_std': 2.0, 'iter': 4}
    """
    basename = os.path.basename(filepath)
    pattern = r'_n(\d+)_p(\d+)_nw(\d+)_kappa([\d.]+)_noise([\d.]+)_iter(\d+)'
    match = re.search(pattern, basename)
    if match is None:
        raise ValueError(f"Could not parse filename: {basename}")
    return {
        'n': int(match.group(1)),
        'p': int(match.group(2)),
        'n_w': int(match.group(3)),
        'kappa': float(match.group(4)),
        'noise_std': float(match.group(5)),
        'iter': int(match.group(6)),
    }

def main():
    parser = argparse.ArgumentParser(description='Aggregate ref_ridge vary_ref results: keep all iterations, filter lambda=10.0.')
    parser.add_argument('--results_dir', type=str,
                        default='/path/to/summary_training/results/ref_ridge_v1',
                        help='Directory containing per-iteration result CSVs')
    parser.add_argument('--output_dir', type=str,
                        default='/path/to/summary_training/results/ref_ridge_v1',
                        help='Directory to save the aggregated CSV')
    parser.add_argument('--lam', type=float, default=10.0,
                        help='Lambda value to filter on (default: 10.0)')
    args = parser.parse_args()

    results_dir = args.results_dir
    output_dir = args.output_dir
    lam_filter = args.lam
    os.makedirs(output_dir, exist_ok=True)

    # ---- Read all results_ind files ----
    ind_pattern = os.path.join(results_dir, 'results_ind_n*_p*_nw*_kappa*_noise*_iter*.csv')
    ind_files = sorted(glob.glob(ind_pattern))
    print(f"Found {len(ind_files)} results_ind files")

    if len(ind_files) == 0:
        print("No results_ind files found. Check the results_dir path.")
        return

    ind_dfs = []
    for f in ind_files:
        meta = parse_filename(f)
        df = pd.read_csv(f)
        # Filter to lambda = lam_filter
        df = df[df['lambda'] == lam_filter].copy()
        if len(df) == 0:
            continue
        # Assign metadata from filename
        df['n'] = meta['n']
        df['p'] = meta['p']
        df['n_w'] = meta['n_w']
        df['kappa'] = meta['kappa']
        df['noise_std'] = meta['noise_std']
        df['iter'] = meta['iter']
        ind_dfs.append(df)
    ind_all = pd.concat(ind_dfs, ignore_index=True)
    print(f"Collected {len(ind_all)} rows from results_ind with lambda={lam_filter}")

    # ---- Read all results_sum files ----
    sum_pattern = os.path.join(results_dir, 'results_sum_n*_p*_nw*_kappa*_noise*_iter*.csv')
    sum_files = sorted(glob.glob(sum_pattern))
    print(f"Found {len(sum_files)} results_sum files")

    if len(sum_files) == 0:
        print("No results_sum files found. Check the results_dir path.")
        return

    sum_dfs = []
    for f in sum_files:
        meta = parse_filename(f)
        df = pd.read_csv(f)
        # Filter to lambda = lam_filter
        df = df[df['lambda'] == lam_filter].copy()
        if len(df) == 0:
            continue
        # Assign metadata from filename
        df['n'] = meta['n']
        df['p'] = meta['p']
        df['n_w'] = meta['n_w']
        df['kappa'] = meta['kappa']
        df['noise_std'] = meta['noise_std']
        df['iter'] = meta['iter']
        sum_dfs.append(df)
    sum_all = pd.concat(sum_dfs, ignore_index=True)
    print(f"Collected {len(sum_all)} rows from results_sum with lambda={lam_filter}")

    # ---- Merge ind and sum results ----
    merge_cols = ['n', 'n_w', 'p', 'kappa', 'noise_std', 'lambda', 'iter']
    merged = pd.merge(ind_all[merge_cols + ['R_squared_ind']],
                       sum_all[merge_cols + ['R_squared_sum']],
                       on=merge_cols, how='outer')

    # Reorder columns
    merged = merged[['n', 'n_w', 'p', 'kappa', 'noise_std', 'lambda', 'iter', 'R_squared_ind', 'R_squared_sum']]
    merged.sort_values(by=['n', 'n_w', 'p', 'kappa', 'noise_std', 'iter'], inplace=True)

    # ---- Save ----
    output_filepath = os.path.join(output_dir, 'aggregated_ref_ridge_vary_ref_all_iter.csv')
    merged.to_csv(output_filepath, index=False)
    print(f"Saved aggregated results ({len(merged)} rows) to {output_filepath}")

if __name__ == '__main__':
    main()
