import numpy as np
import pandas as pd
import argparse
import os

def read_shrinkage_matrix(p, n_w, rank, alpha, lam):
    lam_formatted = f"{lam:.4f}"
    filepath = f'/path/to/summary_training/reference_panels_long_range/reference_panel_long_range_p{p}_nw{n_w}_rank{rank}_alpha{alpha}_lam{lam_formatted}.csv'
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Shrinkage matrix file not found: {filepath}")
    return pd.read_csv(filepath, header=None).values

def main():
    parser = argparse.ArgumentParser(description='Test long-range ref_ridge models using best lambdas from aggregated results.')
    parser.add_argument('--aggregated_csv', type=str,
                        default='/path/to/summary_training/results/ref_ridge_long_range/aggregated_ref_ridge_vary_ref_long_range.csv',
                        help='Path to aggregated CSV with averaged R_squared_ind and R_squared_sum')
    parser.add_argument('--data_dir', type=str,
                        default='/path/to/summary_training/Individual_Data_long_range',
                        help='Directory containing long-range X, y, beta CSVs')
    parser.add_argument('--output_dir', type=str,
                        default='/path/to/summary_training/results/ref_ridge_long_range',
                        help='Directory to save the test results CSV')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ---- Step 1: Read aggregated results and find best lambdas ----
    agg_df = pd.read_csv(args.aggregated_csv)
    print(f"Loaded aggregated CSV with {len(agg_df)} rows")

    # Group by (n, n_w, p, rank, alpha, kappa, noise_std) and find lambda that maximizes each R²
    group_cols = ['n', 'n_w', 'p', 'rank', 'alpha', 'kappa', 'noise_std']

    best_lam_ind = agg_df.loc[agg_df.groupby(group_cols)['R_squared_ind'].idxmax()][group_cols + ['lambda']].rename(columns={'lambda': 'lam_ind'})
    best_lam_sum = agg_df.loc[agg_df.groupby(group_cols)['R_squared_sum'].idxmax()][group_cols + ['lambda']].rename(columns={'lambda': 'lam_sum'})

    best_lams = pd.merge(best_lam_ind, best_lam_sum, on=group_cols, how='outer')
    print(f"Found best lambdas for {len(best_lams)} settings")

    # ---- Step 2-4: For each setting, load data, train, and test ----
    results = []

    for _, row in best_lams.iterrows():
        n = int(row['n'])
        n_w = int(row['n_w'])
        p = int(row['p'])
        rank = int(row['rank'])
        alpha = row['alpha']
        kappa = row['kappa']
        noise_std = row['noise_std']
        lam_ind = row['lam_ind']
        lam_sum = row['lam_sum']

        print(f"\nSetting: n={n}, n_w={n_w}, p={p}, rank={rank}, alpha={alpha}, kappa={kappa}, noise_std={noise_std}")
        print(f"  lam_ind={lam_ind:.4f}, lam_sum={lam_sum:.4f}")

        # Load long-range X and y
        x_filepath = os.path.join(args.data_dir, f"X_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv")
        y_filepath = os.path.join(args.data_dir, f"y_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv")

        X = pd.read_csv(x_filepath, header=None).values
        y = pd.read_csv(y_filepath, header=None).values.flatten()

        # Split: train = all but last 2500, test = last 2500
        X_train = X[:-2500, :]
        y_train = y[:-2500]
        X_test = X[-2500:, :]
        y_test = y[-2500:]

        X_train_T_y_train = X_train.T @ y_train
        X_test_T_y_test = X_test.T @ y_test

        # ---- Model trained with lam_sum ----
        Shrinkage_sum = read_shrinkage_matrix(p, n_w, rank, alpha, lam_sum)
        beta_hat_sum = Shrinkage_sum @ X_train_T_y_train

        dot_product_sum = np.dot(beta_hat_sum, X_test_T_y_test)
        norm_sum = np.linalg.norm(X_test @ beta_hat_sum) * np.linalg.norm(y_test)
        test_R_sum = (dot_product_sum / norm_sum) ** 2

        # ---- Model trained with lam_ind ----
        Shrinkage_ind = read_shrinkage_matrix(p, n_w, rank, alpha, lam_ind)
        beta_hat_ind = Shrinkage_ind @ X_train_T_y_train

        dot_product_ind = np.dot(beta_hat_ind, X_test_T_y_test)
        norm_ind = np.linalg.norm(X_test @ beta_hat_ind) * np.linalg.norm(y_test)
        test_R_ind = (dot_product_ind / norm_ind) ** 2

        print(f"  test_R_sum={test_R_sum:.6f}, test_R_ind={test_R_ind:.6f}")

        results.append({
            'n': n,
            'n_w': n_w,
            'p': p,
            'rank': rank,
            'alpha': alpha,
            'kappa': kappa,
            'noise_std': noise_std,
            'lam_sum': lam_sum,
            'lam_ind': lam_ind,
            'test_R_sum': test_R_sum,
            'test_R_ind': test_R_ind
        })

    # ---- Step 5: Save results ----
    results_df = pd.DataFrame(results)
    results_df.sort_values(by=group_cols, inplace=True)
    output_filepath = os.path.join(args.output_dir, 'test_results_ref_ridge_vary_ref_long_range.csv')
    results_df.to_csv(output_filepath, index=False)
    print(f"\nSaved test results ({len(results_df)} rows) to {output_filepath}")

if __name__ == '__main__':
    main()
