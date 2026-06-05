import numpy as np
from scipy.linalg import block_diag
import argparse
import os
import pandas as pd

def ar1_covariance(block_size, rho):
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov

def covariance_matrix(p, rho, p_blocks):
    if p % p_blocks != 0:
        raise ValueError("p must be divisible by p_blocks.")

    block_size = p // p_blocks
    blocks = []

    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    cov_matrix = block_diag(*blocks)
    return cov_matrix

def read_shrinkage_matrix(p, n_w, df, lam):
    lam_formatted = f"{lam:.4f}"
    filepath = f'/path/to/summary_training/reference_panels_heavy_tail/reference_panel_heavy_tail_p{p}_nw{n_w}_df{df}_lam{lam_formatted}.csv'
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Shrinkage matrix file not found: {filepath}")
    return pd.read_csv(filepath, header=None).values

def main():
    parser = argparse.ArgumentParser(description='Run heavy-tail simulation with specified kappa, noise_std, iteration, p, n, n_w, and df.')
    parser.add_argument('--kappa', type=float, required=True, help='Value of kappa')
    parser.add_argument('--noise_std', type=float, required=True, help='Standard deviation of noise')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration number')
    parser.add_argument('--p', type=int, required=True, help='Number of features')
    parser.add_argument('--n', type=int, required=True, help='Total number of samples in the generated data')
    parser.add_argument('--n_w', type=int, required=True, help='Number of rows in the reference panel')
    parser.add_argument('--df', type=float, required=True, help='Degrees of freedom for multivariate t-distribution')
    args = parser.parse_args()

    kappa = args.kappa
    noise_std = args.noise_std
    iteration = args.iteration
    p = args.p
    n = args.n
    n_w = args.n_w
    df = args.df

    # Initialize random seed using a combination of iteration, kappa, and noise_std
    seed = int((iteration + kappa * 100 + noise_std * 10) * 10000) % (2**(32) - 1)
    np.random.seed(seed)

    # Parameters
    rho = 0.9
    p_blocks = 20
    lambdas = np.linspace(0.01, 10, 100)  # From 0.01 to 10

    covariance = covariance_matrix(p, rho, p_blocks)

    # ---- Load heavy-tail data ----
    base_path = "/path/to/summary_training/Individual_Data_heavy_tail/"
    x_filename = f"X_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    y_filename = f"y_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    beta_filename = f"beta_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"

    x_filepath = os.path.join(base_path, x_filename)
    y_filepath = os.path.join(base_path, y_filename)
    beta_filepath = os.path.join(base_path, beta_filename)

    X = pd.read_csv(x_filepath, header=None).values
    y = pd.read_csv(y_filepath, header=None).values.flatten()
    beta = pd.read_csv(beta_filepath, header=None).values.flatten()

    # Drop last 2500 rows
    X = X[:-2500, :]
    y = y[:-2500]

    # ---- Compute sizes from actual data ----
    n_sum_total = X.shape[0]
    n_ind_total = n_sum_total
    n_sum_train = int(n_sum_total * 4 / 5)
    n_ind_train = n_sum_train
    n_sum_valid = n_sum_total - n_sum_train
    n_ind_valid = n_sum_valid

    print(f"n_sum_total={n_sum_total}, n_sum_train={n_sum_train}, n_sum_valid={n_sum_valid}, n_w={n_w}, df={df}")

    # ---- Save directory ----
    save_dir = "/path/to/summary_training/results/ref_ridge_heavy_tail"
    os.makedirs(save_dir, exist_ok=True)

    # ---- Summary statistics path ----
    X_v = X[n_sum_train:n_sum_total]
    y_v = y[n_sum_train:n_sum_total]

    R_squared_sum_values = []

    X_T_y = X.T @ y
    Cov_X_T_y = np.outer(X_T_y - n_sum_total * covariance @ beta, X_T_y - n_sum_total * covariance @ beta)
    Cov_X_T_y = (Cov_X_T_y + Cov_X_T_y.T) / 2
    mean = (n_sum_train / n_sum_total) * X_T_y

    S_tr = np.random.multivariate_normal(mean, ((n_sum_train * (n_sum_total - n_sum_train)) / n_sum_total**2) * Cov_X_T_y)
    S_v = X_T_y - S_tr

    for lam in lambdas:
        Shrinkage = read_shrinkage_matrix(p, n_w, df, lam)
        ref_ridge_sum = Shrinkage @ S_tr
        dot_product_sum = np.dot(ref_ridge_sum, S_v)
        norm_ref_ridge_sum = np.linalg.norm(X_v @ ref_ridge_sum) * np.linalg.norm(y_v)
        R_squared_sum_values.append((dot_product_sum / norm_ref_ridge_sum)**2)

    # Save R_squared_sum results
    results_sum_df = pd.DataFrame({
        'lambda': lambdas,
        'R_squared_sum': R_squared_sum_values,
        'n_sum_total': [n_sum_total] * len(lambdas),
        'n_sum_train': [n_sum_train] * len(lambdas),
        'n_w': [n_w] * len(lambdas),
        'df': [df] * len(lambdas),
        'p': [p] * len(lambdas),
        'kappa': [kappa] * len(lambdas),
        'noise_std': [noise_std] * len(lambdas)
    })
    results_sum_filename = f'results_sum_heavy_tail_n{n}_p{p}_nw{n_w}_df{df}_kappa{kappa}_noise{noise_std}_iter{iteration}.csv'
    results_sum_filepath = os.path.join(save_dir, results_sum_filename)
    results_sum_df.to_csv(results_sum_filepath, index=False)
    print(f"Saved R_squared_sum results to {results_sum_filepath}")


    indices = np.arange(X.shape[0])
    np.random.shuffle(indices)
    train_indices = indices[:n_ind_train]
    valid_indices = indices[n_ind_train: n_ind_train + n_ind_valid]

    X_tr = X[train_indices]
    y_tr = y[train_indices]
    X_v = X[valid_indices]
    y_v = y[valid_indices]

    X_tr_T_y_tr = X_tr.T @ y_tr
    X_v_T_y_v = X_v.T @ y_v

    R_squared_ind_values = []

    for lam in lambdas:
        Shrinkage = read_shrinkage_matrix(p, n_w, df, lam)
        ref_ridge_ind = Shrinkage @ X_tr_T_y_tr
        dot_product_ind = np.dot(ref_ridge_ind, X_v_T_y_v)
        norm_ref_ridge_ind = np.linalg.norm(X_v @ ref_ridge_ind) * np.linalg.norm(y_v)
        R_squared_ind_values.append((dot_product_ind / norm_ref_ridge_ind)**2)

    # Save R_squared_ind results
    results_ind_df = pd.DataFrame({
        'lambda': lambdas,
        'R_squared_ind': R_squared_ind_values,
        'n_ind_total': [n_ind_total] * len(lambdas),
        'n_ind_train': [n_ind_train] * len(lambdas),
        'n_ind_valid': [n_ind_valid] * len(lambdas),
        'n_w': [n_w] * len(lambdas),
        'df': [df] * len(lambdas),
        'p': [p] * len(lambdas),
        'kappa': [kappa] * len(lambdas),
        'noise_std': [noise_std] * len(lambdas)
    })
    results_ind_filename = f'results_ind_heavy_tail_n{n}_p{p}_nw{n_w}_df{df}_kappa{kappa}_noise{noise_std}_iter{iteration}.csv'
    results_ind_filepath = os.path.join(save_dir, results_ind_filename)
    results_ind_df.to_csv(results_ind_filepath, index=False)
    print(f"Saved R_squared_ind results to {results_ind_filepath}")

if __name__ == '__main__':
    main()
