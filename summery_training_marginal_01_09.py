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

def main():
    parser = argparse.ArgumentParser(description='Run simulation with specified kappa, noise_std, and iteration.')
    parser.add_argument('--kappa', type=float, required=True, help='Value of kappa')
    parser.add_argument('--noise_std', type=float, required=True, help='Standard deviation of noise')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration number')
    args = parser.parse_args()

    kappa = args.kappa
    noise_std = args.noise_std
    iteration = args.iteration

    # Initialize random seed using a combination of iteration, kappa, and noise_std
    seed = int((iteration + kappa * 100 + noise_std * 10) * 10000) % (2**(32) - 1)
    np.random.seed(seed)

    # Parameters
    n_sum_total = 7500
    n_sum_train = 5000
    n_sum_valid = n_sum_total - n_sum_train
    n_ind_total = 7500
    n_ind_train = 5000
    n_ind_valid_list = [2500]
    p_values = [5000]
    rho = 0.6
    p_blocks = 20
    cutoffs = np.linspace(0.0001, 0.4, 100)

    for p in p_values:
        covariance = covariance_matrix(p, rho, p_blocks)
        base_path = "/path/to/summary_training/Individual_Data/"
        x_train_filename = f"X_train_n{n_sum_total}_p{p}_kappa{kappa}_noise{noise_std}.csv"
        y_train_filename = f"y_train_n{n_sum_total}_p{p}_kappa{kappa}_noise{noise_std}.csv"
        beta_filename = f"beta_p{p}_kappa{kappa}_noise{noise_std}.csv"

        x_filepath = os.path.join(base_path, x_train_filename)
        y_filepath = os.path.join(base_path, y_train_filename)
        beta_filepath = os.path.join(base_path, beta_filename)

        # Load data
        X = pd.read_csv(x_filepath, header=None).values
        y = pd.read_csv(y_filepath, header=None).values.flatten()
        beta = pd.read_csv(beta_filepath, header=None).values.flatten()

        X_v = X[n_sum_train: n_sum_total]
        y_v = y[n_sum_train: n_sum_total]

        R_squared_sum_values = []

        X_T_y = X.T @ y
        Cov_X_T_y = np.outer(X_T_y - n_sum_total * covariance @ beta, X_T_y - n_sum_total * covariance @ beta)
        Cov_X_T_y = (Cov_X_T_y + Cov_X_T_y.T) / 2
        mean = (n_sum_train / n_sum_total) * X_T_y

        S_tr = np.random.multivariate_normal(mean, ((n_sum_train * (n_sum_total - n_sum_train)) / n_sum_total**2) * Cov_X_T_y)
        S_v = X_T_y - S_tr

        for cutoff in cutoffs:
            S_tr_cutoff = np.where(np.abs(S_tr / n_sum_train) > cutoff, S_tr, 0)
            dot_product_sum = np.dot(S_tr_cutoff, S_v)
            norm_sum = np.linalg.norm(X_v @ S_tr_cutoff) * np.linalg.norm(y_v)
            R_squared_sum_values.append((dot_product_sum / norm_sum)**2 if norm_sum != 0 else 0)

        # Save R_squared_sum results
        results_sum_df = pd.DataFrame({
            'cutoff': cutoffs,
            'R_squared_sum': R_squared_sum_values,
            'n_sum_total': [n_sum_total] * len(cutoffs),
            'n_sum_train': [n_sum_train] * len(cutoffs),
            'p': [p] * len(cutoffs),
            'kappa': [kappa] * len(cutoffs),
            'noise_std': [noise_std] * len(cutoffs)
        })
        results_sum_filename = f'results_sum_n{n_sum_valid}_p{p}_kappa{kappa}_noise{noise_std}_iter{iteration}.csv'
        results_sum_filepath = os.path.join("/path/to/summary_training/marginal_sc/results_marginal_test_01_24", results_sum_filename)
        results_sum_df.to_csv(results_sum_filepath, index=False)
        print(f"Saved R_squared_sum results to {results_sum_filepath}")

        for n_ind_valid in n_ind_valid_list:
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

            for cutoff in cutoffs:
                X_tr_T_y_tr_cutoff = np.where(np.abs(X_tr_T_y_tr / n_ind_train) > cutoff, X_tr_T_y_tr, 0)
                dot_product_ind = np.dot(X_tr_T_y_tr_cutoff, X_v_T_y_v)
                norm_ind = np.linalg.norm(X_v @ X_tr_T_y_tr_cutoff) * np.linalg.norm(y_v)
                R_squared_ind_values.append((dot_product_ind / norm_ind)**2 if norm_ind != 0 else 0)

            # Save R_squared_ind results
            results_ind_df = pd.DataFrame({
                'cutoff': cutoffs,
                'R_squared_ind': R_squared_ind_values,
                'n_ind_total': [n_ind_total] * len(cutoffs),
                'n_ind_train': [n_ind_train] * len(cutoffs),
                'n_ind_valid': [n_ind_valid] * len(cutoffs),
                'p': [p] * len(cutoffs),
                'kappa': [kappa] * len(cutoffs),
                'noise_std': [noise_std] * len(cutoffs)
            })
            results_ind_filename = f'results_ind_n{n_ind_valid}_p{p}_kappa{kappa}_noise{noise_std}_iter{iteration}.csv'
            results_ind_filepath = os.path.join("/path/to/summary_training/marginal_sc/results_marginal_test_01_24", results_ind_filename)
            results_ind_df.to_csv(results_ind_filepath, index=False)
            print(f"Saved R_squared_ind results to {results_ind_filepath}")

if __name__ == '__main__':
    main()
