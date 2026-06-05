import numpy as np
from scipy.linalg import block_diag
import argparse
import pandas as pd

def ar1_covariance(block_size, rho):
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov

def generate_gaussian_matrix(num_samples, p, rho, p_blocks):
    if p % p_blocks != 0:
        raise ValueError("p must be divisible by p_blocks.")

    block_size = p // p_blocks
    blocks = []

    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    cov_matrix = block_diag(*blocks)
    mean_vector = np.zeros(p)
    X = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, X


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
    seed = int((iteration + kappa * 100 + noise_std * 10) * 1000) % (2**32 - 1)
    np.random.seed(seed)

    # Parameters
    n_values = [100, 500, 2500, 5000]  # Number of samples
    p_values = [5000, 10000]  # Number of features
    mean_beta = 0  # Mean of the Gaussian distribution for beta
    std_beta = 2  # Standard deviation of the Gaussian distribution for beta
    n_w = 1000  # Number of rows for the external dataset W
    rho = 0.9  # Autocorrelation coefficient
    p_blocks = 20  # Number of blocks
    cutoffs = np.linspace(0.0001, 0.4, 100)

    for p in p_values:
        I_p = np.identity(p)
        for n in n_values:
            n_tr = int(7/10 * n)
            n_v = n - n_tr
            h = float(std_beta)**2 / (float(std_beta)**2 + float(noise_std)**2)
            print(f"Heritability (h): {h}")

            beta = np.zeros(p)
            non_zero_mask = np.random.rand(p) > kappa
            beta[non_zero_mask] = np.random.normal(mean_beta, std_beta / np.sqrt(p), size=np.sum(non_zero_mask))

            R_squared_sum_values = []
            R_squared_ind_values = []

            covariance, X = generate_gaussian_matrix(n, p, rho, p_blocks)
            y = X @ beta + np.random.normal(0, noise_std, size=n)

            indices = np.arange(X.shape[0])
            np.random.shuffle(indices)
            train_indices = indices[:n_tr]
            test_indices = indices[n_tr:]

            X_tr = X[train_indices]
            y_tr = y[train_indices]
            X_v = X[test_indices]
            y_v = y[test_indices]

            X_tr_T_y_tr = X_tr.T @ y_tr
            X_v_T_y_v = X_v.T @ y_v
            X_T_y = X.T @ y

            Cov_X_T_y = np.outer(X_T_y - n * covariance @ beta, X_T_y - n * covariance @ beta)
            Cov_X_T_y = (Cov_X_T_y + Cov_X_T_y.T) / 2
            mean = (n_tr / n) * X_T_y
            
            S_tr = np.random.multivariate_normal(mean, ((n_tr * (n - n_tr)) / n**2) * Cov_X_T_y)
            S_v = X_T_y - S_tr
            
            for idx, cutoff in enumerate(cutoffs):

                S_tr_cutoff = np.where(np.abs(S_tr/n_tr) > cutoff, S_tr, 0)
                dot_product_sum = np.dot(S_tr_cutoff, S_v)
                norm_sum = np.linalg.norm(X_v @ S_tr_cutoff) * np.linalg.norm(y_v)
                R_squared_sum_values.append((dot_product_sum / norm_sum)**2)

                X_tr_T_y_tr_cutoff = np.where(np.abs(X_tr_T_y_tr/n_tr) > cutoff, X_tr_T_y_tr, 0)
                dot_product_ind = np.dot(X_tr_T_y_tr_cutoff, X_v_T_y_v)
                norm_ind = np.linalg.norm(X_v @ X_tr_T_y_tr_cutoff) * np.linalg.norm(y_v)
                R_squared_ind_values.append((dot_product_ind / norm_ind)**2)

            # Save results
            results_df = pd.DataFrame({
                'cutoff': cutoffs,
                'R_squared_sum': R_squared_sum_values,
                'R_squared_ind': R_squared_ind_values
            })
            filename = f'/path/to/summary_training/marginal_sc/results_marginal_test/results_n{n}_p{p}_kappa{kappa}_noise{noise_std}_iter{iteration}.csv'
            results_df.to_csv(filename, index=False)
            print(f"Results saved to {filename}")

if __name__ == '__main__':
    main()
