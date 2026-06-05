import numpy as np
import os
from scipy.linalg import block_diag
import argparse

def ar1_covariance(block_size, rho):
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov

def generate_gaussian_matrix(num_samples, p, rho, p_blocks, seed=None):
    if p % p_blocks != 0:
        raise ValueError("p must be divisible by p_blocks.")
    if seed is not None:
        np.random.seed(seed)

    block_size = p // p_blocks
    blocks = []

    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    cov_matrix = block_diag(*blocks)
    mean_vector = np.zeros(p)
    X = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, X

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate X, y, and compute X^Ty.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--n", type=int, default=10000, help="Number of samples")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--kappa", type=float, required=True, help="Value of kappa")
    parser.add_argument("--noise_std", type=float, required=True, help="Noise standard deviation")
    parser.add_argument("--output_dir", type=str, default="results_marginal_test", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    n = args.n
    p = args.p
    rho = args.rho
    kappa = args.kappa
    noise_std = args.noise_std
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    # Generate X
    p_blocks = 20
    _, X = generate_gaussian_matrix(n, p, rho, p_blocks, seed=args.seed)

    # Split X into training and testing sets
    n_train = 7500
    n_test = n - n_train
    X_train = X[:n_train]
    X_test = X[n_train:]

    # Save X_train and X_test
    X_train_filename = f"X_train_n{n_train}_p{p}_kappa{kappa}_noise{noise_std}.csv"
    X_test_filename = f"X_test_n{n_test}_p{p}_kappa{kappa}_noise{noise_std}.csv"
    X_train_filepath = os.path.join(output_dir, X_train_filename)
    X_test_filepath = os.path.join(output_dir, X_test_filename)
    np.savetxt(X_train_filepath, X_train, delimiter=",")
    np.savetxt(X_test_filepath, X_test, delimiter=",")
    print(f"Saved X_train to {X_train_filepath}")
    print(f"Saved X_test to {X_test_filepath}")

    # Generate beta
    beta = np.zeros(p)
    non_zero_mask = np.random.rand(p) > kappa
    beta[non_zero_mask] = np.random.normal(0, 2 / np.sqrt(p), size=np.sum(non_zero_mask))
    # Save beta
    beta_filename = f"beta_p{p}_kappa{kappa}_noise{noise_std}.csv"
    beta_filepath = os.path.join(output_dir, beta_filename)
    np.savetxt(beta_filepath, beta, delimiter=",")
    print(f"Saved beta to {beta_filepath}")

    # Generate y
    y = X @ beta + np.random.normal(0, noise_std, size=n)

    # Split y into training and testing sets
    y_train = y[:n_train]
    y_test = y[n_train:]

    # Save y_train and y_test
    y_train_filename = f"y_train_n{n_train}_p{p}_kappa{kappa}_noise{noise_std}.csv"
    y_test_filename = f"y_test_n{n_test}_p{p}_kappa{kappa}_noise{noise_std}.csv"
    y_train_filepath = os.path.join(output_dir, y_train_filename)
    y_test_filepath = os.path.join(output_dir, y_test_filename)
    np.savetxt(y_train_filepath, y_train, delimiter=",")
    np.savetxt(y_test_filepath, y_test, delimiter=",")
    print(f"Saved y_train to {y_train_filepath}")
    print(f"Saved y_test to {y_test_filepath}")
