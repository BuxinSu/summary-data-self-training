import numpy as np
import os
from scipy.linalg import block_diag
import argparse

def ar1_covariance(block_size, rho):
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov

def long_range_covariance(p, rho, p_blocks, rank, alpha, lr_seed=999):
    """
    Build covariance: Sigma = Sigma_block + alpha * U @ U^T

    - Sigma_block: block-diagonal AR(1) covariance (local/short-range LD)
    - alpha * U @ U^T: low-rank global component (long-range LD across blocks)
    - U is (p x rank), generated with a dedicated fixed seed so the
      long-range structure is identical across all runs with the same (p, rank).

    The result is always positive definite (PD + PSD = PD).
    """
    if p % p_blocks != 0:
        raise ValueError("p must be divisible by p_blocks.")

    # Block-diagonal AR(1)
    block_size = p // p_blocks
    blocks = [ar1_covariance(block_size, rho) for _ in range(p_blocks)]
    sigma_block = block_diag(*blocks)

    # Low-rank global component with a fixed, dedicated seed
    lr_rng = np.random.RandomState(lr_seed)
    U = lr_rng.normal(0, 1.0 / np.sqrt(p), size=(p, rank))
    sigma_lr = alpha * (U @ U.T)

    cov_matrix = sigma_block + sigma_lr
    # Ensure exact symmetry
    cov_matrix = (cov_matrix + cov_matrix.T) / 2

    return cov_matrix

def generate_long_range_matrix(num_samples, p, rho, p_blocks, rank, alpha, lr_seed=999, seed=None):
    """
    Generate samples X ~ N(0, Sigma_long_range).
    """
    if seed is not None:
        np.random.seed(seed)

    cov_matrix = long_range_covariance(p, rho, p_blocks, rank, alpha, lr_seed)
    mean_vector = np.zeros(p)
    X = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, X

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate X, beta, and y with long-range LD covariance.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--n", type=int, default=1000, help="Number of samples")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient for AR(1) blocks")
    parser.add_argument("--rank", type=int, default=5, help="Rank of the low-rank global component")
    parser.add_argument("--alpha", type=float, default=0.1, help="Scaling factor for the low-rank global component")
    parser.add_argument("--kappa", type=float, required=True, help="Value of kappa")
    parser.add_argument("--noise_std", type=float, required=True, help="Noise standard deviation")
    parser.add_argument("--output_dir", type=str, default="Individual_Data_long_range", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--lr_seed", type=int, default=999, help="Fixed seed for generating the low-rank component U (ensures same Sigma across runs)")
    args = parser.parse_args()

    n = args.n
    p = args.p
    rho = args.rho
    rank = args.rank
    alpha = args.alpha
    kappa = args.kappa
    noise_std = args.noise_std
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    # Set global seed ONCE before all random operations
    np.random.seed(args.seed)

    # Generate X from long-range LD covariance
    _, X = generate_long_range_matrix(n, p, rho, p_blocks=20, rank=rank, alpha=alpha, lr_seed=args.lr_seed)

    # Save X
    X_filename = f"X_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, X_filename), X, delimiter=",")
    print(f"Saved X to {os.path.join(output_dir, X_filename)}")

    # Generate beta (continues from same RNG state)
    beta = np.zeros(p)
    non_zero_mask = np.random.rand(p) > kappa
    beta[non_zero_mask] = np.random.normal(0, 2 / np.sqrt(p), size=np.sum(non_zero_mask))

    # Save beta
    beta_filename = f"beta_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, beta_filename), beta, delimiter=",")
    print(f"Saved beta to {os.path.join(output_dir, beta_filename)}")

    # Generate y (noise remains Gaussian)
    y = X @ beta + np.random.normal(0, noise_std, size=n)

    # Save y
    y_filename = f"y_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, y_filename), y, delimiter=",")
    print(f"Saved y to {os.path.join(output_dir, y_filename)}")

    # Compute and save summary stats: X.T @ y  (shape: p,)
    summary_stats = X.T @ y
    summary_filename = f"summary_stats_long_range_n{n}_p{p}_rank{rank}_alpha{alpha}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, summary_filename), summary_stats, delimiter=",")
    print(f"Saved summary stats (X.T @ y) to {os.path.join(output_dir, summary_filename)}")
