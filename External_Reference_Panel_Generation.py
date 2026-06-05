import numpy as np
from scipy.linalg import block_diag
import argparse
import os

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
    parser = argparse.ArgumentParser(description="Generate reference panels.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--output_dir", type=str, default="reference_panels", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    n_w = 1000  # Number of rows for the external dataset W
    p_blocks = 20  # Number of blocks

    os.makedirs(args.output_dir, exist_ok=True)

    I_p = np.identity(args.p)
    covariance, W = generate_gaussian_matrix(n_w, args.p, args.rho, p_blocks, seed=args.seed)
    reference_panel = W.T @ W

    # Save the result
    file_name = f"external_reference_panel_p{args.p}.csv"  # Ensure sufficient precision
    file_path = os.path.join(args.output_dir, file_name)
    np.savetxt(file_path, reference_panel, delimiter=",")
    print(f"Saved reference panel to {file_path}")
