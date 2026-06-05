import numpy as np
from scipy.linalg import block_diag
import argparse
import os

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
    Generate samples W ~ N(0, Sigma_long_range).
    """
    if seed is not None:
        np.random.seed(seed)

    cov_matrix = long_range_covariance(p, rho, p_blocks, rank, alpha, lr_seed)
    mean_vector = np.zeros(p)
    W = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, W

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate reference panels with long-range LD covariance.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--n_w", type=int, required=True, help="Number of rows for the external dataset W")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient for AR(1) blocks")
    parser.add_argument("--rank", type=int, default=5, help="Rank of the low-rank global component")
    parser.add_argument("--alpha", type=float, default=0.1, help="Scaling factor for the low-rank global component")
    parser.add_argument("--output_dir", type=str, default="reference_panels_long_range", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for W generation (depends on p and n_w)")
    parser.add_argument("--lr_seed", type=int, default=999, help="Fixed seed for generating the low-rank component U")
    parser.add_argument("--chunk_index", type=int, required=True, help="Chunk index (0-9), each chunk processes 10 lambdas")
    args = parser.parse_args()

    n_w = args.n_w
    rank = args.rank
    alpha = args.alpha
    p_blocks = 20  # Number of blocks

    # Full lambda grid: 100 values from 0.01 to 10
    all_lambdas = np.linspace(0.01, 10, 100)

    # Select the 10 lambdas for this chunk
    start = args.chunk_index * 10
    end = start + 10
    lambdas = all_lambdas[start:end]
    print(f"Chunk {args.chunk_index}: processing lambdas[{start}:{end}] = [{lambdas[0]:.4f}, ..., {lambdas[-1]:.4f}]")

    os.makedirs(args.output_dir, exist_ok=True)

    # Generate W from long-range LD covariance (same seed for all chunks with same p, n_w)
    I_p = np.identity(args.p)
    covariance, W = generate_long_range_matrix(n_w, args.p, args.rho, p_blocks, rank=rank, alpha=alpha, lr_seed=args.lr_seed, seed=args.seed)
    W_T_W = W.T @ W

    for lam in lambdas:
        regularized_matrix = W_T_W + n_w * lam * I_p
        try:
            reference_panel = np.linalg.inv(regularized_matrix)
        except np.linalg.LinAlgError as e:
            print(f"Matrix inversion failed for p={args.p}, n_w={n_w}, lam={lam}: {e}")
            continue

        # Save the result
        file_name = f"reference_panel_long_range_p{args.p}_nw{n_w}_rank{rank}_alpha{alpha}_lam{lam:.4f}.csv"
        file_path = os.path.join(args.output_dir, file_name)
        np.savetxt(file_path, reference_panel, delimiter=",")
        print(f"Saved reference panel to {file_path}")
