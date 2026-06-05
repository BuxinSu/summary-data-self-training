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
    parser.add_argument("--n_w", type=int, required=True, help="Number of rows for the external dataset W")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--output_dir", type=str, default="reference_panels", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for W generation (depends on p and n_w)")
    parser.add_argument("--chunk_index", type=int, required=True, help="Chunk index (0-9), each chunk processes 10 lambdas")
    args = parser.parse_args()

    n_w = args.n_w
    p_blocks = 20  # Number of blocks

    # Full lambda grid: 100 values from 0.01 to 10
    all_lambdas = np.linspace(0.01, 10, 100)

    # Select the 10 lambdas for this chunk
    start = args.chunk_index * 10
    end = start + 10
    lambdas = all_lambdas[start:end]
    print(f"Chunk {args.chunk_index}: processing lambdas[{start}:{end}] = [{lambdas[0]:.4f}, ..., {lambdas[-1]:.4f}]")

    os.makedirs(args.output_dir, exist_ok=True)

    # Generate W using the seed (same seed for all chunks with same p, n_w)
    I_p = np.identity(args.p)
    covariance, W = generate_gaussian_matrix(n_w, args.p, args.rho, p_blocks, seed=args.seed)
    W_T_W = W.T @ W

    for lam in lambdas:
        regularized_matrix = W_T_W + n_w * lam * I_p
        try:
            reference_panel = np.linalg.inv(regularized_matrix)
        except np.linalg.LinAlgError as e:
            print(f"Matrix inversion failed for p={args.p}, n_w={n_w}, lam={lam}: {e}")
            continue

        # Save the result
        file_name = f"reference_panel_p{args.p}_nw{n_w}_lam{lam:.4f}.csv"
        file_path = os.path.join(args.output_dir, file_name)
        np.savetxt(file_path, reference_panel, delimiter=",")
        print(f"Saved reference panel to {file_path}")
