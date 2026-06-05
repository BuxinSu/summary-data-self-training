import numpy as np
from scipy.linalg import block_diag
import argparse
import os

def ar1_covariance(block_size, rho):
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov

def generate_heavy_tail_matrix(num_samples, p, rho, p_blocks, df, seed=None):
    """
    Generate samples from a multivariate t-distribution with `df` degrees of freedom,
    zero mean, and covariance matrix equal to the AR(1) block-diagonal Sigma.

    Method:
        Z ~ N(0, Sigma)               (n x p)
        V ~ chi-squared(df)            (n,)
        X = Z * sqrt((df-2) / V)

    The rescaling ensures Cov(X) = Sigma (not df/(df-2)*Sigma).
    Requires df > 2 for the covariance to be finite.
    """
    if p % p_blocks != 0:
        raise ValueError("p must be divisible by p_blocks.")
    if df <= 2:
        raise ValueError("df must be > 2 for finite covariance.")
    if seed is not None:
        np.random.seed(seed)

    block_size = p // p_blocks
    blocks = []
    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    cov_matrix = block_diag(*blocks)
    mean_vector = np.zeros(p)

    # Step 1: Z ~ N(0, Sigma)
    Z = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)

    # Step 2: V ~ chi-squared(df), one per sample
    V = np.random.chisquare(df, size=num_samples)

    # Step 3: X = Z * sqrt((df-2) / V)  so that Cov(X) = Sigma
    X = Z * np.sqrt((df - 2) / V)[:, np.newaxis]

    return cov_matrix, X

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate reference panels from heavy-tailed (multivariate t) W.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--n_w", type=int, required=True, help="Number of rows for the external dataset W")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--df", type=float, default=5.0, help="Degrees of freedom for multivariate t-distribution (must be > 2)")
    parser.add_argument("--output_dir", type=str, default="reference_panels_heavy_tail", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for W generation (depends on p and n_w)")
    parser.add_argument("--chunk_index", type=int, required=True, help="Chunk index (0-9), each chunk processes 10 lambdas")
    args = parser.parse_args()

    n_w = args.n_w
    df = args.df
    p_blocks = 20  # Number of blocks

    # Full lambda grid: 100 values from 0.01 to 10
    all_lambdas = np.linspace(0.01, 10, 100)

    # Select the 10 lambdas for this chunk
    start = args.chunk_index * 10
    end = start + 10
    lambdas = all_lambdas[start:end]
    print(f"Chunk {args.chunk_index}: processing lambdas[{start}:{end}] = [{lambdas[0]:.4f}, ..., {lambdas[-1]:.4f}]")

    os.makedirs(args.output_dir, exist_ok=True)

    # Generate W from multivariate t-distribution (same seed for all chunks with same p, n_w, df)
    I_p = np.identity(args.p)
    covariance, W = generate_heavy_tail_matrix(n_w, args.p, args.rho, p_blocks, df=df, seed=args.seed)
    W_T_W = W.T @ W

    for lam in lambdas:
        regularized_matrix = W_T_W + n_w * lam * I_p
        try:
            reference_panel = np.linalg.inv(regularized_matrix)
        except np.linalg.LinAlgError as e:
            print(f"Matrix inversion failed for p={args.p}, n_w={n_w}, df={df}, lam={lam}: {e}")
            continue

        # Save the result
        file_name = f"reference_panel_heavy_tail_p{args.p}_nw{n_w}_df{df}_lam{lam:.4f}.csv"
        file_path = os.path.join(args.output_dir, file_name)
        np.savetxt(file_path, reference_panel, delimiter=",")
        print(f"Saved reference panel to {file_path}")
