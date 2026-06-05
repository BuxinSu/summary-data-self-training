import numpy as np
import os
from scipy.linalg import block_diag
import argparse

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
        X = Z * sqrt(df / V) * sqrt((df-2)/df)
          = Z * sqrt((df-2) / V)

    The rescaling by sqrt((df-2)/df) ensures Cov(X) = Sigma (not df/(df-2)*Sigma).
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
    # V has shape (n,), broadcast to (n, 1)
    X = Z * np.sqrt((df - 2) / V)[:, np.newaxis]

    return cov_matrix, X

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate X, beta, and y with heavy-tailed (multivariate t) design.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--n", type=int, default=1000, help="Number of samples")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--df", type=float, default=5.0, help="Degrees of freedom for multivariate t-distribution (must be > 2)")
    parser.add_argument("--kappa", type=float, required=True, help="Value of kappa")
    parser.add_argument("--noise_std", type=float, required=True, help="Noise standard deviation")
    parser.add_argument("--output_dir", type=str, default="Individual_Data_heavy_tail", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    n = args.n
    p = args.p
    rho = args.rho
    df = args.df
    kappa = args.kappa
    noise_std = args.noise_std
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    # Set global seed ONCE before all random operations
    np.random.seed(args.seed)

    # Generate X from multivariate t-distribution
    _, X = generate_heavy_tail_matrix(n, p, rho, p_blocks=20, df=df)

    # Save X
    X_filename = f"X_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, X_filename), X, delimiter=",")
    print(f"Saved X to {os.path.join(output_dir, X_filename)}")

    # Generate beta (continues from same RNG state)
    beta = np.zeros(p)
    non_zero_mask = np.random.rand(p) > kappa
    beta[non_zero_mask] = np.random.normal(0, 2 / np.sqrt(p), size=np.sum(non_zero_mask))

    # Save beta
    beta_filename = f"beta_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, beta_filename), beta, delimiter=",")
    print(f"Saved beta to {os.path.join(output_dir, beta_filename)}")

    # Generate y (noise remains Gaussian)
    y = X @ beta + np.random.normal(0, noise_std, size=n)

    # Save y
    y_filename = f"y_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, y_filename), y, delimiter=",")
    print(f"Saved y to {os.path.join(output_dir, y_filename)}")

    # Compute and save summary stats: X.T @ y  (shape: p,)
    summary_stats = X.T @ y
    summary_filename = f"summary_stats_heavy_tail_n{n}_p{p}_df{df}_kappa{kappa}_noise{noise_std}.csv"
    np.savetxt(os.path.join(output_dir, summary_filename), summary_stats, delimiter=",")
    print(f"Saved summary stats (X.T @ y) to {os.path.join(output_dir, summary_filename)}")
