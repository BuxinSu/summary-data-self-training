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
    blocks = [ar1_covariance(block_size, rho) for _ in range(p_blocks)]
    return block_diag(*blocks)

def generate_gaussian_matrix(num_samples, p, rho, p_blocks, seed=None):
    if seed is not None:
        np.random.seed(seed)
    cov_matrix = covariance_matrix(p, rho, p_blocks)
    mean_vector = np.zeros(p)
    W = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, W

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Measure how well W^T W / n_w recovers Sigma.")
    parser.add_argument("--p", type=int, required=True, help="Number of features")
    parser.add_argument("--rho", type=float, default=0.6, help="Autocorrelation coefficient")
    parser.add_argument("--n_reps", type=int, default=50, help="Number of repetitions per n_w")
    parser.add_argument("--output_dir", type=str, default="results/covariance_recovery", help="Directory to save results")
    parser.add_argument("--seed", type=int, default=42, help="Base random seed")
    args = parser.parse_args()

    p = args.p
    rho = args.rho
    p_blocks = 20
    n_reps = args.n_reps

    os.makedirs(args.output_dir, exist_ok=True)

    # True covariance
    Sigma = covariance_matrix(p, rho, p_blocks)

    # Range of n_w values to test
    n_w_values = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]

    results = []

    for n_w in n_w_values:
        for rep in range(n_reps):
            seed = args.seed * 1000 + n_w + rep
            np.random.seed(seed)

            # Generate W
            mean_vector = np.zeros(p)
            W = np.random.multivariate_normal(mean_vector, Sigma, size=n_w)

            # Sample covariance
            sample_cov = W.T @ W / n_w

            # Difference
            diff = sample_cov - Sigma

            # Operator norm (largest singular value = spectral norm)
            operator_norm = np.linalg.norm(diff, ord=2)

            # Frobenius norm
            frobenius_norm = np.linalg.norm(diff, ord='fro')

            results.append({
                'p': p,
                'rho': rho,
                'n_w': n_w,
                'rep': rep,
                'operator_norm': operator_norm,
                'frobenius_norm': frobenius_norm
            })

            print(f"p={p}, n_w={n_w}, rep={rep}: operator_norm={operator_norm:.6f}, frobenius_norm={frobenius_norm:.6f}")

    df = pd.DataFrame(results)

    # Save raw results (all reps)
    raw_filepath = os.path.join(args.output_dir, f"covariance_recovery_raw_p{p}_rho{rho}.csv")
    df.to_csv(raw_filepath, index=False)
    print(f"\nSaved raw results to {raw_filepath}")

    # Save averaged results
    df_avg = df.groupby(['p', 'rho', 'n_w']).agg(
        operator_norm_mean=('operator_norm', 'mean'),
        operator_norm_std=('operator_norm', 'std'),
        frobenius_norm_mean=('frobenius_norm', 'mean'),
        frobenius_norm_std=('frobenius_norm', 'std'),
        n_reps=('rep', 'count')
    ).reset_index()

    avg_filepath = os.path.join(args.output_dir, f"covariance_recovery_avg_p{p}_rho{rho}.csv")
    df_avg.to_csv(avg_filepath, index=False)
    print(f"Saved averaged results to {avg_filepath}")
