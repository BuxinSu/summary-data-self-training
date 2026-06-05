import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import glob
from scipy.linalg import block_diag


!apt-get install -y texlive texlive-latex-extra texlive-fonts-recommended dvipng cm-super

# Set the text properties to use LaTeX fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30  # 30 for two pictures in a row
plt.rcParams['figure.figsize'] = (13, 8)  # width, height in inches
plt.rcParams['figure.dpi'] = 300  # dots per inch
plt.rcParams['font.family'] = 'serif'  # Use serif fonts to match LaTeX's default


def ar1_covariance(block_size, rho):
    """
    Generate an AR(1) covariance matrix of size block_size x block_size.

    Parameters:
    - block_size (int): Size of the AR(1) block.
    - rho (float): Autocorrelation coefficient (-1 < rho < 1).

    Returns:
    - cov (ndarray): AR(1) covariance matrix.
    """
    indices = np.arange(block_size)
    cov = rho ** np.abs(np.subtract.outer(indices, indices))
    return cov



def generate_gaussian_vector(p, rho=0.5, p_blocks=10):
    """
    Generate an n-dimensional Gaussian vector with mean zero and a blockwise AR(1) covariance matrix.

    Parameters:
    - n (int): Total dimension of the Gaussian vector.
    - rho (float): Autocorrelation coefficient for AR(1) blocks.
    - n_blocks (int): Number of blocks along the diagonal.

    Returns:
    - x (ndarray): Generated Gaussian vector of size n.
    """
    # Ensure that n is divisible by n_blocks
    if p % p_blocks != 0:
        raise ValueError("n must be divisible by n_blocks.")

    block_size = p // p_blocks
    blocks = []

    # Generate AR(1) covariance matrices for each block
    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    # Create the block diagonal covariance matrix
    cov_matrix = block_diag(*blocks)

    # Mean vector of zeros
    mean_vector = np.zeros(p)

    # Generate the Gaussian vector
    x = np.random.multivariate_normal(mean_vector, cov_matrix)
    return x



def generate_gaussian_matrix(num_samples, p, rho, p_blocks):
    """
    Generate a matrix where each row is an n-dimensional Gaussian vector with mean zero
    and a blockwise AR(1) covariance matrix.

    Parameters:
    - num_samples (int): Number of samples (rows) in the matrix.
    - n (int): Total dimension of each Gaussian vector (number of columns).
    - rho (float): Autocorrelation coefficient for AR(1) blocks.
    - n_blocks (int): Number of blocks along the diagonal.

    Returns:
    - X (ndarray): Generated Gaussian matrix of size (num_samples, n).
    """
    # Ensure that n is divisible by n_blocks
    if p % p_blocks != 0:
        raise ValueError("n must be divisible by n_blocks.")

    block_size = p // p_blocks
    blocks = []

    # Generate AR(1) covariance matrices for each block
    for _ in range(p_blocks):
        cov_block = ar1_covariance(block_size, rho)
        blocks.append(cov_block)

    # Create the block diagonal covariance matrix
    cov_matrix = block_diag(*blocks)

    # Mean vector of zeros
    mean_vector = np.zeros(p)

    # Generate the Gaussian matrix
    X = np.random.multivariate_normal(mean_vector, cov_matrix, size=num_samples)
    return cov_matrix, X


# Load the results from the previous code
results_df = pd.read_csv('/path/to/summary_training/max_cutoff_marginal.csv')

results_df = results_df[results_df['n'] != 5000]

# Initialize parameters
mean_beta = 0
std_beta = 2
n_t = 4000
rho = 0.9
p_blocks = 20

# Arrays to store test R^2 values
R_squared_sum_test_values = []
R_squared_ind_test_values = []

# Loop over each combination of n, kappa, noise_std, and cutoff from filtered results
for _, row in results_df.iterrows():
    n = int(row['n'])
    p = int(row['p'])
    kappa = row['kappa']
    noise_std = row['noise_std']
    max_cutoff_sum = row['cutoff_max_R_squared_sum']
    max_cutoff_ind = row['cutoff_max_R_squared_ind']

    # Define beta as per the specified model
    beta = np.zeros(p)
    non_zero_mask = np.random.rand(p) > kappa
    beta[non_zero_mask] = np.random.normal(mean_beta, std_beta / np.sqrt(p), size=np.sum(non_zero_mask))

    # Generate the training data X_tr, y_tr
    covariance, X_tr = generate_gaussian_matrix(n, p, rho, p_blocks)
    y_tr = X_tr @ beta + np.random.normal(0, noise_std, size=n)

    # Compute S_tr for the sum method
    S_tr = X_tr.T @ y_tr

    # Generate the test data T, y_t
    _, T = generate_gaussian_matrix(n_t, p, rho, p_blocks)
    y_t = T @ beta + np.random.normal(0, noise_std, size=n_t)

    # Compute S_test = T^T y_t
    S_test = T.T @ y_t

    # Apply the cutoff for the sum method
    S_tr_cutoff = np.where(np.abs(S_tr/n) > max_cutoff_sum, S_tr, 0)
    dot_product_sum = np.dot(S_tr_cutoff, S_test)
    norm_sum = np.linalg.norm(T @ S_tr_cutoff) * np.linalg.norm(y_t)
    R_squared_sum_test_values.append((dot_product_sum / norm_sum)**2)

    # Apply the cutoff for the individual method
    X_tr_T_y_tr_cutoff = np.where(np.abs(S_tr/n) > max_cutoff_ind, S_tr, 0)
    dot_product_ind = np.dot(X_tr_T_y_tr_cutoff, S_test)
    norm_ind = np.linalg.norm(T @ X_tr_T_y_tr_cutoff) * np.linalg.norm(y_t)
    R_squared_ind_test_values.append((dot_product_ind / norm_ind)**2)

# Define markers for each kappa value
kappa_markers = {0.05: 'o', 0.5: 's', 0.9: '^'}  # Circle, Square, Triangle

# Plotting R_squared_sum_test_values versus R_squared_ind_test_values with different markers for each kappa
plt.figure(figsize=(8, 8))

# Loop over each unique kappa and plot with its designated marker
for kappa, marker in kappa_markers.items():
    # Filter values corresponding to the current kappa
    kappa_filter = results_df['kappa'] == kappa
    R_squared_sum_values = np.array(R_squared_sum_test_values)[kappa_filter]
    R_squared_ind_values = np.array(R_squared_ind_test_values)[kappa_filter]

    plt.scatter(
        R_squared_sum_values,
        R_squared_ind_values,
        marker=marker,
        label=f'$\kappa$ = {kappa}',
        alpha=0.8
    )

# Add a dashed diagonal line y = x
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')

# Set limits for x and y axes
plt.xlim(0, 0.65)
plt.ylim(0, 0.65)

# Add labels and title
plt.xlabel('Testing $R^2_{sum}$')
plt.ylabel('Testing $R^2_{ind}$')
plt.legend()
plt.grid(True)

# Define save path and save the plot as a PDF
save_path = os.path.join('/path/to/summary_training/optimal_R2_marginal_w_retain.pdf')
plt.savefig(save_path, format='pdf')

# Show plot
plt.show()
