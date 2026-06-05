import numpy as np
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool, cpu_count

# Set up command-line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--n_value', type=int, required=True, help='Number of samples (n)')
parser.add_argument('--kappa_value', type=float, required=True, help='Probability that a coordinate of beta is zero')
args = parser.parse_args()

n = args.n_value
kappa = args.kappa_value

# Set figure properties
plt.rcParams['figure.figsize'] = (10, 10)  # Width, height in inches
plt.rcParams['figure.dpi'] = 300  # Dots per inch
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Georgia']
plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern for MathText
plt.rcParams['font.size'] = 30

# Parameters
p = 2000  # Number of features
mean_beta = 0  # Mean of the Gaussian distribution for beta
std_beta = 2  # Standard deviation of the Gaussian distribution for beta
noise_std_values = [1]  # Standard deviation of the noise
n_w = 1000  # Number of rows for the external dataset W
Sigma = np.identity(p)
I_p = np.identity(p)
iterations = 72

# Parallelize this function using multiprocessing
def run_single_iteration(lam, n, kappa, Sigma, I_p, p, mean_beta, std_beta, noise_std, n_w):
    n_tr = int(6/10 * n)  # Number of training samples
    n_t = int(n - n_tr)  # Number of testing samples
    h = float(std_beta)**2 / (float(std_beta)**2 + float(noise_std)**2)
    
    # Generate the signal beta
    beta = np.zeros(p)
    non_zero_mask = np.random.rand(p) > kappa
    beta[non_zero_mask] = np.random.normal(mean_beta, std_beta/np.sqrt(p), size=np.sum(non_zero_mask))
    
    # Generate the external dataset W
    W = np.random.multivariate_normal(mean=np.zeros(p), cov=Sigma, size=n_w)
    W_T_W = W.T @ W

    # Generate the design matrix X and response variable y
    X = np.random.multivariate_normal(mean=np.zeros(p), cov=Sigma, size=n)
    y = X @ beta + np.random.normal(0, noise_std, size=n)

    # Split the data into training and testing sets
    indices = np.arange(X.shape[0])
    np.random.shuffle(indices)
    train_indices = indices[:n_tr]
    test_indices = indices[n_tr:n]
    X_tr = X[train_indices]
    y_tr = y[train_indices]
    X_t = X[test_indices]
    y_t = y[test_indices]

    # Compute summary statistics X^T y for training and testing sets
    X_tr_T_y_tr = X_tr.T @ y_tr
    X_T_y = X.T @ y
    Cov_X_T_y = np.outer(X_T_y - n * Sigma @ beta, X_T_y - n * Sigma @ beta)
    Cov_X_T_y = 0.5 * (Cov_X_T_y + Cov_X_T_y.T)

    # Compute S^(tr)
    mean = (n_tr / n) * X_T_y
    S_tr = np.random.multivariate_normal(mean, ((n_tr * (n - n_tr)) / n**2) * Cov_X_T_y)

    # Compute S^(t)
    S_t = X_T_y - S_tr

    # Compute shrinkage
    Skrinkage = np.linalg.inv(W_T_W + n_w * lam * I_p)

    # Compute the ridge estimator for S_tr
    ref_ridge_sum = Skrinkage @ S_tr
    dot_product_sum = np.dot(ref_ridge_sum, S_t)
    norm_ref_ridge_sum = (np.linalg.norm(X_t @ ref_ridge_sum) * np.linalg.norm(y_t))
    R_squared_sum_value = (dot_product_sum / norm_ref_ridge_sum)**2

    # Compute the ridge estimator for S_tr (independent)
    ref_ridge_ind = Skrinkage @ X_tr_T_y_tr
    dot_product_ind = np.dot(ref_ridge_ind, X_t.T @ y_t)
    norm_ref_ridge_ind = (np.linalg.norm(X_t @ ref_ridge_ind) * np.linalg.norm(y_t))
    R_squared_ind_value = (dot_product_ind / norm_ref_ridge_ind)**2

    return R_squared_sum_value, R_squared_ind_value

# Run the parallelized loop
for noise_std in noise_std_values:
    lambdas = np.linspace(0.001, 1, 360)
    R_squared_sum = []
    R_squared_ind = []

    # Parallelize across lambda values using multiprocessing
    with Pool(cpu_count()) as pool:
        results = pool.starmap(run_single_iteration, [(lam, n, kappa, Sigma, I_p, p, mean_beta, std_beta, noise_std, n_w) for lam in lambdas])

    R_squared_sum, R_squared_ind = zip(*results)

    # Plot the results
    plt.figure()
    plt.scatter(R_squared_sum, R_squared_ind, color='royalblue', s=7)
    plt.plot([0, max(R_squared_sum)], [0, max(R_squared_ind)], 'k--', lw=1)
    plt.xlabel(r'$R^2$ (sum approach)', fontsize=12)
    plt.ylabel(r'$R^2$ (individual approach)', fontsize=12)
    plt.title(f'Comparison of $R^2$ values for n={n}, kappa={kappa}, noise_std={noise_std}', fontsize=14)
    plt.xlim([0, max(R_squared_sum)])
    plt.ylim([0, max(R_squared_ind)])
    plt.savefig(f'/path/to/summary_training/summery_n_{n}_p_{p}_kappa_{kappa}_noise_{noise_std}.png')
    plt.show()
