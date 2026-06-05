import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse

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

# Disable LaTeX rendering, use MathText instead
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Georgia']
plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern for MathText
plt.rcParams['font.size'] = 30


# Parameters
# n_values = [5000, 10000, 50000]  # Number of samples
p = 2000  # Number of features
# kappa_values = [0.1, 0.2, 0.5]   # Probability that a coordinate of beta is zero
mean_beta = 0  # Mean of the Gaussian distribution for beta
std_beta = 2  # Standard deviation of the Gaussian distribution for beta
noise_std_values =  [1]  # Standard deviation of the noise
n_w = 1000  # Number of rows for the external dataset W

# Generate the covariance matrix Sigma for X
Sigma = np.identity(p)

# Create an identity matrix I_p of size p x p
I_p = np.identity(p)

# Number of iterations for averaging
iterations = 100

# Loop through lambda values and compute average R^2 for each
# Iterate over n, kappa, and noise_std
for noise_std in noise_std_values:
            # Parameters
            n_tr = int(6/10 * n)  # Number of training samples
            n_t = int(n - n_tr)  # Number of testing samples

            h = float(std_beta)**2/(float(std_beta)**2 + float(noise_std)**2)
            print("Herb:", h)

            # Generate the signal beta
            beta = np.zeros(p)
            non_zero_mask = np.random.rand(p) > kappa
            beta[non_zero_mask] = np.random.normal(mean_beta, std_beta/np.sqrt(p), size=np.sum(non_zero_mask))

            # Define a range of lambda values
            lambdas = np.linspace(0.001, 1, 100)  # From 0.01 to 10
            R_squared_ind = []
            R_squared_sum = []

            # Loop through lambda values and compute average R^2 for each
            for lam in lambdas:
                R_squared_sum_values = []
                R_squared_ind_values = []

                for _ in range(iterations):
                    # Generate the external dataset W
                    W = np.random.multivariate_normal(mean=np.zeros(p), cov=Sigma, size=n_w)

                    # Compute W^T W
                    W_T_W = W.T @ W

                    # Generate the design matrix X, where each row is Gaussian with covariance Sigma
                    X = np.random.multivariate_normal(mean=np.zeros(p), cov=Sigma, size=n)
                    # print(X[0])

                    # Generate the response variable y as y = X @ beta + noise
                    y = X @ beta + np.random.normal(0, noise_std, size=n)

                    # Split the data into training and testing sets
                    # X_tr, X_t, y_tr, y_t = train_test_split(X, y, train_size=n_tr, test_size=n_t, random_state=42)

                    # Shuffle indices
                    indices = np.arange(X.shape[0])
                    np.random.shuffle(indices)

                    # Split indices based on the given sizes
                    train_indices = indices[:n_tr]
                    test_indices = indices[n_tr:n]

                    # Select the corresponding data
                    X_tr = X[train_indices]
                    y_tr = y[train_indices]
                    X_t = X[test_indices]
                    y_t = y[test_indices]

                    # Compute summary statistics X^T y for training and testing sets
                    X_tr_T_y_tr = X_tr.T @ y_tr
                    X_t_T_y_t = X_t.T @ y_t
                    X_T_y = X.T @ y
                    # print(np.linalg.norm(X_T_y))

                    Cov_X_T_y = np.outer(X_T_y - n * Sigma @ beta, X_T_y - n * Sigma @ beta)
                    Cov_X_T_y = 0.5 * (Cov_X_T_y + Cov_X_T_y.T)
                    # print(Cov_X_T_y.shape)
                    mean = (n_tr / n) * X_T_y  # Mean vector of zeros

                    # Compute S^(tr)
                    S_tr = np.random.multivariate_normal(mean, ( (n_tr * (n - n_tr)) / n**2 ) * Cov_X_T_y)

                    # Compute S^(t)
                    S_t = X_T_y - S_tr

                    # Output the results
                    # print("S^(tr):", S_tr)
                    # print("S^(t):", S_t)

                    Skrinkage = np.linalg.inv(W_T_W + n_w * lam * I_p)

                    # Compute the ridge estimator for S_tr
                    ref_ridge_sum = Skrinkage @ S_tr

                    # Compute the dot product <ref_ridge_ind, X_t_T_y_t>
                    dot_product_sum = np.dot(ref_ridge_sum, S_t)
                    # print("sum num:", dot_product_sum)

                    # Compute the norm of ref_ridge_ind
                    norm_ref_ridge_sum = ( np.linalg.norm(X_t @ ref_ridge_sum) * np.linalg.norm(y_t) )
                    # print("sum den:", norm_ref_ridge_sum)

                    # Compute the R^2 value for sum and store it
                    R_squared_sum_values.append( (dot_product_sum / norm_ref_ridge_sum)**2 )



                    # Compute the ridge estimator for S_tr
                    ref_ridge_ind = Skrinkage @ X_tr_T_y_tr

                    # Compute the dot product <ref_ridge_ind, X_t_T_y_t>
                    dot_product_ind = np.dot(ref_ridge_ind, X_t_T_y_t)
                    # print("ind num:", dot_product_ind)

                    # Compute the norm of ref_ridge_ind
                    norm_ref_ridge_ind = ( np.linalg.norm(X_t @ ref_ridge_ind) * np.linalg.norm(y_t) )
                    # print("ind den:", norm_ref_ridge_ind)

                    # Compute the R^2 value for ind and store it
                    R_squared_ind_values.append( (dot_product_ind / norm_ref_ridge_ind)**2 )


                # Compute the average R^2 for this lambda
                R_squared_sum.append(np.mean(R_squared_sum_values))
                print("sum R^2:", np.mean(R_squared_sum_values))

                R_squared_ind.append(np.mean(R_squared_ind_values))
                print("ind R^2:", np.mean(R_squared_ind_values))

            
            # Ensure the lists are numpy arrays for easier manipulation
            R_squared_sum = np.array(R_squared_sum)
            R_squared_ind = np.array(R_squared_ind)

            # Create the plot
            plt.figure()

            # Scatter plot of R_squared_sum vs R_squared_ind
            plt.scatter(R_squared_sum, R_squared_ind, color='royalblue', s=7)  # You can adjust color and size as needed

            # Plot the 45-degree reference line (y=x)
            max_val = max(np.max(R_squared_sum), np.max(R_squared_ind))
            plt.plot([0, max_val], [0, max_val], 'k--', lw=1)

            # Labels and title
            plt.xlabel(r'$R^2$ (sum approach)', fontsize=12)
            plt.ylabel(r'$R^2$ (individual approach)', fontsize=12)
            plt.title('Comparison of $R^2$ values', fontsize=14)

            # Set equal scaling for both axes
            plt.xlim([0, max_val])
            plt.ylim([0, max_val])

            # Show the plot
            plt.show()  

            plt.figure()
            plt.plot(lambdas, R_squared_sum, label='$R^2_{sum}$', color='royalblue', linestyle='-', marker='o')
            plt.plot(lambdas, R_squared_ind, label='$R^2_{ind}$', color='red', linestyle='--', marker='x')
            plt.xscale('linear')
            plt.xlabel('$\lambda$')
            plt.ylabel('$R^2$')
            plt.title(f'$R^2$ vs. $\lambda$ for n={n}, p={p}, $\kappa$={kappa}, \n noise_variance={noise_std}, beta_variance={std_beta}')
            plt.legend()
            plt.grid(True)
            filename = f'/path/to/summary_training/summery_n_{n}_p_{p}_kappa_{kappa}_sigmabeta_{std_beta}_noise_{noise_std}.png'
            plt.savefig(filename)
            plt.show()