import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Set figure properties
plt.rcParams['figure.figsize'] = (13, 8)  # Width, height in inches
plt.rcParams['figure.dpi'] = 300  # Dots per inch

# Disable LaTeX rendering, use MathText instead
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Georgia']
plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern for MathText
plt.rcParams['font.size'] = 30

# Parameters
n_values = [5000, 10000]  # Number of samples
p = 46000  # Number of features
kappa_values = [0.05]  # Probability that a coordinate of beta is zero
mean_beta = 0  # Mean of the Gaussian distribution for beta
std_beta = 5  # Standard deviation of the Gaussian distribution for beta
noise_std_values = [1]  # Standard deviation of the noise
n_w = 2000  # Number of rows for the external dataset W

# Generate the covariance matrix Sigma for X
Sigma = np.identity(p)

# Create an identity matrix I_p of size p x p
I_p = np.eye(p)

# Define a range of lambda values
lambdas = np.linspace(0.01, 5, 100)  # From 0.01 to 10
R_squared_ind = []
R_squared_sum = []

# Number of iterations for averaging
iterations = 10

# Loop through lambda values and compute average R^2 for each
# Iterate over n, kappa, and noise_std
for n in n_values:
    for kappa in kappa_values:
        for noise_std in noise_std_values:
            # Parameters
            n_tr = int(6/10 * n)  # Number of training samples
            n_t = int(n - n_tr)  # Number of testing samples

            # Generate the signal beta
            beta = np.zeros(p)
            non_zero_mask = np.random.rand(p) > kappa
            beta[non_zero_mask] = np.random.normal(mean_beta, std_beta/p, size=np.sum(non_zero_mask))

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
                    dot_product_sum = np.dot(ref_ridge_sum, S_t/n_t)
                    print(dot_product_sum)

                    # Compute the norm of ref_ridge_ind
                    norm_ref_ridge_sum = ( np.linalg.norm(Sigma @ ref_ridge_sum) * np.linalg.norm(y_t) )
                    print(norm_ref_ridge_sum)

                    # Compute the R^2 value for sum and store it
                    R_squared_sum_values.append( (dot_product_sum / norm_ref_ridge_sum)**2 )

                    # Compute the ridge estimator for S_tr
                    ref_ridge_ind = Skrinkage @ X_tr_T_y_tr

                    # Compute the dot product <ref_ridge_ind, X_t_T_y_t>
                    dot_product_ind = np.dot(ref_ridge_ind, X_t_T_y_t/n_t)
                    print("ind num:", dot_product_ind)

                    # Compute the norm of ref_ridge_ind
                    norm_ref_ridge_ind = ( np.linalg.norm(Sigma @ ref_ridge_ind) * np.linalg.norm(y_t) )
                    print("ind den:", norm_ref_ridge_ind)

                    # Compute the R^2 value for ind and store it
                    R_squared_ind_values.append( (dot_product_ind / norm_ref_ridge_ind)**2 )

                # break

                # Compute the average R^2 for this lambda
                R_squared_sum.append(np.mean(R_squared_sum_values))
                R_squared_ind.append(np.mean(R_squared_ind_values))


            plt.figure()
            plt.plot(lambdas, R_squared_sum, label='$R^2_{sum}$', color='royalblue', linestyle='-', marker='o')
            plt.plot(lambdas, R_squared_ind, label='$R^2_{ind}$', color='red', linestyle='--', marker='x')
            plt.xscale('linear')
            plt.xlabel('$\lambda$')
            plt.ylabel('$R^2$')
            plt.title(f'$R^2$ vs. $\lambda$ for n={n}, p={p}, $\kappa$={kappa}, noise_variance={noise_std}, beta_variance={std_beta}')
            plt.legend()
            plt.grid(True)
            filename = f'/path/to/summary_training/summery_n{n}_p_4k6_kappa{kappa}_sigmabeta_5_noise{noise_std}.png'
            plt.savefig(filename)
            plt.show()

            # Save the R_squared_sum_values to a text file
            txt_filename = f'/path/to/summary_training/R_squared_sum_n{n}_p_4k6_kappa{kappa}_sigmabeta_5_noise{noise_std}.txt'
            np.savetxt(txt_filename, R_squared_sum_values)

            # Save the R_squared_sum_values to a text file
            txt_filename = f'/path/to/summary_training/R_squared_ind_n{n}_p_4k6_kappa{kappa}_sigmabeta_5_noise{noise_std}.txt'
            np.savetxt(txt_filename, R_squared_ind_values)