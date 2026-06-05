import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Set the text properties to use LaTeX fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30  # 30 for two pictures in a row
plt.rcParams['figure.figsize'] = (10, 10)  # width, height in inches
plt.rcParams['figure.dpi'] = 300  # dots per inch
plt.rcParams['font.family'] = 'serif'  # Use serif fonts to match LaTeX's default

# Load data
data_path = os.path.join(os.path.dirname(__file__), '..', 'results', 'covariance_recovery_avg_p5000_rho0.6.csv')
df = pd.read_csv(data_path)

output_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(output_dir, exist_ok=True)

n_w = df['n_w'].values

# --- Plot 1: Operator Norm ---
fig, ax = plt.subplots()
ax.errorbar(n_w, df['operator_norm_mean'], yerr=df['operator_norm_std'],
            fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8, color='tab:blue')
ax.set_xscale('log')
ax.set_xlabel(r'$n_w$')
ax.set_ylabel(r'Operator Norm')
ax.set_title(r'$\|\hat{\Sigma} - \Sigma\|_{\mathrm{op}}$ \ ($p=5000,\ \rho=0.6$)')
ax.tick_params(axis='both', which='major')
fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'covariance_recovery_operator_norm_p5000_rho0.6.pdf'))
plt.close(fig)

# --- Plot 2: Frobenius Norm ---
fig, ax = plt.subplots()
ax.errorbar(n_w, df['frobenius_norm_mean'], yerr=df['frobenius_norm_std'],
            fmt='s-', capsize=5, capthick=2, linewidth=2, markersize=8, color='tab:red')
ax.set_xscale('log')
ax.set_xlabel(r'$n_w$')
ax.set_ylabel(r'Frobenius Norm')
ax.set_title(r'$\|\hat{\Sigma} - \Sigma\|_{F}$ \ ($p=5000,\ \rho=0.6$)')
ax.tick_params(axis='both', which='major')
fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'covariance_recovery_frobenius_norm_p5000_rho0.6.pdf'))
plt.close(fig)

print("Figures saved to:", output_dir)
