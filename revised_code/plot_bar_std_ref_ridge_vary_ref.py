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

# Load data (all iterations, lambda=10.0)
data_path = os.path.join(os.path.dirname(__file__), '..', 'results', 'ref_ridge',
                         'aggregated_ref_ridge_vary_ref_all_iter.csv')
df = pd.read_csv(data_path)

output_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(output_dir, exist_ok=True)

# Group by (n, n_w, p, kappa, noise_std) and compute std over iterations
group_cols = ['n', 'n_w', 'p', 'kappa', 'noise_std']
std_df = df.groupby(group_cols).agg(
    R_squared_ind_std=('R_squared_ind', 'std'),
    R_squared_sum_std=('R_squared_sum', 'std'),
).reset_index()

# Generate one plot per n_w
n_w_values = sorted(std_df['n_w'].unique())

for n_w in n_w_values:
    sub = std_df[std_df['n_w'] == n_w].copy()
    sub.sort_values(by=['n', 'p', 'kappa', 'noise_std'], inplace=True)

    # Create group labels
    sub['label'] = sub.apply(
        lambda row: (
            f"$n\\!={int(row['n'])}$\n"
            f"$p\\!={int(row['p'])}$\n"
            f"$\\kappa\\!={row['kappa']}$\n"
            f"$\\sigma\\!={row['noise_std']}$"
        ), axis=1
    )

    x = np.arange(len(sub))
    width = 0.35

    fig, ax = plt.subplots()
    bars_ind = ax.bar(x - width / 2, sub['R_squared_ind_std'].values, width,
                      label=r'$R^2_{\rm ind}$', color='royalblue', alpha=0.85)
    bars_sum = ax.bar(x + width / 2, sub['R_squared_sum_std'].values, width,
                      label=r'$R^2_{\rm sum}$', color='red', alpha=0.85)

    ax.set_ylabel(r'Standard Deviation')
    ax.set_title(r'Std of $R^2$ across iterations ($n_w = %d$)' % n_w)
    ax.set_xticks(x)
    ax.set_xticklabels(sub['label'].values, fontsize=14)
    ax.legend(fontsize=22)
    ax.grid(axis='y', linestyle='--', alpha=0.5)

    fig.tight_layout()
    save_path = os.path.join(output_dir, f'bar_std_ref_ridge_nw{n_w}.pdf')
    fig.savefig(save_path, format='pdf')
    plt.close(fig)
    print(f"Saved: {save_path}")

print("All plots saved.")
