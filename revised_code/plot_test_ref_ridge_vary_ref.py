import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

# Set the text properties to use LaTeX fonts
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 30  # 30 for two pictures in a row
plt.rcParams['figure.figsize'] = (10, 10)  # width, height in inches
plt.rcParams['figure.dpi'] = 300  # dots per inch
plt.rcParams['font.family'] = 'serif'  # Use serif fonts to match LaTeX's default

# Load data
data_path = os.path.join(os.path.dirname(__file__), '..', 'results', 'ref_ridge',
                         'test_results_ref_ridge_vary_ref.csv')
results_df = pd.read_csv(data_path)

output_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(output_dir, exist_ok=True)

# Define markers for each kappa value
kappa_markers = {0.05: 'o', 0.5: 's', 0.9: '^'}  # Circle, Square, Triangle

# Define colors for each n value
n_colors = {7500: 'red', 12500: 'royalblue'}

# Plotting test_R_sum versus test_R_ind
plt.figure()

for i, row in results_df.iterrows():
    n = int(row['n'])
    kappa = row['kappa']

    color = n_colors.get(n, 'black')
    marker = kappa_markers.get(kappa, 'x')

    plt.scatter(
        row['test_R_sum'],
        row['test_R_ind'],
        marker=marker,
        color=color,
        s=100,
        alpha=0.8
    )

# Add a dashed diagonal line y = x
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')

# Set limits for x and y axes
plt.xlim(0, 0.65)
plt.ylim(0, 0.65)

# Add labels
plt.xlabel(r'$R^2_{\rm sum, R}$')
plt.ylabel(r'$R^2_{\rm ind, R}$')

# Legend for kappa values
kappa_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
           label=r'$\kappa=0.05$', markersize=10),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='black',
           label=r'$\kappa=0.5$', markersize=10),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='black',
           label=r'$\kappa=0.9$', markersize=10),
]

# Legend for n values
n_legend_elements = [
    Line2D([0], [0], color='red', lw=4, label=r'$n = 7500,\ p = 5000$'),
    Line2D([0], [0], color='royalblue', lw=4, label=r'$n = 12500,\ p = 5000$'),
]

# Add both legends
leg1 = plt.legend(handles=kappa_legend_elements, loc='upper left')
plt.gca().add_artist(leg1)
plt.legend(handles=n_legend_elements, loc='lower right')

plt.grid(True)
plt.tight_layout()

save_path = os.path.join(output_dir, 'testing_R_squared_ref_ridge_vary_ref.pdf')
plt.savefig(save_path, format='pdf')
plt.close()

print(f"Figure saved to: {save_path}")
