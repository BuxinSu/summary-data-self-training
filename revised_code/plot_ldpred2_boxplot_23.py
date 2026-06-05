from pathlib import Path
import os
import tempfile

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "Results" / "Summery_Validation"
FIGURES_DIR = BASE_DIR / "Figures"

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg-cache"))

import matplotlib.pyplot as plt

INPUT_FILES = {
    "LDpred2_ind_23": RESULTS_DIR / "LDpred2_ind_23.csv",
    "LDpred2_sum_23": RESULTS_DIR / "LDpred2_sum_23.csv",
}

OUTPUT_FILE = FIGURES_DIR / "LDpred2_23_r_squared_boxplot.pdf"


def load_r_squared(label: str, path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    return pd.DataFrame(
        {
            "method": label,
            "r_squared": pd.to_numeric(df["r_squared"], errors="coerce"),
        }
    ).dropna(subset=["r_squared"])


def main() -> None:
    plot_df = pd.concat(
        [load_r_squared(label, path) for label, path in INPUT_FILES.items()],
        ignore_index=True,
    )

    print("Sample SD of r_squared:")
    for method, group in plot_df.groupby("method", sort=False):
        print(f"{method}: {group['r_squared'].std(ddof=1):.10f}")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    labels = list(INPUT_FILES)
    values = [
        plot_df.loc[plot_df["method"] == label, "r_squared"].to_numpy()
        for label in labels
    ]

    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(values, tick_labels=labels, patch_artist=True, showmeans=True)
    ax.set_ylabel(r"$R^2$")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTPUT_FILE, bbox_inches="tight")
    print(f"Boxplot saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
