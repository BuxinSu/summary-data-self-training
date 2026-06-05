from pathlib import Path
import os
import re
import tempfile

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "Results" / "Summery_Validation"
FIGURES_DIR = BASE_DIR / "Figures"

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg-cache"))

import matplotlib.pyplot as plt

EXPECTED_PAIRS = 71
OUTPUT_FILE = FIGURES_DIR / "LDpred2_top10_sd_paired_barplot.pdf"


def find_numeric_suffix_files(prefix: str) -> dict[int, Path]:
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)\.csv$")
    files = {}

    for path in RESULTS_DIR.glob(f"{prefix}_*.csv"):
        match = pattern.match(path.name)
        if match:
            files[int(match.group(1))] = path

    return files


def r_squared_sd(path: Path) -> float:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    values = pd.to_numeric(df["r_squared"], errors="coerce").dropna()
    return values.std(ddof=1)


def compute_sd_table() -> pd.DataFrame:
    ind_files = find_numeric_suffix_files("LDpred2_ind")
    sum_files = find_numeric_suffix_files("LDpred2_sum")
    pair_ids = sorted(set(ind_files) & set(sum_files))

    if len(pair_ids) != EXPECTED_PAIRS:
        missing_ind = sorted(set(sum_files) - set(ind_files))
        missing_sum = sorted(set(ind_files) - set(sum_files))
        raise ValueError(
            f"Expected {EXPECTED_PAIRS} numeric LDpred2 pairs, found {len(pair_ids)}. "
            f"Missing ind for: {missing_ind}; missing sum for: {missing_sum}"
        )

    rows = []
    for pair_id in pair_ids:
        ind_sd = r_squared_sd(ind_files[pair_id])
        sum_sd = r_squared_sd(sum_files[pair_id])
        rows.append(
            {
                "pair": pair_id,
                "LDpred2_ind_sd": ind_sd,
                "LDpred2_sum_sd": sum_sd,
                "pair_max_sd": max(ind_sd, sum_sd),
            }
        )

    return pd.DataFrame(rows)


def print_sd_table(sd_table: pd.DataFrame) -> None:
    print("Sample SD of r_squared for all numeric LDpred2 pairs:")
    print(
        sd_table[["pair", "LDpred2_ind_sd", "LDpred2_sum_sd"]].to_string(
            index=False,
            formatters={
                "LDpred2_ind_sd": "{:.10f}".format,
                "LDpred2_sum_sd": "{:.10f}".format,
            },
        )
    )


def plot_top10(sd_table: pd.DataFrame) -> None:
    top10 = sd_table.nlargest(10, "pair_max_sd").sort_values("pair_max_sd", ascending=False)

    labels = [str(pair_id) for pair_id in top10["pair"]]
    x_positions = range(len(top10))
    bar_width = 0.38

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 15

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(
        [x - bar_width / 2 for x in x_positions],
        top10["LDpred2_ind_sd"],
        width=bar_width,
        label="LDpred2_ind",
    )
    ax.bar(
        [x + bar_width / 2 for x in x_positions],
        top10["LDpred2_sum_sd"],
        width=bar_width,
        label="LDpred2_sum",
    )

    ax.set_xticks(list(x_positions))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Pair")
    ax.set_ylabel("Sample SD of validation $R^2$")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTPUT_FILE, bbox_inches="tight")

    print("\nTop 10 pairs included in paired bar plot:")
    print(
        top10[["pair", "LDpred2_ind_sd", "LDpred2_sum_sd", "pair_max_sd"]].to_string(
            index=False,
            formatters={
                "LDpred2_ind_sd": "{:.10f}".format,
                "LDpred2_sum_sd": "{:.10f}".format,
                "pair_max_sd": "{:.10f}".format,
            },
        )
    )
    print(f"\nPaired bar plot saved to: {OUTPUT_FILE}")


def main() -> None:
    sd_table = compute_sd_table()
    print_sd_table(sd_table)
    plot_top10(sd_table)


if __name__ == "__main__":
    main()
