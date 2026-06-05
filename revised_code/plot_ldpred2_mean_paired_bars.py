from pathlib import Path
import os
import re
import tempfile

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "Results" / "Summery_Validation"
FIGURES_DIR = BASE_DIR / "Figures"

EXPECTED_PAIR_COUNT = 71
TOP_N = 10

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg-cache"))

import matplotlib.pyplot as plt


def collect_files(pattern: str) -> dict[int, Path]:
    regex = re.compile(pattern)
    files = {}

    for path in RESULTS_DIR.glob("*.csv"):
        match = regex.fullmatch(path.name)
        if match is None:
            continue

        trait = int(match.group(1))
        if trait in files:
            raise ValueError(f"Duplicate file for trait {trait}: {files[trait]} and {path}")
        files[trait] = path

    return files


def r_squared_mean(path: Path) -> float:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    values = pd.to_numeric(df["r_squared"], errors="coerce").dropna()
    return values.mean()


def build_mean_table(
    comparison_label: str,
    ind_files: dict[int, Path],
    sum_files: dict[int, Path],
    sum_label: str,
) -> pd.DataFrame:
    ind_traits = set(ind_files)
    sum_traits = set(sum_files)

    missing_sum = sorted(ind_traits - sum_traits)
    missing_ind = sorted(sum_traits - ind_traits)
    if missing_sum or missing_ind:
        raise ValueError(
            f"{comparison_label} trait mismatch. "
            f"Missing sum traits: {missing_sum}; missing ind traits: {missing_ind}"
        )

    if len(ind_traits) != EXPECTED_PAIR_COUNT:
        raise ValueError(
            f"{comparison_label} has {len(ind_traits)} pairs, expected {EXPECTED_PAIR_COUNT}"
        )

    rows = []
    for trait in sorted(ind_traits):
        ind_mean = r_squared_mean(ind_files[trait])
        sum_mean = r_squared_mean(sum_files[trait])
        rows.append(
            {
                "trait": trait,
                "LDpred2_ind_mean": ind_mean,
                f"{sum_label}_mean": sum_mean,
                "pair_max_mean": max(ind_mean, sum_mean),
            }
        )

    return pd.DataFrame(rows)


def print_mean_table(comparison_label: str, table: pd.DataFrame, sum_label: str) -> None:
    sum_mean_col = f"{sum_label}_mean"
    print(f"\nAll mean pairs: {comparison_label}")
    print(
        table[["trait", "LDpred2_ind_mean", sum_mean_col, "pair_max_mean"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )


def plot_top_mean_pairs(
    table: pd.DataFrame,
    comparison_label: str,
    sum_label: str,
    output_file: Path,
) -> None:
    sum_mean_col = f"{sum_label}_mean"
    top_df = table.nlargest(TOP_N, "pair_max_mean").sort_values(
        "pair_max_mean", ascending=False
    )

    x = range(len(top_df))
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(
        [position - width / 2 for position in x],
        top_df["LDpred2_ind_mean"],
        width=width,
        label="LDpred2-pesudo",
    )
    ax.bar(
        [position + width / 2 for position in x],
        top_df[sum_mean_col],
        width=width,
        label="LDpred2",
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels([str(trait) for trait in top_df["trait"]])
    ax.set_xlabel("Trait ID")
    ax.set_ylabel("Mean R-Squared")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

    print(f"\nTop {TOP_N} mean pairs for {comparison_label}:")
    print(
        top_df[["trait", "LDpred2_ind_mean", sum_mean_col, "pair_max_mean"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )
    print(f"Saved mean paired bar plot to: {output_file}")


def main() -> None:
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 15

    ind_files = collect_files(r"LDpred2_ind_(\d+)\.csv")
    plain_sum_files = collect_files(r"LDpred2_sum_(\d+)\.csv")
    val_100_sum_files = collect_files(r"LDpred2_sum_(\d+)_val_100\.csv")

    plain_table = build_mean_table(
        comparison_label="LDpred2_ind vs LDpred2_sum",
        ind_files=ind_files,
        sum_files=plain_sum_files,
        sum_label="LDpred2_sum",
    )
    print_mean_table("LDpred2_ind vs LDpred2_sum", plain_table, "LDpred2_sum")
    plot_top_mean_pairs(
        table=plain_table,
        comparison_label="LDpred2_ind vs LDpred2_sum",
        sum_label="LDpred2_sum",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_top10_mean_paired_barplot.pdf",
    )

    val_100_table = build_mean_table(
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        ind_files=ind_files,
        sum_files=val_100_sum_files,
        sum_label="LDpred2_sum_val_100",
    )
    print_mean_table(
        "LDpred2_ind vs LDpred2_sum_val_100",
        val_100_table,
        "LDpred2_sum_val_100",
    )
    plot_top_mean_pairs(
        table=val_100_table,
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        sum_label="LDpred2_sum_val_100",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_val_100_top10_mean_paired_barplot.pdf",
    )


if __name__ == "__main__":
    main()
