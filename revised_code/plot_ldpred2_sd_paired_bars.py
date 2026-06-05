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


def r_squared_sd(path: Path) -> float:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    values = pd.to_numeric(df["r_squared"], errors="coerce").dropna()
    return values.std(ddof=1)


def build_sd_table(
    comparison_name: str,
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
            f"{comparison_name} trait mismatch. "
            f"Missing sum traits: {missing_sum}; missing ind traits: {missing_ind}"
        )

    if len(ind_traits) != EXPECTED_PAIR_COUNT:
        raise ValueError(
            f"{comparison_name} has {len(ind_traits)} pairs, expected {EXPECTED_PAIR_COUNT}"
        )

    rows = []
    for trait in sorted(ind_traits):
        ind_sd = r_squared_sd(ind_files[trait])
        sum_sd = r_squared_sd(sum_files[trait])
        rows.append(
            {
                "trait": trait,
                "LDpred2_ind_SD": ind_sd,
                f"{sum_label}_SD": sum_sd,
                "pair_max_SD": max(ind_sd, sum_sd),
            }
        )

    return pd.DataFrame(rows)


def print_sd_table(title: str, df: pd.DataFrame) -> None:
    print(f"\n{title}")
    print(df.to_string(index=False, float_format=lambda value: f"{value:.10f}"))


def plot_top_pairs(
    df: pd.DataFrame,
    comparison_label: str,
    sum_label: str,
    output_file: Path,
    sum_legend_label: str,
) -> None:
    sum_sd_col = f"{sum_label}_SD"
    top_df = df.nlargest(TOP_N, "pair_max_SD").sort_values("pair_max_SD", ascending=False)

    x = range(len(top_df))
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(
        [position - width / 2 for position in x],
        top_df["LDpred2_ind_SD"],
        width=width,
        label="LDpred2-pesudo",
    )
    ax.bar(
        [position + width / 2 for position in x],
        top_df[sum_sd_col],
        width=width,
        label=sum_legend_label,
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels([str(trait) for trait in top_df["trait"]])
    ax.set_xlabel("Trait ID")
    ax.set_ylabel("SD for R-Squared")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

    print(f"\nTop {TOP_N} plotted pairs for {comparison_label}:")
    print(
        top_df[["trait", "LDpred2_ind_SD", sum_sd_col, "pair_max_SD"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )
    print(f"Saved paired bar plot to: {output_file}")


def main() -> None:
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 15

    ind_files = collect_files(r"LDpred2_ind_(\d+)\.csv")
    plain_sum_files = collect_files(r"LDpred2_sum_(\d+)\.csv")
    val_100_sum_files = collect_files(r"LDpred2_sum_(\d+)_val_100\.csv")

    plain_df = build_sd_table(
        comparison_name="LDpred2_ind vs LDpred2_sum",
        ind_files=ind_files,
        sum_files=plain_sum_files,
        sum_label="LDpred2_sum",
    )
    print_sd_table("All SD pairs: LDpred2_ind vs LDpred2_sum", plain_df)
    plot_top_pairs(
        df=plain_df,
        comparison_label="LDpred2_ind vs LDpred2_sum",
        sum_label="LDpred2_sum",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_top10_sd_paired_barplot.pdf",
        sum_legend_label="LDpred2",
    )

    val_100_df = build_sd_table(
        comparison_name="LDpred2_ind vs LDpred2_sum_val_100",
        ind_files=ind_files,
        sum_files=val_100_sum_files,
        sum_label="LDpred2_sum_val_100",
    )
    print_sd_table("All SD pairs: LDpred2_ind vs LDpred2_sum_val_100", val_100_df)
    plot_top_pairs(
        df=val_100_df,
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        sum_label="LDpred2_sum_val_100",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_val_100_top10_sd_paired_barplot.pdf",
        sum_legend_label="LDpred2",
    )


if __name__ == "__main__":
    main()
