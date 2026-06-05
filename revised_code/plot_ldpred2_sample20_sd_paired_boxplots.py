from pathlib import Path
import hashlib
import os
import re
import tempfile

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "Results" / "Summery_Validation"
FIGURES_DIR = BASE_DIR / "Figures"

EXPECTED_PAIR_COUNT = 71
TOP_N = 10
SAMPLE_N = 20
RANDOM_SEED = 20260422

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg-cache"))

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


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


def sample_random_state(path: Path) -> int:
    seed_text = f"{RANDOM_SEED}:{path.name}"
    digest = hashlib.sha256(seed_text.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def read_sampled_r_squared(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    values = pd.to_numeric(df["r_squared"], errors="coerce").dropna()
    if len(values) < SAMPLE_N:
        raise ValueError(f"{path} has only {len(values)} valid r_squared rows; need {SAMPLE_N}")

    return values.sample(n=SAMPLE_N, replace=False, random_state=sample_random_state(path))


def build_sampled_pair_table(
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
        ind_values = read_sampled_r_squared(ind_files[trait])
        sum_values = read_sampled_r_squared(sum_files[trait])
        ind_sd = ind_values.std(ddof=1)
        sum_sd = sum_values.std(ddof=1)

        rows.append(
            {
                "trait": trait,
                "LDpred2_ind_SD": ind_sd,
                f"{sum_label}_SD": sum_sd,
                "pair_max_SD": max(ind_sd, sum_sd),
                "ind_values": ind_values,
                "sum_values": sum_values,
            }
        )

    return pd.DataFrame(rows)


def print_sd_table(comparison_label: str, table: pd.DataFrame, sum_label: str) -> None:
    sum_sd_col = f"{sum_label}_SD"
    print(f"\nAll sample-20 SD pairs: {comparison_label}")
    print(
        table[["trait", "LDpred2_ind_SD", sum_sd_col, "pair_max_SD"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )


def plot_top_pair_boxplots(
    table: pd.DataFrame,
    comparison_label: str,
    sum_label: str,
    output_file: Path,
) -> None:
    sum_sd_col = f"{sum_label}_SD"
    top_df = table.nlargest(TOP_N, "pair_max_SD").sort_values("pair_max_SD", ascending=False)

    centers = list(range(1, len(top_df) + 1))
    offset = 0.18
    box_width = 0.30
    ind_positions = [center - offset for center in centers]
    sum_positions = [center + offset for center in centers]

    ind_data = [values.to_numpy() for values in top_df["ind_values"]]
    sum_data = [values.to_numpy() for values in top_df["sum_values"]]

    fig, ax = plt.subplots(figsize=(12, 6))
    ind_box = ax.boxplot(
        ind_data,
        positions=ind_positions,
        widths=box_width,
        patch_artist=True,
        showfliers=False,
    )
    sum_box = ax.boxplot(
        sum_data,
        positions=sum_positions,
        widths=box_width,
        patch_artist=True,
        showfliers=False,
    )

    ind_color = "#4C78A8"
    sum_color = "#F58518"
    for patch in ind_box["boxes"]:
        patch.set_facecolor(ind_color)
        patch.set_alpha(0.80)
    for patch in sum_box["boxes"]:
        patch.set_facecolor(sum_color)
        patch.set_alpha(0.80)

    ax.set_xticks(centers)
    ax.set_xticklabels([str(trait) for trait in top_df["trait"]])
    ax.set_xlabel("Trait ID")
    ax.set_ylabel("R-Squared")
    ax.legend(
        handles=[
            Patch(facecolor=ind_color, alpha=0.80, label="LDpred2-pesudo"),
            Patch(facecolor=sum_color, alpha=0.80, label="LDpred2"),
        ]
    )
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

    print(f"\nTop {TOP_N} sample-20 boxplot pairs for {comparison_label}:")
    print(
        top_df[["trait", "LDpred2_ind_SD", sum_sd_col, "pair_max_SD"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )
    print(f"Saved sample-20 paired boxplot to: {output_file}")


def main() -> None:
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"

    ind_files = collect_files(r"LDpred2_ind_(\d+)\.csv")
    plain_sum_files = collect_files(r"LDpred2_sum_(\d+)\.csv")
    val_100_sum_files = collect_files(r"LDpred2_sum_(\d+)_val_100\.csv")

    plain_table = build_sampled_pair_table(
        comparison_label="LDpred2_ind vs LDpred2_sum",
        ind_files=ind_files,
        sum_files=plain_sum_files,
        sum_label="LDpred2_sum",
    )
    print_sd_table("LDpred2_ind vs LDpred2_sum", plain_table, "LDpred2_sum")
    plot_top_pair_boxplots(
        table=plain_table,
        comparison_label="LDpred2_ind vs LDpred2_sum",
        sum_label="LDpred2_sum",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_top10_sd_sample20_paired_boxplot.pdf",
    )

    val_100_table = build_sampled_pair_table(
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        ind_files=ind_files,
        sum_files=val_100_sum_files,
        sum_label="LDpred2_sum_val_100",
    )
    print_sd_table(
        "LDpred2_ind vs LDpred2_sum_val_100",
        val_100_table,
        "LDpred2_sum_val_100",
    )
    plot_top_pair_boxplots(
        table=val_100_table,
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        sum_label="LDpred2_sum_val_100",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_val_100_top10_sd_sample20_paired_boxplot.pdf",
    )


if __name__ == "__main__":
    main()
