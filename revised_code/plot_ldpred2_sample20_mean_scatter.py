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
SAMPLE_N = 20
RANDOM_SEED = 20260422

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


def sample_random_state(path: Path) -> int:
    seed_text = f"{RANDOM_SEED}:{path.name}"
    digest = hashlib.sha256(seed_text.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def sampled_r_squared_mean(path: Path) -> float:
    df = pd.read_csv(path)
    if "r_squared" not in df.columns:
        raise ValueError(f"{path} does not contain an 'r_squared' column")

    values = pd.to_numeric(df["r_squared"], errors="coerce").dropna()
    if len(values) < SAMPLE_N:
        raise ValueError(f"{path} has only {len(values)} valid r_squared rows; need {SAMPLE_N}")

    sampled_values = values.sample(
        n=SAMPLE_N,
        replace=False,
        random_state=sample_random_state(path),
    )
    return sampled_values.mean()


def build_sampled_mean_table(
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
        rows.append(
            {
                "trait": trait,
                "LDpred2_ind_mean": sampled_r_squared_mean(ind_files[trait]),
                f"{sum_label}_mean": sampled_r_squared_mean(sum_files[trait]),
            }
        )

    return pd.DataFrame(rows)


def print_mean_table(comparison_label: str, table: pd.DataFrame, sum_label: str) -> None:
    sum_mean_col = f"{sum_label}_mean"
    print(f"\nAll sample-20 mean pairs: {comparison_label}")
    print(
        table[["trait", sum_mean_col, "LDpred2_ind_mean"]].to_string(
            index=False,
            float_format=lambda value: f"{value:.10f}",
        )
    )


def plot_mean_scatter(
    table: pd.DataFrame,
    comparison_label: str,
    sum_label: str,
    output_file: Path,
    x_axis_label: str,
) -> None:
    # Labels stay in the original orientation, but points are swapped as requested.
    x_col = f"{sum_label}_mean"
    y_col = "LDpred2_ind_mean"

    x = table[y_col]
    y = table[x_col]
    axis_max = max(x.max(), y.max()) * 1.08

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, s=14, color="royalblue", alpha=0.9, edgecolors="none")
    ax.plot([0, axis_max], [0, axis_max], linestyle="--", color="black", linewidth=1)

    ax.set_xlim(0, axis_max)
    ax.set_ylim(0, axis_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(r"LDpred2 $R^2_{\mathrm{ind}}$")
    ax.grid(False)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved sample-20 mean scatter plot for {comparison_label} to: {output_file}")


def main() -> None:
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.size"] = 25
    plt.rcParams["mathtext.fontset"] = "cm"

    ind_files = collect_files(r"LDpred2_ind_(\d+)\.csv")
    plain_sum_files = collect_files(r"LDpred2_sum_(\d+)\.csv")
    val_100_sum_files = collect_files(r"LDpred2_sum_(\d+)_val_100\.csv")

    plain_table = build_sampled_mean_table(
        comparison_label="LDpred2_ind vs LDpred2_sum",
        ind_files=ind_files,
        sum_files=plain_sum_files,
        sum_label="LDpred2_sum",
    )
    print_mean_table("LDpred2_ind vs LDpred2_sum", plain_table, "LDpred2_sum")
    plot_mean_scatter(
        table=plain_table,
        comparison_label="LDpred2_ind vs LDpred2_sum",
        sum_label="LDpred2_sum",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_mean_sample20_scatter.pdf",
        x_axis_label=r"LDpred2-pseudo $R^2_{\mathrm{sum}}$",
    )

    val_100_table = build_sampled_mean_table(
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
    plot_mean_scatter(
        table=val_100_table,
        comparison_label="LDpred2_ind vs LDpred2_sum_val_100",
        sum_label="LDpred2_sum_val_100",
        output_file=FIGURES_DIR / "LDpred2_ind_vs_sum_val_100_mean_sample20_scatter.pdf",
        x_axis_label=r"LDpred2-pseudo $R^2_{\mathrm{sum,val100}}$",
    )


if __name__ == "__main__":
    main()
