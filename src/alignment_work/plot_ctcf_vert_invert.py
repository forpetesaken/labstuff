import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Human CTCF feature coordinates in human residue numbering.
YDF = (226, 228)
ZINC_FINGERS = [
    (266, 288, "ZF1"),
    (294, 316, "ZF2"),
    (322, 344, "ZF3"),
    (350, 372, "ZF4"),
    (379, 401, "ZF5"),
    (416, 438, "ZF6"),
    (444, 466, "ZF7"),
    (472, 494, "ZF8"),
    (502, 524, "ZF9"),
    (530, 552, "ZF10"),
    (558, 580, "ZF11"),
]
LYSINE_LINKER = (229, 265)


def parse_consurf_grades(grades_file: Path) -> pd.DataFrame:
    rows = []
    in_table = False
    with grades_file.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not in_table:
                if line.strip().startswith("POS"):
                    in_table = True
                continue
            if not line.strip():
                continue
            if line.lstrip().startswith("*"):
                break
            parts = line.split()
            if len(parts) < 4:
                continue
            if not parts[0].isdigit():
                continue
            rows.append(
                {
                    "POS": int(parts[0]),
                    "CONSURF_GRADE": int(parts[3].replace("*", "")),
                }
            )
    if not rows:
        raise ValueError(f"Could not parse rows from {grades_file}")
    return pd.DataFrame(rows)


def conservation_to_grade(cons: np.ndarray) -> np.ndarray:
    # Map 0..100 conservation to 1..9 grade scale (9 = more conserved).
    vals = np.clip(cons, 0.0, 100.0)
    grades = np.rint(1.0 + 8.0 * (vals / 100.0))
    return np.clip(grades, 1, 9)


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    series = pd.Series(values)
    return series.rolling(window=window, center=True, min_periods=max(3, window // 3)).mean().to_numpy()


def add_feature_spans(ax: plt.Axes) -> None:
    ax.axvspan(YDF[0], YDF[1], color="#f39c12", alpha=0.28, label="YDF")
    ax.axvspan(LYSINE_LINKER[0], LYSINE_LINKER[1], color="#8e44ad", alpha=0.14, label="Lysine linker")

    for idx, (start, end, name) in enumerate(ZINC_FINGERS):
        ax.axvspan(start, end, color="#2980b9", alpha=0.10, label="Zinc fingers" if idx == 0 else None)
        if idx in {0, 5, 10}:
            ax.text((start + end) / 2.0, ax.get_ylim()[1] * 0.96, name, ha="center", va="top", fontsize=8)


def make_plot(df: pd.DataFrame, output_path: Path, title: str) -> None:
    x = df["human_pos"].to_numpy(dtype=float)
    y_v = df["cons_vert_gapaware"].to_numpy(dtype=float)
    y_i = df["cons_invert_gapaware"].to_numpy(dtype=float)
    y_d = df["delta_vert_minus_invert"].to_numpy(dtype=float)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
        constrained_layout=True,
    )

    # Top panel: class conservation curves.
    axes[0].plot(x, y_v, color="#2ca25f", lw=1.1, alpha=0.8, label="Vertebrates (raw)")
    axes[0].plot(x, y_i, color="#de2d26", lw=1.1, alpha=0.8, label="Invertebrates (raw)")
    axes[0].plot(x, moving_average(y_v, 25), color="#006d2c", lw=1.7, label="Vertebrates (smoothed)")
    axes[0].plot(x, moving_average(y_i, 25), color="#a50f15", lw=1.7, label="Invertebrates (smoothed)")
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("Conservation (%)")
    axes[0].set_title(title)
    axes[0].grid(True, alpha=0.2)
    add_feature_spans(axes[0])
    axes[0].legend(frameon=False, ncol=3, loc="upper right")

    # Bottom panel: vertebrate minus invertebrate delta.
    axes[1].plot(x, y_d, color="#6baed6", lw=0.9, alpha=0.75, label="Delta raw")
    axes[1].plot(x, moving_average(y_d, 25), color="#08519c", lw=1.8, label="Delta smoothed")
    axes[1].axhline(0.0, color="#555", lw=0.9, alpha=0.9)
    axes[1].set_ylabel("V - I (%)")
    axes[1].set_xlabel("Human CTCF residue position")
    axes[1].grid(True, alpha=0.2)
    add_feature_spans(axes[1])
    axes[1].legend(frameon=False, loc="upper right")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def make_grade_plot(
    df: pd.DataFrame,
    output_path: Path,
    title: str,
    consurf_df: pd.DataFrame | None = None,
) -> None:
    x = df["human_pos"].to_numpy(dtype=float)
    vert_grade = conservation_to_grade(df["cons_vert_gapaware"].to_numpy(dtype=float))
    invert_grade = conservation_to_grade(df["cons_invert_gapaware"].to_numpy(dtype=float))
    delta_grade = vert_grade - invert_grade

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
        constrained_layout=True,
    )

    axes[0].plot(x, vert_grade, color="#006d2c", lw=1.4, label="Vertebrates (grade-scaled)")
    axes[0].plot(x, invert_grade, color="#a50f15", lw=1.4, label="Invertebrates (grade-scaled)")

    if consurf_df is not None:
        merged = pd.merge(
            df[["human_pos"]],
            consurf_df,
            left_on="human_pos",
            right_on="POS",
            how="left",
        )
        axes[0].plot(
            merged["human_pos"],
            merged["CONSURF_GRADE"],
            color="#1f77b4",
            lw=1.1,
            alpha=0.9,
            label="ConSurf grade (full set)",
        )

    axes[0].set_ylim(0.5, 9.5)
    axes[0].set_yticks(range(1, 10))
    axes[0].set_ylabel("Grade (1-9)")
    axes[0].set_title(title)
    axes[0].grid(True, alpha=0.2)
    add_feature_spans(axes[0])
    axes[0].legend(frameon=False, ncol=2, loc="upper right")

    axes[1].plot(x, delta_grade, color="#08519c", lw=1.2, label="Grade delta (V - I)")
    axes[1].axhline(0.0, color="#555", lw=0.9, alpha=0.9)
    axes[1].set_ylabel("V - I")
    axes[1].set_xlabel("Human CTCF residue position")
    axes[1].grid(True, alpha=0.2)
    add_feature_spans(axes[1])
    axes[1].legend(frameon=False, loc="upper right")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot vertebrate vs invertebrate CTCF conservation")
    parser.add_argument("input_csv", type=Path, help="ctcf_chat_followup_gapaware_human727.csv path")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: ctcf_vert_vs_invert_annotated.png next to CSV)",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="CTCF Vertebrate vs Invertebrate Conservation (Gap-Aware)",
        help="Figure title",
    )
    parser.add_argument(
        "--use-grade-scale",
        action="store_true",
        help="Plot vertebrate/invertebrate comparison on 1-9 grade scale.",
    )
    parser.add_argument(
        "--consurf-grades",
        type=Path,
        default=None,
        help="Optional ConSurf .grades file to overlay actual grades.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    required = {"human_pos", "cons_vert_gapaware", "cons_invert_gapaware", "delta_vert_minus_invert"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    df = df[df["human_pos"].notna()].copy()
    df["human_pos"] = df["human_pos"].astype(int)

    if args.use_grade_scale:
        out = args.output or args.input_csv.with_name("ctcf_vert_vs_invert_grade_scaled.png")
        consurf_df = parse_consurf_grades(args.consurf_grades) if args.consurf_grades else None
        make_grade_plot(df, out, args.title, consurf_df)
    else:
        out = args.output or args.input_csv.with_name("ctcf_vert_vs_invert_annotated.png")
        make_plot(df, out, args.title)

    print(f"Saved plot: {out}")
    print(f"Rows plotted: {len(df)}")


if __name__ == "__main__":
    main()
