import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd


CTCF_ZF_RANGES = [
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

# Lysine-rich linker between the YDF motif and the first zinc finger.
CTCF_LYSINE_LINKER = (229, 265)
DEFAULT_HIGHLIGHT_COLORS = ["#c0392b", "#16a085", "#7f8c8d", "#d35400"]


def find_motif_positions(df: pd.DataFrame, motif: str) -> list[tuple[int, int]]:
    seq = "".join(df["AA"].tolist())
    hits = []
    start = 0
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        pos_start = int(df.iloc[idx]["POS"])
        pos_end = int(df.iloc[idx + len(motif) - 1]["POS"])
        hits.append((pos_start, pos_end))
        start = idx + 1
    return hits


def find_c2h2_like_regions(df: pd.DataFrame) -> list[tuple[int, int, str]]:
    # Canonical C2H2-like motif: C-X(2-4)-C-X(10-14)-H-X(2-5)-H
    seq = "".join(df["AA"].tolist())
    pos = df["POS"].astype(int).tolist()
    pattern = re.compile(r"C.{2,4}C.{10,14}H.{2,5}H")

    zf = []
    for i, match in enumerate(pattern.finditer(seq), start=1):
        s_idx, e_idx = match.start(), match.end() - 1
        zf.append((pos[s_idx], pos[e_idx], f"ZF{i}"))
    return zf


def infer_sequence_features(df: pd.DataFrame):
    ydf_hits = find_motif_positions(df, "YDF")
    zf_ranges = find_c2h2_like_regions(df)

    linker = None
    if ydf_hits and zf_ranges:
        ydf_end = max(e for _, e in ydf_hits)
        zf1_start = zf_ranges[0][0]
        if ydf_end + 1 <= zf1_start - 1:
            linker = (ydf_end + 1, zf1_start - 1)

    return ydf_hits, zf_ranges, linker


def parse_grades(grades_path: Path) -> pd.DataFrame:
    rows = []
    with grades_path.open("r", encoding="utf-8") as handle:
        in_table = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
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

            pos = int(parts[0])
            aa = parts[1]
            score = float(parts[2])
            color = int(parts[3].replace("*", ""))
            low_conf = "*" in parts[3]
            rows.append(
                {
                    "POS": pos,
                    "AA": aa,
                    "SCORE": score,
                    "COLOR": color,
                    "LOW_CONF": low_conf,
                }
            )

    if not rows:
        raise ValueError("No data rows were parsed from the ConSurf grades file.")

    df = pd.DataFrame(rows).sort_values("POS").reset_index(drop=True)
    return df


def parse_highlight_regions(raw_regions: list[str]) -> list[dict]:
    regions = []
    for index, raw in enumerate(raw_regions):
        parts = raw.split("|")
        if len(parts) not in {3, 4}:
            raise ValueError(
                "Each --highlight-region value must use 'label|start|end' or 'label|start|end|color'"
            )

        label = parts[0].strip()
        start = int(parts[1])
        end = int(parts[2])
        color = parts[3].strip() if len(parts) == 4 else DEFAULT_HIGHLIGHT_COLORS[index % len(DEFAULT_HIGHLIGHT_COLORS)]

        if start > end:
            raise ValueError(f"Invalid highlight region '{raw}': start must be <= end")

        regions.append(
            {
                "label": label,
                "start": start,
                "end": end,
                "color": color,
            }
        )

    return regions


def build_plot(
    df: pd.DataFrame,
    output_path: Path,
    title: str,
    annotation_mode: str,
    linker_override: tuple[int, int] | None,
    highlight_regions: list[dict],
) -> None:
    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(16, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.5, 1]},
        constrained_layout=True,
    )

    # Top panel: normalized score across sequence.
    ax1.plot(df["POS"], df["SCORE"], linewidth=1.5)
    ax1.axhline(0.0, linestyle="--", linewidth=1)
    ax1.set_ylabel("Normalized Score")
    ax1.set_title(title)
    ax1.grid(True, alpha=0.25)

    # Highlight features either from fixed human coordinates or inferred from sequence.
    if annotation_mode == "sequence":
        ydf_hits, zf_ranges, linker = infer_sequence_features(df)
    else:
        ydf_hits = find_motif_positions(df, "YDF")
        zf_ranges = CTCF_ZF_RANGES
        linker = CTCF_LYSINE_LINKER

    if linker_override is not None:
        linker = linker_override

    for i, (start, end) in enumerate(ydf_hits):
        label = "YDF motif" if i == 0 else None
        ax1.axvspan(start - 0.5, end + 0.5, alpha=0.30, color="#f39c12", label=label)

    for i, (start, end, name) in enumerate(zf_ranges):
        label = "Zinc fingers" if i == 0 else None
        ax1.axvspan(start - 0.5, end + 0.5, alpha=0.10, color="#2980b9", label=label)
        if i in {0, 5, 10}:
            ax1.text(
                (start + end) / 2,
                ax1.get_ylim()[1] * 0.92,
                name,
                ha="center",
                va="top",
                fontsize=8,
                color="#1f4e79",
            )

    linker_start = None
    linker_end = None
    if linker is not None:
        linker_start, linker_end = linker
        ax1.axvspan(
            linker_start - 0.5,
            linker_end + 0.5,
            alpha=0.16,
            color="#8e44ad",
            label="Lysine linker",
        )

    low_conf = df[df["LOW_CONF"]]
    if not low_conf.empty:
        ax1.scatter(
            low_conf["POS"],
            low_conf["SCORE"],
            marker="x",
            s=28,
            linewidths=0.9,
            label="Low confidence",
        )

    for idx, region in enumerate(highlight_regions):
        ax1.axvspan(
            region["start"] - 0.5,
            region["end"] + 0.5,
            alpha=0.18,
            color=region["color"],
            label=region["label"],
        )
        ax2.axvspan(
            region["start"] - 0.5,
            region["end"] + 0.5,
            alpha=0.18,
            color=region["color"],
        )
        ax1.text(
            (region["start"] + region["end"]) / 2,
            ax1.get_ylim()[1] * (0.84 - idx * 0.06),
            region["label"],
            ha="center",
            va="top",
            fontsize=8,
            color=region["color"],
        )

    ax1.legend(loc="best", frameon=False, ncol=2)

    # Bottom panel: ConSurf color grade (1 variable -> 9 conserved).
    cmap = plt.get_cmap("viridis", 9)
    colors = [cmap(v - 1) for v in df["COLOR"]]
    ax2.scatter(df["POS"], df["COLOR"], c=colors, s=18)
    for start, end in ydf_hits:
        ax2.axvspan(start - 0.5, end + 0.5, alpha=0.30, color="#f39c12")
    for start, end, _ in zf_ranges:
        ax2.axvspan(start - 0.5, end + 0.5, alpha=0.10, color="#2980b9")
    if linker_start is not None and linker_end is not None:
        ax2.axvspan(linker_start - 0.5, linker_end + 0.5, alpha=0.16, color="#8e44ad")
    ax2.set_ylabel("Grade")
    ax2.set_xlabel("Residue Position")
    ax2.set_yticks(range(1, 10))
    ax2.set_ylim(0.5, 9.5)
    ax2.grid(True, axis="y", alpha=0.25)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot ConSurf .grades output")
    parser.add_argument("grades_file", type=Path, help="Path to Human_CTCF_consurf.grades")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output image path (default: <grades_file_stem>_plot.png)",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="ConSurf Conservation Profile",
        help="Plot title",
    )
    parser.add_argument(
        "--annotation-mode",
        choices=["human", "sequence"],
        default="human",
        help="Use fixed human coordinates or infer annotations from plotted sequence.",
    )
    parser.add_argument(
        "--linker-start",
        type=int,
        default=None,
        help="Optional linker start residue (overrides inferred/default linker).",
    )
    parser.add_argument(
        "--linker-end",
        type=int,
        default=None,
        help="Optional linker end residue (overrides inferred/default linker).",
    )
    parser.add_argument(
        "--highlight-region",
        action="append",
        default=[],
        help="Add a highlighted region using 'label|start|end' or 'label|start|end|color'. Repeatable.",
    )
    args = parser.parse_args()

    linker_override = None
    if args.linker_start is not None or args.linker_end is not None:
        if args.linker_start is None or args.linker_end is None:
            raise ValueError("Provide both --linker-start and --linker-end together")
        if args.linker_start > args.linker_end:
            raise ValueError("--linker-start must be <= --linker-end")
        linker_override = (args.linker_start, args.linker_end)

    out = args.output or args.grades_file.with_name(args.grades_file.stem + "_plot.png")
    df = parse_grades(args.grades_file)
    highlight_regions = parse_highlight_regions(args.highlight_region)
    build_plot(df, out, args.title, args.annotation_mode, linker_override, highlight_regions)
    print(f"Saved plot: {out}")
    print(f"Plotted residues: {len(df)}")


if __name__ == "__main__":
    main()
