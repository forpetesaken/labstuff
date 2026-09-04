from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


FULL_ALIGNMENT = Path(r"C:\Users\Nat\Downloads\updated_alignment_0611.fas")
VERT_GRADES = Path(
    r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output\ctcf_consurf_vertebrates\Human_CTCF_consurf.grades"
)
INVERT_GRADES = Path(
    r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output\ctcf_consurf_invertebrates\Ciona_intestinalis_consurf.grades"
)
OUT_DIR = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output")

YDF = (226, 228)
LYSINE_LINKER = (229, 265)
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

GAP = {"-", ".", "X"}


def read_fasta(path: Path):
    records = []
    header = None
    seq_parts = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip()
            seq_parts = []
        else:
            if header is None:
                raise ValueError("Invalid FASTA: sequence before header")
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def get_seq(records, key: str):
    for h, s in records:
        if key in h:
            return s
    raise ValueError(f"Could not find sequence containing '{key}'")


def parse_grades(path: Path) -> pd.DataFrame:
    rows = []
    in_table = False
    with path.open("r", encoding="utf-8") as handle:
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
            if len(parts) < 4 or not parts[0].isdigit():
                continue

            rows.append(
                {
                    "pos": int(parts[0]),
                    "aa": parts[1],
                    "score": float(parts[2]),
                    "grade": int(parts[3].replace("*", "")),
                    "low_conf": "*" in parts[3],
                }
            )

    if not rows:
        raise ValueError(f"No rows parsed from {path}")

    return pd.DataFrame(rows)


def build_position_map(human_seq: str, ciona_seq: str) -> pd.DataFrame:
    if len(human_seq) != len(ciona_seq):
        raise ValueError("Human and Ciona sequences must have equal alignment length")

    rows = []
    human_pos = 0
    ciona_pos = 0

    for aln_col, (h_aa, c_aa) in enumerate(zip(human_seq, ciona_seq), start=1):
        if h_aa not in GAP:
            human_pos += 1
        if c_aa not in GAP:
            ciona_pos += 1

        if h_aa in GAP or c_aa in GAP:
            continue

        rows.append(
            {
                "alignment_col": aln_col,
                "human_pos": human_pos,
                "human_aa": h_aa,
                "ciona_pos": ciona_pos,
                "ciona_aa": c_aa,
            }
        )

    return pd.DataFrame(rows)


def add_features(ax):
    ax.axvspan(YDF[0], YDF[1], color="#f39c12", alpha=0.28, label="YDF")
    ax.axvspan(LYSINE_LINKER[0], LYSINE_LINKER[1], color="#8e44ad", alpha=0.14, label="Lysine linker")

    for idx, (start, end, name) in enumerate(ZINC_FINGERS):
        ax.axvspan(start, end, color="#2980b9", alpha=0.10, label="Zinc fingers" if idx == 0 else None)
        if idx in {0, 5, 10}:
            ax.text((start + end) / 2.0, ax.get_ylim()[1] * 0.96, name, ha="center", va="top", fontsize=8)


def main() -> None:
    records = read_fasta(FULL_ALIGNMENT)
    human_seq = get_seq(records, "Homo_sapiens")
    ciona_seq = get_seq(records, "Ciona_intestinalis")

    map_df = build_position_map(human_seq, ciona_seq)
    v_df = parse_grades(VERT_GRADES).rename(
        columns={
            "pos": "human_pos",
            "aa": "human_grade_aa",
            "score": "human_score",
            "grade": "human_grade",
            "low_conf": "human_low_conf",
        }
    )
    i_df = parse_grades(INVERT_GRADES).rename(
        columns={
            "pos": "ciona_pos",
            "aa": "ciona_grade_aa",
            "score": "ciona_score",
            "grade": "ciona_grade",
            "low_conf": "ciona_low_conf",
        }
    )

    merged = map_df.merge(v_df, on="human_pos", how="inner").merge(i_df, on="ciona_pos", how="inner")
    merged["grade_delta_vert_minus_invert"] = merged["human_grade"] - merged["ciona_grade"]
    merged["score_delta_vert_minus_invert"] = merged["human_score"] - merged["ciona_score"]

    merged = merged.sort_values("human_pos").reset_index(drop=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUT_DIR / "ctcf_split_consurf_grade_comparison.csv"
    merged.to_csv(out_csv, index=False)

    x = merged["human_pos"].to_numpy()
    y_v = merged["human_grade"].to_numpy()
    y_i = merged["ciona_grade"].to_numpy()
    y_d = merged["grade_delta_vert_minus_invert"].to_numpy()

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
        constrained_layout=True,
    )

    axes[0].plot(x, y_v, color="#006d2c", lw=1.3, label="Vertebrates ConSurf grade (Human query)")
    axes[0].plot(x, y_i, color="#a50f15", lw=1.3, label="Invertebrates ConSurf grade (Ciona mapped)")
    axes[0].set_ylim(0.5, 9.5)
    axes[0].set_yticks(range(1, 10))
    axes[0].set_ylabel("ConSurf grade")
    axes[0].set_title("CTCF split ConSurf comparison: vertebrates vs invertebrates")
    axes[0].grid(True, alpha=0.2)
    add_features(axes[0])
    axes[0].legend(frameon=False, ncol=2, loc="upper right")

    axes[1].plot(x, y_d, color="#08519c", lw=1.2, label="Grade delta (V - I)")
    axes[1].axhline(0.0, color="#555", lw=0.9, alpha=0.9)
    axes[1].set_ylabel("V - I")
    axes[1].set_xlabel("Human CTCF residue position")
    axes[1].grid(True, alpha=0.2)
    add_features(axes[1])
    axes[1].legend(frameon=False, loc="upper right")

    out_png = OUT_DIR / "ctcf_split_consurf_grade_comparison.png"
    fig.savefig(out_png, dpi=220)
    plt.close(fig)

    print(f"Wrote CSV: {out_csv}")
    print(f"Wrote PNG: {out_png}")
    print(f"Rows compared: {len(merged)}")


if __name__ == "__main__":
    main()
