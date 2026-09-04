#!/usr/bin/env python3
"""
Chat follow-up conservation analysis for CTCF using updated_alignment_0611.fas.

Implements the requested changes:
1) Class conservation uses all species in class as denominator (gap-aware).
2) Human-mapped table (1..N human residues) for cleaner visualization.
3) Zinc-finger and YXF highlights.
4) Vertebrate - invertebrate conservation delta plots in:
   - full MSA coordinate space
   - human coordinate space
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IN_FASTA = Path(
    r"C:\Users\Nat\Downloads\updated_alignment_0611.fas"
)
OUT_DIR = Path(
    r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\alignments_out\CTCF\07"
)

FEATURES = [
    (226, 228, "YXF"),
    (266, 288, "ZF1"),
    (294, 316, "ZF2"),
    (322, 345, "ZF3"),
    (351, 373, "ZF4"),
    (379, 401, "ZF5"),
    (407, 430, "ZF6"),
    (437, 460, "ZF7"),
    (467, 489, "ZF8"),
    (495, 517, "ZF9"),
    (523, 546, "ZF10"),
]

# Use a smaller smoothing window for the human-mapped figure to preserve local variation.
HUMAN_SMOOTH_WIN = 11

GAP_CHARS = {"-", ".", "X"}



def read_fasta(path: Path) -> List[Tuple[str, str]]:
    items: List[Tuple[str, str]] = []
    sid = None
    buf: List[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                items.append((sid, "".join(buf).upper()))
            sid = line[1:].strip()
            buf = []
        else:
            if sid is None:
                raise ValueError("Invalid FASTA: sequence encountered before first header")
            buf.append(line)
    if sid is not None:
        items.append((sid, "".join(buf).upper()))
    return items


def find_human(records: Sequence[Tuple[str, str]]) -> Optional[Tuple[str, str]]:
    for h, s in records:
        if "Homo_sapiens" in h or "Human" in h:
            return h, s
    return None


def class_conservation_gap_aware(col_residues: Sequence[str], class_size: int) -> float:
    """
    Percent conservation with denominator = total species in class.
    Numerator = count of most frequent non-gap residue.
    """
    nongap = [aa for aa in col_residues if aa not in GAP_CHARS]
    if not nongap:
        return 0.0
    max_count = max(Counter(nongap).values())
    return 100.0 * max_count / class_size


def moving_average(y: np.ndarray, win: int = 25) -> np.ndarray:
    s = pd.Series(y)
    return s.rolling(window=win, center=True, min_periods=max(3, win // 3)).mean().to_numpy()


def add_feature_spans_human(ax):
    for start, end, label in FEATURES:
        ax.axvspan(start, end, alpha=0.13, color="#fdae6b")
    # Label only once each to keep figure readable
    y_top = ax.get_ylim()[1]
    for start, end, label in FEATURES:
        x = (start + end) / 2.0
        ax.text(x, y_top * 0.98, label, ha="center", va="top", fontsize=7)

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    records = read_fasta(IN_FASTA)
    if not records:
        raise ValueError("No sequences found")

    lengths = {len(seq) for _, seq in records if seq}
    if len(lengths) != 1:
        raise ValueError(f"Alignment has unequal lengths: {sorted(lengths)}")
    aln_len = lengths.pop()

    headers = [h for h, _ in records]
    
    # Find Ciona_intestinalis to split vertebrate/invertebrate
    split_idx = None
    for idx, header in enumerate(headers):
        if "Ciona_intestinalis" in header:
            split_idx = idx
            break
    
    if split_idx is None:
        raise ValueError("Cannot find Ciona_intestinalis split point")
    
    vertebrate_records = records[:split_idx]
    invertebrate_records = records[split_idx:]
    if not vertebrate_records or not invertebrate_records:
        raise ValueError("Expected non-empty vertebrate and invertebrate groups")

    human_item = find_human(vertebrate_records) or find_human(records)
    if human_item is None:
        raise ValueError("Human sequence not found")
    human_header, human_seq = human_item

    v_seqs = [s for _, s in vertebrate_records]
    i_seqs = [s for _, s in invertebrate_records]
    all_seqs = v_seqs + i_seqs

    n_v = len(v_seqs)
    n_i = len(i_seqs)
    n_all = len(all_seqs)

    rows: List[Dict[str, object]] = []
    human_pos = 0
    human_pos_to_aln_col: Dict[int, int] = {}

    for col_idx in range(aln_len):
        h_aa = human_seq[col_idx]
        if h_aa not in GAP_CHARS:
            human_pos += 1
            human_pos_to_aln_col[human_pos] = col_idx + 1
            h_pos_val: Optional[int] = human_pos
        else:
            h_pos_val = None

        v_col = [s[col_idx] for s in v_seqs]
        i_col = [s[col_idx] for s in i_seqs]
        a_col = [s[col_idx] for s in all_seqs]

        v_cons = class_conservation_gap_aware(v_col, n_v)
        i_cons = class_conservation_gap_aware(i_col, n_i)
        a_cons = class_conservation_gap_aware(a_col, n_all)

        rows.append(
            {
                "alignment_col": col_idx + 1,
                "human_pos": h_pos_val,
                "human_residue": h_aa if h_aa not in GAP_CHARS else "",
                "cons_vert_gapaware": round(v_cons, 2),
                "cons_invert_gapaware": round(i_cons, 2),
                "cons_all_gapaware": round(a_cons, 2),
                "delta_vert_minus_invert": round(v_cons - i_cons, 2),
            }
        )

    msa_df = pd.DataFrame(rows)
    human_df = msa_df[msa_df["human_pos"].notna()].copy()
    human_df["human_pos"] = human_df["human_pos"].astype(int)

    msa_csv = OUT_DIR / "ctcf_chat_followup_gapaware_msa.csv"
    human_csv = OUT_DIR / "ctcf_chat_followup_gapaware_human727.csv"
    msa_df.to_csv(msa_csv, index=False)
    human_df.to_csv(human_csv, index=False)

    # Plot 1: Delta in full MSA coordinates
    x_msa = msa_df["alignment_col"].to_numpy(dtype=float)
    y_delta_msa = msa_df["delta_vert_minus_invert"].to_numpy(dtype=float)
    y_delta_msa_s = moving_average(y_delta_msa, win=41)

    fig, ax = plt.subplots(figsize=(15, 4.8))
    ax.plot(x_msa, y_delta_msa, lw=0.6, color="#9ecae1", alpha=0.7, label="Raw delta")
    ax.plot(x_msa, y_delta_msa_s, lw=1.5, color="#08519c", label="Smoothed delta")
    ax.axhline(0, color="#555", lw=0.9, alpha=0.8)

    # Add feature spans mapped from human positions to alignment columns
    for start, end, label in FEATURES:
        aln_cols = [human_pos_to_aln_col[p] for p in range(start, end + 1) if p in human_pos_to_aln_col]
        if not aln_cols:
            continue
        left, right = min(aln_cols), max(aln_cols)
        ax.axvspan(left, right, alpha=0.12, color="#fdae6b")

    ax.set_title("Vertebrate - Invertebrate conservation (gap-aware), full MSA")
    ax.set_xlabel("Alignment column")
    ax.set_ylabel("Delta conservation (%)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ctcf_chat_followup_delta_full_msa.png", dpi=190)
    plt.close(fig)

    # Plot 2: Human-mapped conservation + delta with feature highlights
    x_h = human_df["human_pos"].to_numpy(dtype=float)
    y_v = human_df["cons_vert_gapaware"].to_numpy(dtype=float)
    y_i = human_df["cons_invert_gapaware"].to_numpy(dtype=float)
    y_d = human_df["delta_vert_minus_invert"].to_numpy(dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7.2), sharex=True, gridspec_kw={"height_ratios": [2, 1]})

    axes[0].plot(x_h, y_v, color="#2ca25f", lw=1.2, label="Vertebrates")
    axes[0].plot(x_h, y_i, color="#de2d26", lw=1.2, label="Invertebrates")
    axes[0].plot(x_h, moving_average(y_v, HUMAN_SMOOTH_WIN), color="#006d2c", lw=1.7, alpha=0.9)
    axes[0].plot(x_h, moving_average(y_i, HUMAN_SMOOTH_WIN), color="#a50f15", lw=1.7, alpha=0.9)
    axes[0].set_ylabel("Conservation (%)")
    axes[0].set_ylim(0, 100)
    axes[0].set_title("Gap-aware class conservation mapped to human CTCF")
    axes[0].legend(frameon=False, loc="upper right")

    axes[1].plot(x_h, y_d, color="#3182bd", lw=1.0, alpha=0.75, label="Raw delta")
    axes[1].plot(x_h, moving_average(y_d, HUMAN_SMOOTH_WIN), color="#08519c", lw=1.7, label="Smoothed delta")
    axes[1].axhline(0, color="#555", lw=0.9, alpha=0.8)
    axes[1].set_ylabel("V - I (%)")
    axes[1].set_xlabel("Human CTCF residue position")
    axes[1].legend(frameon=False, loc="upper right")

    for ax in axes:
        add_feature_spans_human(ax)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "ctcf_chat_followup_human727_conservation_and_delta.png", dpi=190)
    plt.close(fig)

    # Horizontal summary table similar to your previous format, but gap-aware
    cols = [str(int(p)) for p in human_df["human_pos"].tolist()]
    summary = pd.DataFrame(
        [
            human_df["human_pos"].astype(int).tolist(),
            human_df["cons_vert_gapaware"].round(0).astype(int).tolist(),
            human_df["cons_invert_gapaware"].round(0).astype(int).tolist(),
            human_df["cons_all_gapaware"].round(0).astype(int).tolist(),
            human_df["delta_vert_minus_invert"].round(0).astype(int).tolist(),
        ],
        index=[
            "human_position",
            "pct_conserved_vertebrates_gapaware",
            "pct_conserved_invertebrates_gapaware",
            "pct_conserved_all_gapaware",
            "delta_vert_minus_invert_gapaware",
        ],
        columns=cols,
    )
    summary.to_csv(OUT_DIR / "ctcf_chat_followup_human727_summary_horizontal.csv")

    print("Wrote:")
    print(f"  {msa_csv}")
    print(f"  {human_csv}")
    print(f"  {OUT_DIR / 'ctcf_chat_followup_delta_full_msa.png'}")
    print(f"  {OUT_DIR / 'ctcf_chat_followup_human727_conservation_and_delta.png'}")
    print(f"  {OUT_DIR / 'ctcf_chat_followup_human727_summary_horizontal.csv'}")
    print(f"Counts: vertebrates={n_v}, invertebrates={n_i}, all={n_all}, human_residues={len(human_df)}")
    print(f"Human header used: {human_header}")


if __name__ == "__main__":
    main()
