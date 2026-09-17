#!/usr/bin/env python3
"""
automated_conservationscore_perresidue.py
-----------------------------------------
Reads an aligned FASTA with an >INVERTEBRATES marker row and produces a
horizontal CSV (one column per alignment column) containing 4 rows:

    Row 1  human_position        : human protein residue number (1-indexed);
                                                                    blank where human has a gap in that column
  Row 2  pct_conserved_vert    : % conservation excluding gaps, vertebrates
  Row 3  pct_conserved_invert  : % conservation excluding gaps, invertebrates
  Row 4  pct_conserved_all     : % conservation excluding gaps, all species

All % values are rounded to the nearest whole number.

Usage
-----
Edit IN_FASTA and OUT_DIR below, then run:
    python automated_conservationscore_perresidue.py
"""

from pathlib import Path
from collections import Counter

IN_FASTA = Path(
    r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work"
    r"\alignments_out\CTCF\06\ctcfaln_witharctic.fas"
)
OUT_DIR = Path(
    r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work"
    r"\alignments_out\CTCF\06"
)

# ── helpers ──────────────────────────────────────────────────────────────────

def read_fasta(path: Path):
    """Return list of (header, sequence) in file order."""
    items = []
    sid = None
    buf = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                items.append((sid, "".join(buf)))
            sid = line[1:].strip()
            buf = []
        else:
            if sid is None:
                raise ValueError("Sequence data before first header")
            buf.append(line.upper())
    if sid is not None:
        items.append((sid, "".join(buf)))
    return items


def find_human(records):
    """Return the human sequence string, or None if not found."""
    for h, s in records:
        if "Homo_sapiens" in h or "Human" in h:
            return s
    return None


def pct_conserved(sequences, col_idx):
    """
    % conservation at alignment column col_idx across given sequences,
    omitting gaps (-, ., X).  Returns None if no non-gap residues present.
    """
    col = [s[col_idx] for s in sequences]
    residues = [aa for aa in col if aa not in {"-", ".", "X"}]
    if len(residues) < 2:
        return None
    counts = Counter(residues)
    return round(100.0 * max(counts.values()) / len(residues))


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if not IN_FASTA.exists():
        raise FileNotFoundError(f"FASTA not found: {IN_FASTA}")

    records = read_fasta(IN_FASTA)
    headers = [h for h, _ in records]

    # Split at INVERTEBRATES marker
    if "INVERTEBRATES" not in headers:
        raise ValueError(
            "No >INVERTEBRATES marker found in FASTA. "
            "Please add a blank sequence with that header to separate groups."
        )
    split_idx = headers.index("INVERTEBRATES")
    vert_seqs  = [s for _, s in records[:split_idx]]
    invert_seqs = [s for _, s in records[split_idx + 1:]]
    all_seqs   = vert_seqs + invert_seqs

    # Validate alignment lengths
    all_lengths = {len(s) for s in all_seqs}
    if len(all_lengths) != 1:
        raise ValueError(f"Sequences have unequal lengths: {sorted(all_lengths)}")
    aln_len = all_lengths.pop()

    # Find human sequence
    human_seq = find_human(records[:split_idx])  # human should be a vertebrate
    if human_seq is None:
        human_seq = find_human(records)
    if human_seq is None:
        raise ValueError("No human sequence found (looked for 'Homo_sapiens' or 'Human')")

    # Build one entry per alignment column.
    # Human numbering advances only at non-gap human residues.
    human_pos = 0               # running count of human residues seen
    col_labels = []             # alignment columns (1-indexed)
    human_positions = []        # human protein positions or None where human is gap
    v_scores = []               # vertebrate % conserved
    i_scores = []               # invertebrate % conserved
    a_scores = []               # all-species % conserved

    for col_idx in range(aln_len):
        col_labels.append(col_idx + 1)
        human_aa = human_seq[col_idx]

        if human_aa in {"-", ".", "X"}:
            human_positions.append(None)
        else:
            human_pos += 1
            human_positions.append(human_pos)

        v_scores.append(pct_conserved(vert_seqs,   col_idx))
        i_scores.append(pct_conserved(invert_seqs, col_idx))
        a_scores.append(pct_conserved(all_seqs,    col_idx))

    # Write CSV — horizontal layout: one column per alignment column
    stem = IN_FASTA.stem
    out_path = OUT_DIR / f"{stem}_conservation_perresidue.csv"

    def fmt(val):
        return str(val) if val is not None else ""

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        # Header row uses alignment-column numbers so every alignment position is represented.
        fh.write("metric," + ",".join(str(c) for c in col_labels) + "\n")
        fh.write("human_position," + ",".join(fmt(p) for p in human_positions) + "\n")
        fh.write("pct_conserved_vertebrates,"   + ",".join(fmt(v) for v in v_scores)   + "\n")
        fh.write("pct_conserved_invertebrates," + ",".join(fmt(i) for i in i_scores) + "\n")
        fh.write("pct_conserved_all,"           + ",".join(fmt(a) for a in a_scores)   + "\n")

    print(f"Wrote: {out_path}")
    print(f"  Vertebrates  : {len(vert_seqs)} sequences")
    print(f"  Invertebrates: {len(invert_seqs)} sequences")
    print(f"  Alignment length: {aln_len} columns")
    print(f"  Columns written: {len(col_labels)}")
    print(f"  Human residues mapped: {human_pos}")


if __name__ == "__main__":
    main()
