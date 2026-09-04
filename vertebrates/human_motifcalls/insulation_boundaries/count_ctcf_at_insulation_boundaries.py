#!/usr/bin/env python3
"""
Count CTCF motifs around insulation-score boundaries (BED).

Expected boundary file format (tab-delimited):
- BED3+: chrom, start, end, ...
Boundary position is defined as midpoint: (start + end) // 2.

Output:
- Summary printed to stdout
- Optional per-boundary TSV written with --out
"""

from __future__ import annotations

import argparse
import bisect
import collections
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count CTCF motifs near insulation boundaries."
    )
    parser.add_argument("motifs_bed", help="CTCF motif BED file (chrom, start, end, strand?)")
    parser.add_argument("boundaries_bed", help="Insulation boundary BED file (chrom, start, end)")
    parser.add_argument(
        "--window-bp",
        type=int,
        default=10_000,
        help="Half-window around each boundary midpoint (default: 10000)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output TSV path with per-boundary counts",
    )
    return parser.parse_args()


def load_motif_centers(motifs_bed: Path) -> dict[str, list[int]]:
    motifs_by_chrom: dict[str, list[int]] = collections.defaultdict(list)
    with motifs_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            center = (start + end) // 2
            motifs_by_chrom[chrom].append(center)

    for chrom in motifs_by_chrom:
        motifs_by_chrom[chrom].sort()
    return motifs_by_chrom


def count_in_window(motifs_by_chrom: dict[str, list[int]], chrom: str, left: int, right: int) -> int:
    motifs = motifs_by_chrom.get(chrom)
    if not motifs:
        return 0
    left_idx = bisect.bisect_left(motifs, left)
    right_idx = bisect.bisect_right(motifs, right)
    return right_idx - left_idx


def main() -> None:
    args = parse_args()
    motifs_bed = Path(args.motifs_bed)
    boundaries_bed = Path(args.boundaries_bed)

    motifs_by_chrom = load_motif_centers(motifs_bed)

    total_boundaries = 0
    boundaries_with_ctcf = 0
    count_distribution: collections.Counter[int] = collections.Counter()
    per_boundary_rows: list[tuple[str, int, int, int, int, int]] = []

    with boundaries_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            midpoint = (start + end) // 2
            window_start = midpoint - args.window_bp
            window_end = midpoint + args.window_bp

            ctcf_count = count_in_window(motifs_by_chrom, chrom, window_start, window_end)

            total_boundaries += 1
            count_distribution[ctcf_count] += 1
            if ctcf_count > 0:
                boundaries_with_ctcf += 1

            per_boundary_rows.append((chrom, start, end, midpoint, args.window_bp, ctcf_count))

    print(f"=== CTCF Coverage at Insulation Boundaries (+/-{args.window_bp}bp) ===")
    print(f"Total boundaries: {total_boundaries}")
    if total_boundaries == 0:
        print("Boundaries with >=1 CTCF: 0 (0.0%)")
        return

    pct = 100.0 * boundaries_with_ctcf / total_boundaries
    print(f"Boundaries with >=1 CTCF: {boundaries_with_ctcf} ({pct:.2f}%)")

    print("\n=== Distribution of CTCF counts per boundary ===")
    for ctcf_count in sorted(count_distribution):
        n = count_distribution[ctcf_count]
        p = 100.0 * n / total_boundaries
        print(f"{ctcf_count} CTCFs: {n} boundaries ({p:.2f}%)")

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as handle:
            handle.write("chrom\tstart\tend\tmidpoint\twindow_bp\tctcf_count\n")
            for row in per_boundary_rows:
                handle.write("\t".join(str(value) for value in row) + "\n")
        print(f"\nPer-boundary counts written to: {out_path}")


if __name__ == "__main__":
    main()
