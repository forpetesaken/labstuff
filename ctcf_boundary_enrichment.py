#!/usr/bin/env python3
"""
Quantify CTCF motif enrichment at domain boundaries relative to shift controls.
For each boundary, count motifs within WINDOW_BP. Compare to shifted controls.
"""
import sys
import os
import random
import numpy as np

if len(sys.argv) < 3:
    print("Usage: ctcf_boundary_enrichment.py motifs.bed domains.bedpe [WINDOW_BP] [NUM_SHIFTS] [SHIFT_DIST] [MODE]", file=sys.stderr)
    print("  WINDOW_BP: distance window around boundary (default: 10000)", file=sys.stderr)
    print("  NUM_SHIFTS: number of random shift controls (default: 100)", file=sys.stderr)
    print("  SHIFT_DIST: max shift distance in bp (default: 500000)", file=sys.stderr)
    print("  MODE: 'global' (same shift per iteration) or 'independent' (each boundary shifts independently) (default: global)", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]
WINDOW_BP = int(sys.argv[3]) if len(sys.argv) >= 4 else 10000
NUM_SHIFTS = int(sys.argv[4]) if len(sys.argv) >= 5 else 100
SHIFT_DIST = int(sys.argv[5]) if len(sys.argv) >= 6 else 500000
MODE = sys.argv[6].lower() if len(sys.argv) >= 7 else 'global'

# Load motifs by chromosome
from collections import defaultdict
motifs_by_chrom = defaultdict(list)
with open(bed_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3:
            continue
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        center = (start + end) // 2
        motifs_by_chrom[chrom].append(center)

for chrom in motifs_by_chrom:
    motifs_by_chrom[chrom].sort()

# Load boundaries (use start and end of each domain interval)
boundaries = []
with open(domains_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3:
            continue
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        boundaries.append((chrom, start))
        boundaries.append((chrom, end))

print(f"Loaded {len(boundaries)} domain boundaries", file=sys.stderr)

def count_motifs_near(chrom, pos, window):
    """Count motifs within [pos-window, pos+window]"""
    if chrom not in motifs_by_chrom:
        return 0
    motifs = motifs_by_chrom[chrom]
    count = 0
    for m in motifs:
        if abs(m - pos) <= window:
            count += 1
    return count

# Count at actual boundaries
actual_counts = []
for chrom, pos in boundaries:
    count = count_motifs_near(chrom, pos, WINDOW_BP)
    actual_counts.append(count)

actual_mean = np.mean(actual_counts)
actual_total = sum(actual_counts)

print(f"\nActual boundaries:")
print(f"  Total motifs within {WINDOW_BP}bp: {actual_total}")
print(f"  Mean motifs per boundary: {actual_mean:.2f}")

# Generate shift controls
print(f"\nRunning {NUM_SHIFTS} shift controls (mode: {MODE})...", file=sys.stderr)
shift_means = []

if MODE == 'independent':
    # Each boundary gets independent random shift
    for shift_idx in range(NUM_SHIFTS):
        shift_counts = []
        for chrom, pos in boundaries:
            shift_offset = random.randint(-SHIFT_DIST, SHIFT_DIST)
            shifted_pos = pos + shift_offset
            if shifted_pos < 0:
                shifted_pos = 0
            count = count_motifs_near(chrom, shifted_pos, WINDOW_BP)
            shift_counts.append(count)
        shift_means.append(np.mean(shift_counts))
else:
    # Global shift: all boundaries shifted by same amount per iteration
    for shift_idx in range(NUM_SHIFTS):
        shift_offset = random.randint(-SHIFT_DIST, SHIFT_DIST)
        shift_counts = []
        for chrom, pos in boundaries:
            shifted_pos = pos + shift_offset
            if shifted_pos < 0:
                shifted_pos = 0
            count = count_motifs_near(chrom, shifted_pos, WINDOW_BP)
            shift_counts.append(count)
        shift_means.append(np.mean(shift_counts))

control_mean = np.mean(shift_means)
control_std = np.std(shift_means)

enrichment = actual_mean / control_mean if control_mean > 0 else float('inf')
z_score = (actual_mean - control_mean) / control_std if control_std > 0 else float('inf')

print(f"\nShift controls ({NUM_SHIFTS} iterations):")
print(f"  Mean motifs per boundary: {control_mean:.2f} ± {control_std:.2f}")
print(f"\nEnrichment: {enrichment:.2f}x")
print(f"Z-score: {z_score:.2f}")
