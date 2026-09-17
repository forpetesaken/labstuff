#!/usr/bin/env python3
"""
Analyze CTCF motif pairs stratified by chromosome.
Similar to ctcfbeds.py but analyzes each chromosome separately and creates stacked bar charts.
"""
import sys
import os
import collections
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: ctcfbeds_by_chromosome.py motifs.bed MAX_DIST [LABEL] [OUT_DIR]", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
MAX_DIST = int(sys.argv[2]) if len(sys.argv) >= 3 else 500
LABEL = sys.argv[3] if len(sys.argv) >= 4 else os.path.basename(bed_file)
OUT_DIR = sys.argv[4] if len(sys.argv) >= 5 else None

MIN_DIST = 120

def natural_sort_key(chrom):
    """Extract numeric part from chromosome names for natural sorting."""
    match = re.search(r'\d+', chrom)
    if match:
        return (0, int(match.group()))
    if 'X' in chrom.upper():
        return (1, 0)
    elif 'Y' in chrom.upper():
        return (1, 1)
    elif 'M' in chrom.upper() or 'MT' in chrom.upper():
        return (1, 2)
    else:
        return (2, 0, chrom)

# Load motifs
motifs = []
with open(bed_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        center = (start + end) // 2
        motifs.append((chrom, start, end, strand, center))

motifs.sort(key=lambda x: (natural_sort_key(x[0]), x[1]))

def classify_pair(s1, s2):
    if s1 == '+' and s2 == '-':
        return "convergent ><"
    elif s1 == '-' and s2 == '+':
        return "divergent <>"
    elif s1 == '+' and s2 == '+':
        return "tandem_plus >>"
    elif s1 == '-' and s2 == '-':
        return "tandem_minus <<"
    else:
        return "other"

# Count pairs per chromosome
chrom_counts = collections.defaultdict(collections.Counter)

prev = None
for m in motifs:
    if prev is None:
        prev = m
        continue
    
    chrom, start, end, strand, center = m
    pchrom, pstart, pend, pstrand, pcenter = prev
    
    if chrom != pchrom:
        prev = m
        continue
    
    dist = center - pcenter
    
    if MIN_DIST < dist <= MAX_DIST:
        pair_type = classify_pair(pstrand, strand)
        chrom_counts[chrom][pair_type] += 1
    
    prev = m

# Print summary
print(f"\n# Pair counts by chromosome (>{MIN_DIST}bp, ≤{MAX_DIST}bp):", file=sys.stderr)
total = 0
for chrom in sorted(chrom_counts.keys(), key=natural_sort_key):
    chrom_total = sum(chrom_counts[chrom].values())
    total += chrom_total
    print(f"#   {chrom}: {chrom_total} pairs", file=sys.stderr)
print(f"#   TOTAL: {total} pairs", file=sys.stderr)

# Create stacked bar chart
if chrom_counts:
    # Get chromosomes with at least 3 pairs
    chromosomes = [c for c in sorted(chrom_counts.keys(), key=natural_sort_key) 
                   if sum(chrom_counts[c].values()) >= 3]
    
    if not chromosomes:
        print("# No chromosomes with sufficient data (≥3 pairs)", file=sys.stderr)
        sys.exit(0)
    
    # Define pair types and colors
    pair_types = ['convergent ><', 'divergent <>', 'tandem_plus >>', 'tandem_minus <<']
    colors = {
        'convergent ><': '#FF6B6B',
        'divergent <>': '#4ECDC4',
        'tandem_plus >>': '#45B7D1',
        'tandem_minus <<': '#FFA07A'
    }
    
    # Calculate percentages
    percentages = {pt: [] for pt in pair_types}
    for chrom in chromosomes:
        total_chrom = sum(chrom_counts[chrom].values())
        for pt in pair_types:
            count = chrom_counts[chrom][pt]
            pct = 100 * count / total_chrom if total_chrom > 0 else 0
            percentages[pt].append(pct)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(chromosomes) * 0.4), 6))
    
    x = np.arange(len(chromosomes))
    width = 0.6
    
    # Plot stacked bars
    bottom = np.zeros(len(chromosomes))
    for pt in pair_types:
        bars = ax.bar(x, percentages[pt], width, label=pt, 
                     color=colors[pt], alpha=0.8, bottom=bottom)
        bottom += np.array(percentages[pt])
    
    # Customize plot
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'CTCF Pair Orientations by Chromosome\n{LABEL}, window ≤{MAX_DIST}bp, >{MIN_DIST}bp (n={sum(sum(chrom_counts[c].values()) for c in chromosomes)})', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right')
    ax.set_ylim(0, 100)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save plot
    if OUT_DIR:
        os.makedirs(OUT_DIR, exist_ok=True)
        png_path = os.path.join(OUT_DIR, f"ctcf_pairs_by_chrom_{LABEL}_{MAX_DIST}bp.png")
    else:
        png_path = f"ctcf_pairs_by_chrom_{LABEL}_{MAX_DIST}bp.png"
    
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    print(f"# Saved plot: {png_path}", file=sys.stderr)
    plt.close()
else:
    print("# No pairs found", file=sys.stderr)
