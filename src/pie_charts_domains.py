#!/usr/bin/env python3
"""
Generate pie charts for domain-filtered CTCF motifs
Shows strand pair orientation distributions
"""
import sys
import os
import collections
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import re

if len(sys.argv) < 3:
    print("Usage: pie_charts_domains.py motifs.bed MAX_DIST [LABEL] [OUT_DIR]", file=sys.stderr)
    print("  - motifs.bed: CTCF motif BED file (chrom, start, end, strand)", file=sys.stderr)
    print("  - MAX_DIST: Maximum distance between consecutive motifs (bp)", file=sys.stderr)
    print("  - LABEL: optional short label for file names (e.g. '5.5_domains')", file=sys.stderr)
    print("  - OUT_DIR: optional output directory for PNG files", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
MAX_DIST = int(sys.argv[2])

LABEL = os.path.basename(bed_file).replace('.bed', '')
if len(sys.argv) >= 4:
    LABEL = sys.argv[3]

OUT_DIR = None
if len(sys.argv) >= 5:
    OUT_DIR = sys.argv[4]

# Column indices for the bed file
CHROM_COL = 0
START_COL = 1
END_COL   = 2
STRAND_COL = 3

print(f"Loading CTCF motifs from {bed_file}...", file=sys.stderr)

motifs = []
with open(bed_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4:
            continue
        chrom = fields[CHROM_COL]
        start = int(fields[START_COL])
        end   = int(fields[END_COL])
        strand = fields[STRAND_COL]
        center = (start + end) // 2
        # Motif tuple format: (chrom, start, end, strand, center)
        motifs.append((chrom, start, end, strand, center))

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

pair_counts = collections.Counter()

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

    if dist <= MAX_DIST and dist > 120:
        pair_type = classify_pair(pstrand, strand)
        pair_counts[pair_type] += 1

    prev = m

# Create pie chart
if pair_counts:
    color_map = {
        "convergent ><": '#00BFFF',  # bright sky blue
        "divergent <>": '#CCCCCC',   # light grey
        "tandem_plus >>": '#999999', # medium grey
        "tandem_minus <<": '#666666' # dark grey
    }
    
    labels = []
    sizes = []
    colors = []
    for k, v in pair_counts.items():
        if v > 0:
            labels.append(k)
            sizes.append(v)
            colors.append(color_map.get(k, '#CCCCCC'))

    total = sum(sizes)
    
    print(f"\nPair type counts (window ≤ {MAX_DIST} bp):", file=sys.stderr)
    for k, v in sorted(pair_counts.items(), key=lambda x: x[1], reverse=True):
        if v > 0:
            pct = 100.0 * v / total
            print(f"  {k}: {v} ({pct:.1f}%)", file=sys.stderr)
    
    plt.figure(figsize=(10, 8))
    wedges, texts, autotexts = plt.pie(
        sizes, 
        labels=labels, 
        autopct='%1.1f%%',
        startangle=90,
        colors=colors,
        textprops={'fontsize': 11, 'weight': 'bold'}
    )
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(10)
        autotext.set_weight('bold')
    
    plt.title(f"CTCF motif pairs {LABEL}\nwindow ≤ {MAX_DIST} bp (n={total})", 
              fontsize=14, weight='bold', pad=20)
    plt.tight_layout()

    # Create output directory if needed
    if OUT_DIR:
        out_dir = OUT_DIR
    else:
        out_dir = os.path.dirname(bed_file) if os.path.dirname(bed_file) else '.'

    pie_chart_dir = os.path.join(out_dir, "pie_charts", "domains_filter")
    os.makedirs(pie_chart_dir, exist_ok=True)

    safe_label = str(LABEL).replace('/', '_').replace('\\', '_')
    out_path = os.path.join(pie_chart_dir, f"{safe_label}_pairs.png")
    
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Pie chart saved to {out_path}", file=sys.stderr)
    plt.close()
else:
    print("No pairs found within distance threshold", file=sys.stderr)
