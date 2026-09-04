#!/usr/bin/env python3
"""
Filter CTCF motifs to only include those overlapping loop anchors,
then analyze motif pairs.
"""
import sys
import os
import collections
import matplotlib.pyplot as plt
import re

if len(sys.argv) < 4:
    print("Usage: ctcfbeds_loopanchors.py motifs.sorted.bed loop_anchors.bed MAX_DIST [LABEL] [OUT_DIR]", file=sys.stderr)
    print("  - motifs.sorted.bed: CTCF motif BED file (chrom, start, end, strand)", file=sys.stderr)
    print("  - loop_anchors.bed: Loop anchor BED file", file=sys.stderr)
    print("  - MAX_DIST: Maximum distance between consecutive motifs (bp)", file=sys.stderr)
    print("  - LABEL: optional short label for file names (e.g. '4.0')", file=sys.stderr)
    print("  - OUT_DIR: optional output directory for PNGs/TSVs (absolute or relative)", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
loop_anchor_file = sys.argv[2]
MAX_DIST = int(sys.argv[3])

LABEL = os.path.basename(bed_file)
OUT_DIR = None
if len(sys.argv) >= 5:
    LABEL = sys.argv[4]
if len(sys.argv) >= 6:
    OUT_DIR = sys.argv[5]

# Column indices for the bed file
CHROM_COL = 0
START_COL = 1
END_COL   = 2
STRAND_COL = 3 

print("Loading loop anchors...", file=sys.stderr)
# Load loop anchors - parse both anchor1 and anchor2 from the 7-column format
loop_anchors = []
with open(loop_anchor_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        # Anchor 1
        chrom1 = fields[0]
        start1 = int(fields[1])
        end1 = int(fields[2])
        loop_anchors.append((chrom1, start1, end1))
        # Anchor 2
        if len(fields) >= 7:
            chrom2 = fields[4]
            start2 = int(fields[5])
            end2 = int(fields[6])
            loop_anchors.append((chrom2, start2, end2))

print(f"Loaded {len(loop_anchors)} loop anchor regions", file=sys.stderr)

# Build an index for fast lookup: chrom -> list of (start, end)
loop_index = collections.defaultdict(list)
for chrom, start, end in loop_anchors:
    loop_index[chrom].append((start, end))

# Sort each chromosome's anchors by start position for efficient searching
for chrom in loop_index:
    loop_index[chrom].sort()

def overlaps_loop_anchor(chrom, motif_start, motif_end):
    """Check if motif overlaps any loop anchor on the same chromosome."""
    if chrom not in loop_index:
        return False
    
    # Use binary search to find the first anchor that could overlap
    anchors = loop_index[chrom]
    left, right = 0, len(anchors)
    
    # Binary search for first anchor where anchor_end >= motif_start
    while left < right:
        mid = (left + right) // 2
        if anchors[mid][1] < motif_start:  # anchor_end < motif_start
            left = mid + 1
        else:
            right = mid
    
    # Check a few anchors starting from 'left'
    for i in range(left, min(left + 10, len(anchors))):
        anchor_start, anchor_end = anchors[i]
        # If anchor starts after motif ends, no more overlaps possible
        if anchor_start > motif_end:
            break
        # Check for overlap
        if not (motif_end < anchor_start or motif_start > anchor_end):
            return True
    return False

print("Loading CTCF motifs...", file=sys.stderr)
motifs = []
total_motifs = 0
filtered_motifs = 0

with open(bed_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        chrom = fields[CHROM_COL]
        start = int(fields[START_COL])
        end   = int(fields[END_COL])
        strand = fields[STRAND_COL]
        center = (start + end) // 2
        
        total_motifs += 1
        
        # Filter: only keep motifs that overlap loop anchors
        if overlaps_loop_anchor(chrom, start, end):
            motifs.append((chrom, start, end, strand, center))
            filtered_motifs += 1

print(f"Loaded {total_motifs} total motifs", file=sys.stderr)
print(f"Filtered to {filtered_motifs} motifs at loop anchors ({100.0*filtered_motifs/total_motifs:.1f}%)", file=sys.stderr)

def natural_sort_key(chrom):
    """
    Extract numeric part from chromosome names for natural sorting.
    E.g., 'chr10' -> 10, 'chrX' -> 100, 'chrY' -> 101, 'chrM' -> 102
    """
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

print("\t".join([
    "chrom",
    "motif1_start", "motif1_end", "motif1_strand",
    "motif2_start", "motif2_end", "motif2_strand",
    "distance_bp",
    "pair_type"
]))

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

        print("\t".join(map(str, [
            chrom,
            pstart, pend, pstrand,
            start, end, strand,
            round(dist, 1),
            pair_type
        ])))

    prev = m

if pair_counts:
    # Print counts
    print(f"\n# Pair type counts (window ≤ {MAX_DIST} bp, at loop anchors):", file=sys.stderr)
    for k, v in sorted(pair_counts.items(), key=lambda x: x[1], reverse=True):
        if v > 0:
            pct = 100.0 * v / sum(pair_counts.values())
            print(f"#   {k}: {v} ({pct:.1f}%)", file=sys.stderr)
    
    # Create pie chart
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
    
    total = sum(sizes)
    plt.title(f"CTCF motifs at loop anchors, {LABEL}, window ≤ {MAX_DIST} bp (n={total})",
              fontsize=14, weight='bold', pad=20)
    plt.axis('equal')
    
    # Save plot
    if OUT_DIR:
        os.makedirs(OUT_DIR, exist_ok=True)
        pie_charts_dir = os.path.join(OUT_DIR, "pie_charts", "loopanchors")
        os.makedirs(pie_charts_dir, exist_ok=True)
        png_path = os.path.join(pie_charts_dir, f"ctcf_loopanchors_{LABEL}_{MAX_DIST}bp.png")
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        print(f"# Saved plot: {png_path}", file=sys.stderr)
    else:
        png_path = f"ctcf_loopanchors_{LABEL}_{MAX_DIST}bp.png"
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        print(f"# Saved plot: {png_path}", file=sys.stderr)
    
    plt.close()
else:
    print("# No pairs found within the specified distance.", file=sys.stderr)
