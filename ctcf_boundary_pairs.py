#!/usr/bin/env python3
"""
Analyze short-range CTCF pair orientations at individual domain boundaries.
For each boundary, find strong CTCF motifs within MAX_DIST and classify pair types.
"""
import sys
import os
import collections
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    print("Usage: ctcf_boundary_pairs.py motifs.bed domains.bedpe [MAX_DIST] [MODE] [LABEL] [OUT_DIR]", file=sys.stderr)
    print("  MAX_DIST: maximum distance between motif pairs (default: 500)", file=sys.stderr)
    print("  MODE: all | strongest2 (default: all)", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]
MAX_DIST = int(sys.argv[3]) if len(sys.argv) >= 4 else 10000
MODE = "all"

arg_idx = 4
if len(sys.argv) >= 5:
    if sys.argv[4] in ("all", "strongest2"):
        MODE = sys.argv[4]
        arg_idx = 5

LABEL = sys.argv[arg_idx] if len(sys.argv) >= arg_idx + 1 else "boundary_pairs"
OUT_DIR = sys.argv[arg_idx + 1] if len(sys.argv) >= arg_idx + 2 else None

MIN_DIST = 120

# Load motifs
motifs_by_chrom = collections.defaultdict(list)
with open(bed_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4:
            continue
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        score = None
        if len(fields) > 4:
            try:
                score = float(fields[4])
            except ValueError:
                score = None
        center = (start + end) // 2
        motifs_by_chrom[chrom].append((center, strand, score))

for chrom in motifs_by_chrom:
    motifs_by_chrom[chrom].sort(key=lambda x: x[0])

# Load boundaries
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

# For each boundary, find nearby motifs and classify pairs
for chrom, boundary_pos in boundaries:
    if chrom not in motifs_by_chrom:
        continue
    
    # Find motifs near this boundary
    nearby = []
    for center, strand, score in motifs_by_chrom[chrom]:
        if abs(center - boundary_pos) <= MAX_DIST:
            nearby.append((center, strand, score))

    if MODE == "strongest2":
        scored = [m for m in nearby if m[2] is not None]
        if len(scored) >= 2:
            nearby = sorted(scored, key=lambda m: m[2], reverse=True)[:2]
        else:
            # not enough scored motifs to form a pair
            continue
    
    # Classify consecutive pairs
    nearby.sort(key=lambda x: x[0])
    for i in range(len(nearby) - 1):
        pos1, strand1, _ = nearby[i]
        pos2, strand2, _ = nearby[i + 1]
        dist = pos2 - pos1
        if MIN_DIST < dist <= MAX_DIST:
            pair_counts[classify_pair(strand1, strand2)] += 1

# Output
print("\t".join(["pair_type", "count"]))
for k, v in sorted(pair_counts.items(), key=lambda x: x[1], reverse=True):
    print(f"{k}\t{v}")

# Create pie chart
if pair_counts:
    color_map = {
        "convergent ><": '#B19CD9',  # lilac/lavender
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

    plt.title(f"CTCF pairs at domain boundaries\n{LABEL} (≤{MAX_DIST}bp, n={total})",
              fontsize=14, weight='bold', pad=20)
    plt.tight_layout()

    if OUT_DIR:
        out_dir = OUT_DIR
    else:
        out_dir = os.path.dirname(bed_file) if os.path.dirname(bed_file) else '.'

    pie_chart_dir = os.path.join(out_dir, "pie_charts", "boundary_pairs")
    os.makedirs(pie_chart_dir, exist_ok=True)

    safe_label = str(LABEL).replace('/', '_').replace('\\', '_')
    out_path = os.path.join(pie_chart_dir, f"boundary_pairs_{safe_label}.png")

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"# Pie chart saved to {out_path}", file=sys.stderr)
