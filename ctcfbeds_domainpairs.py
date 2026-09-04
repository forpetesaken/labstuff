#!/usr/bin/env python3
"""
Determine CTCF pair types using domain boundary pairs (BEDPE).
For each domain boundary pair (anchor1, anchor2), find CTCF motifs within 10kb
on the *interior side* of each boundary coordinate and classify a single pair
formed by the nearest motif to each boundary. No distance filter between motifs.

Usage:
    ctcfbeds_domainpairs.py motifs.bed domains.bedpe [WINDOW_BP] [MODE] [LABEL] [OUT_DIR]
    MODE: nearest | interior | discard | strongest (default: nearest)
"""
import sys
import os
import collections
import bisect
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    print("Usage: ctcfbeds_domainpairs.py motifs.bed domains.bedpe [WINDOW_BP] [MODE] [LABEL] [OUT_DIR]", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]
WINDOW_BP = 10000
MODE = "nearest"

arg_idx = 3
if len(sys.argv) >= 4:
    try:
        WINDOW_BP = int(sys.argv[3])
        arg_idx = 4
    except ValueError:
        WINDOW_BP = 10000
        arg_idx = 3

if len(sys.argv) >= arg_idx + 1:
    if sys.argv[arg_idx] in ("nearest", "interior", "discard", "strongest"):
        MODE = sys.argv[arg_idx]
        arg_idx += 1

LABEL = os.path.basename(bed_file)
if len(sys.argv) >= arg_idx + 1:
    LABEL = sys.argv[arg_idx]

OUT_DIR = None
if len(sys.argv) >= arg_idx + 2:
    OUT_DIR = sys.argv[arg_idx + 1]

# Load motifs (store centers for fast nearest-boundary lookup)
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

# Sort motifs by center for each chrom
for chrom in motifs_by_chrom:
    motifs_by_chrom[chrom].sort(key=lambda x: x[0])


def get_motifs_in_window(chrom, window_start, window_end):
    """Return motif tuples within [window_start, window_end]."""
    motifs = motifs_by_chrom.get(chrom)
    if not motifs:
        return []
    centers = [m[0] for m in motifs]
    left = bisect.bisect_left(centers, window_start)
    right = bisect.bisect_right(centers, window_end)
    if left >= right:
        return []
    return motifs[left:right]


def select_motif_strand_in_window(chrom, window_start, window_end, target_pos, mode):
    """Return strand of selected motif within [window_start, window_end] using mode."""
    candidates = get_motifs_in_window(chrom, window_start, window_end)
    if not candidates:
        return None
    if mode == "discard" and len(candidates) > 1:
        return None
    if mode == "interior":
        # Choose the most interior motif: farthest from boundary into the domain
        if target_pos == window_start:
            best = max(candidates, key=lambda m: m[0])  # rightmost
        else:
            best = min(candidates, key=lambda m: m[0])  # leftmost
        return best[1]
    if mode == "strongest":
        # Choose highest motif score; fall back to nearest if no scores present
        scored = [m for m in candidates if m[2] is not None]
        if scored:
            best = max(scored, key=lambda m: m[2])
            return best[1]
    # nearest (default)
    best = min(candidates, key=lambda m: abs(m[0] - target_pos))
    return best[1]


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

# Process domain pairs
with open(domains_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 6:
            continue

        chrom1 = fields[0]
        start1 = int(fields[1])
        end1 = int(fields[2])
        chrom2 = fields[3]
        start2 = int(fields[4])
        end2 = int(fields[5])

        # Skip if anchors are on different chromosomes
        if chrom1 != chrom2:
            continue
        
        chrom = chrom1

        # Use anchor regions as the boundaries
        if chrom not in motifs_by_chrom:
            continue

        # For left boundary (anchor1): search on interior side (extending rightward from anchor1.start)
        # This assumes anchor1.start is the left domain boundary
        left_boundary_pos = start1  # Left edge of left anchor
        left_strand = select_motif_strand_in_window(
            chrom,
            window_start=left_boundary_pos,
            window_end=left_boundary_pos + WINDOW_BP,
            target_pos=left_boundary_pos,
            mode=MODE
        )
        
        # For right boundary (anchor2): search on interior side (extending leftward from anchor2.end)
        # This assumes anchor2.end is the right domain boundary
        right_boundary_pos = end2  # Right edge of right anchor
        right_strand = select_motif_strand_in_window(
            chrom,
            window_start=right_boundary_pos - WINDOW_BP,
            window_end=right_boundary_pos,
            target_pos=right_boundary_pos,
            mode=MODE
        )

        if MODE == "all":
            left_candidates = get_motifs_in_window(
                chrom,
                window_start=left_boundary_pos,
                window_end=left_boundary_pos + WINDOW_BP,
            )
            right_candidates = get_motifs_in_window(
                chrom,
                window_start=right_boundary_pos - WINDOW_BP,
                window_end=right_boundary_pos,
            )
            if not left_candidates or not right_candidates:
                continue
            for l in left_candidates:
                for r in right_candidates:
                    pair_counts[classify_pair(l[1], r[1])] += 1
        else:
            if left_strand is None or right_strand is None:
                continue

            pair_counts[classify_pair(left_strand, right_strand)] += 1

# Output TSV to stdout
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

    plt.title(f"CTCF motif pairs by domain boundaries\n{LABEL} (n={total})",
              fontsize=14, weight='bold', pad=20)
    plt.tight_layout()

    if OUT_DIR:
        out_dir = OUT_DIR
    else:
        out_dir = os.path.dirname(bed_file) if os.path.dirname(bed_file) else '.'

    pie_chart_dir = os.path.join(out_dir, "pie_charts", "domainpairs")
    os.makedirs(pie_chart_dir, exist_ok=True)

    safe_label = str(LABEL).replace('/', '_').replace('\\', '_')
    out_path = os.path.join(pie_chart_dir, f"ctcf_pairs_{safe_label}_domainpairs_pie.png")

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"# Pie chart saved to {out_path}", file=sys.stderr)
