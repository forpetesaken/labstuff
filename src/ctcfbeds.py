#!/usr/bin/env python3
import sys
import os
import collections
import matplotlib.pyplot as plt
import re

if len(sys.argv) < 3:
    print("Usage: ctcfbeds.py motifs.sorted.bed MAX_DIST [LABEL] [OUT_DIR] > pairs.tsv", file=sys.stderr)
    print("  - LABEL: optional short label for file names (e.g. '4.0')", file=sys.stderr)
    print("  - OUT_DIR: optional output directory for PNGs/TSVs (absolute or relative)", file=sys.stderr)
    print("Note: Make sure you're running from the ctcf_bedmotif_analysis directory", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
MAX_DIST = int(sys.argv[2])

LABEL = os.path.basename(bed_file)
OUT_DIR = None
if len(sys.argv) >= 4:
    LABEL = sys.argv[3]
if len(sys.argv) >= 5:
    OUT_DIR = sys.argv[4]

# Column indices for the bed file
CHROM_COL = 0
START_COL = 1
END_COL   = 2
STRAND_COL = 3 

motifs = []
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
        # Motif tuple format: (chrom, start, end, strand, center)
        motifs.append((chrom, start, end, strand, center))

def natural_sort_key(chrom):
    """
    Extract numeric part from chromosome names for natural sorting.
    E.g., 'chr10' -> 10, 'chrX' -> 100, 'chrY' -> 101, 'chrM' -> 102
    """
    # Try to extract the numeric part
    match = re.search(r'\d+', chrom)
    if match:
        return (0, int(match.group()))
    # Handle special chromosomes
    if 'X' in chrom.upper():
        return (1, 0)  # X comes after numbers
    elif 'Y' in chrom.upper():
        return (1, 1)  # Y after X
    elif 'M' in chrom.upper() or 'MT' in chrom.upper():
        return (1, 2)  # M/MT after Y
    else:
        # For scaffolds, contigs, etc., sort alphabetically after standard chroms
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

    dist = center - pcenter  # motif2 - motif1 (bp)

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
    # Define consistent colors for each pair type
    color_map = {
        "convergent ><": '#B19CD9',  # lilac/lavender
        "divergent <>": '#87CEEB',   # sky blue
        "tandem_plus >>": '#98D8C8', # pearl aqua
        "tandem_minus <<": '#FFB347' # pastel orange
    }
    
    labels = []
    sizes = []
    colors = []
    for k, v in pair_counts.items():
        if v > 0:
            labels.append(k)
            sizes.append(v)
            colors.append(color_map.get(k, '#CCCCCC'))  # default gray for unknown types

    total = sum(sizes)
    
    # Print counts
    print(f"\n# Pair type counts (window ≤ {MAX_DIST} bp):", file=sys.stderr)
    for k, v in sorted(pair_counts.items(), key=lambda x: x[1], reverse=True):
        if v > 0:
            pct = 100.0 * v / total
            print(f"#   {k}: {v} ({pct:.1f}%)", file=sys.stderr)
    
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
    
    plt.title(f"CTCF motif {LABEL}, window ≤ {MAX_DIST} bp (n={total})", 
              fontsize=14, weight='bold', pad=20)
    plt.tight_layout()

    # safe-ish filename
    safe_label = str(LABEL).replace('/', '_').replace('\\', '_')

    # If the user supplied an output directory on the command line, use it.
    # Otherwise, keep the previous behavior and infer human vs ciona output locations.
    if OUT_DIR:
        out_dir = OUT_DIR
    else:
        bed_file_lower = bed_file.replace('/', '\\').lower()
        if ('hs_fimo_out' in bed_file_lower) or ('hs_with_hs_motif' in bed_file_lower) or ('\\hs\\' in bed_file_lower):
            out_dir = r"C:\Users\Nat\Downloads\AIDEN Lab\Code\ctcf_bedmotif_analysis\human_motifcalls"
        else:
            out_dir = r"C:\Users\Nat\Downloads\AIDEN Lab\Code\ctcf_bedmotif_analysis\cion_intestines_motifcalls"

    # Create subdirectories for organized output
    pie_chart_dir = os.path.join(out_dir, "pie_charts", "ctcfbeds")
    tsv_dir = os.path.join(out_dir, "tsvs")
    os.makedirs(pie_chart_dir, exist_ok=True)
    os.makedirs(tsv_dir, exist_ok=True)

    out_png = f"ctcf_pairs_{safe_label}_{MAX_DIST}bp_pie.png"
    out_path = os.path.join(pie_chart_dir, out_png)
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"# Pie chart saved to {out_path}", file=sys.stderr)
    print(f"# TSV output should be redirected to {os.path.join(tsv_dir, f'ctcf_pairs_{safe_label}_{MAX_DIST}bp.tsv')}", file=sys.stderr)
