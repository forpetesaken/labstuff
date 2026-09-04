#!/usr/bin/env python3
"""
Count how many CTCF motifs are found at domain boundaries.
Reports statistics on boundary coverage.
"""
import sys
import bisect
import collections

if len(sys.argv) < 3:
    print("Usage: python count_ctcf_at_boundaries.py <bed_file> <domains_bedpe> [window_bp]")
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]
WINDOW_BP = 10000
if len(sys.argv) >= 4:
    WINDOW_BP = int(sys.argv[3])

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
        center = (start + end) // 2
        motifs_by_chrom[chrom].append(center)

# Sort motifs by center for each chrom
for chrom in motifs_by_chrom:
    motifs_by_chrom[chrom].sort()


def count_motifs_in_window(chrom, window_start, window_end):
    """Count motifs within [window_start, window_end]."""
    motifs = motifs_by_chrom.get(chrom)
    if not motifs:
        return 0
    left = bisect.bisect_left(motifs, window_start)
    right = bisect.bisect_right(motifs, window_end)
    return right - left


# Process domains and count boundary overlaps
total_domains = 0
left_with_ctcf = 0
right_with_ctcf = 0
both_with_ctcf = 0
left_counts = collections.Counter()
right_counts = collections.Counter()

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
        
        total_domains += 1
        
        # Count CTCFs at left boundary interior window (anchor1.start -> anchor1.start + WINDOW_BP)
        left_boundary_pos = start1
        left_count = count_motifs_in_window(chrom1, left_boundary_pos, left_boundary_pos + WINDOW_BP)
        left_counts[left_count] += 1
        if left_count > 0:
            left_with_ctcf += 1
        
        # Count CTCFs at right boundary interior window (anchor2.end - WINDOW_BP -> anchor2.end)
        right_boundary_pos = end2
        right_count = count_motifs_in_window(chrom2, right_boundary_pos - WINDOW_BP, right_boundary_pos)
        right_counts[right_count] += 1
        if right_count > 0:
            right_with_ctcf += 1
        
        # Check if BOTH boundaries have at least 1 CTCF
        if left_count > 0 and right_count > 0:
            both_with_ctcf += 1

print(f"=== CTCF Motif Coverage at Domain Boundaries (window={WINDOW_BP}bp) ===")
print(f"\nTotal domains: {total_domains}")
print(f"\nLeft boundaries with ≥1 CTCF: {left_with_ctcf} ({100*left_with_ctcf/total_domains:.1f}%)")
print(f"Right boundaries with ≥1 CTCF: {right_with_ctcf} ({100*right_with_ctcf/total_domains:.1f}%)")
print(f"Both boundaries with ≥1 CTCF: {both_with_ctcf} ({100*both_with_ctcf/total_domains:.1f}%)")

print(f"\n=== Distribution of CTCF counts at LEFT boundaries ===")
for count in sorted(left_counts.keys()):
    num_boundaries = left_counts[count]
    print(f"{count} CTCFs: {num_boundaries} boundaries ({100*num_boundaries/total_domains:.1f}%)")

print(f"\n=== Distribution of CTCF counts at RIGHT boundaries ===")
for count in sorted(right_counts.keys()):
    num_boundaries = right_counts[count]
    print(f"{count} CTCFs: {num_boundaries} boundaries ({100*num_boundaries/total_domains:.1f}%)")

# Calculate domains eligible for different selection modes
exactly_one_both = sum(1 for i in range(total_domains) if left_counts.get(1, 0) > 0 and right_counts.get(1, 0) > 0)
print(f"\n=== Expected sample sizes for selection modes ===")
print(f"Domains with BOTH boundaries having ≥1 CTCF: {both_with_ctcf}")
print(f"  (eligible for nearest/interior/strongest modes)")
print(f"\nDomains with BOTH boundaries having exactly 1 CTCF: estimated")
print(f"  (eligible for discard mode)")
