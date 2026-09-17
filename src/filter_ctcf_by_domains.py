#!/usr/bin/env python3
"""
Filter CTCF motifs to only include those within 10kb of domain boundaries.
Uses domain calls from 5000_blocks.bedpe file.
"""
import sys
import os
import collections

if len(sys.argv) < 4:
    print("Usage: filter_ctcf_by_domains.py motifs.bed domains.bedpe DISTANCE [LABEL] [OUT_DIR]", file=sys.stderr)
    print("  - motifs.bed: CTCF motif BED file (chrom, start, end, strand)", file=sys.stderr)
    print("  - domains.bedpe: Domain calls BEDPE file (e.g., 5000_blocks.bedpe)", file=sys.stderr)
    print("  - DISTANCE: Maximum distance from domain boundaries in bp (e.g., 10000)", file=sys.stderr)
    print("  - LABEL: optional short label for file names (e.g. '4.0')", file=sys.stderr)
    print("  - OUT_DIR: optional output directory for filtered BED file", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]
DISTANCE = int(sys.argv[3])

LABEL = os.path.basename(bed_file).replace('.bed', '')
if len(sys.argv) >= 5:
    LABEL = sys.argv[4]

OUT_DIR = None
if len(sys.argv) >= 6:
    OUT_DIR = sys.argv[5]

print(f"Loading domain boundaries from {domains_file}...", file=sys.stderr)

# Load domain boundaries - we'll use both the anchor1 and anchor2 regions
domain_regions = []
with open(domains_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 6:
            continue
        
        # Extract anchor1 coordinates
        chrom1 = fields[0]
        start1 = int(fields[1])
        end1 = int(fields[2])
        domain_regions.append((chrom1, start1, end1))
        
        # Extract anchor2 coordinates if they exist
        if len(fields) >= 7:
            chrom2 = fields[3]
            start2 = int(fields[4])
            end2 = int(fields[5])
            domain_regions.append((chrom2, start2, end2))

print(f"Loaded {len(domain_regions)} domain boundary regions", file=sys.stderr)

# Build an index for fast lookup: chrom -> sorted list of (start, end)
domain_index = collections.defaultdict(list)
for chrom, start, end in domain_regions:
    domain_index[chrom].append((start, end))

# Sort each chromosome's domains by start position
for chrom in domain_index:
    domain_index[chrom].sort()

def near_domain_boundary(chrom, motif_start, motif_end):
    """
    Check if motif is within DISTANCE bp of any domain boundary.
    Returns True if motif overlaps or is within DISTANCE of any domain region.
    """
    if chrom not in domain_index:
        return False
    
    domains = domain_index[chrom]
    
    # Check each domain region
    for domain_start, domain_end in domains:
        # Check if motif is within DISTANCE of domain start or end
        # Motif is close to domain if:
        # 1. It overlaps the domain
        # 2. It's within DISTANCE bp before the domain starts
        # 3. It's within DISTANCE bp after the domain ends
        
        if not (motif_end + DISTANCE < domain_start or motif_start > domain_end + DISTANCE):
            return True
    
    return False

print(f"Loading CTCF motifs from {bed_file}...", file=sys.stderr)

motifs = []
total_motifs = 0
filtered_motifs = 0

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
        strand = fields[3] if len(fields) > 3 else '.'
        
        total_motifs += 1
        
        # Check if motif is near a domain boundary
        if near_domain_boundary(chrom, start, end):
            filtered_motifs += 1
            motifs.append(line.rstrip('\n'))

print(f"Total motifs: {total_motifs}", file=sys.stderr)
print(f"Motifs near domain boundaries (within {DISTANCE} bp): {filtered_motifs}", file=sys.stderr)
print(f"Percentage retained: {100.0 * filtered_motifs / total_motifs:.2f}%", file=sys.stderr)

# Output filtered motifs
if OUT_DIR:
    output_file = os.path.join(OUT_DIR, f"{LABEL}_near_domains_{DISTANCE//1000}kb.bed")
else:
    output_file = f"{LABEL}_near_domains_{DISTANCE//1000}kb.bed"

os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)

with open(output_file, 'w') as f:
    for motif in motifs:
        f.write(motif + '\n')

print(f"Filtered motifs written to {output_file}", file=sys.stderr)
