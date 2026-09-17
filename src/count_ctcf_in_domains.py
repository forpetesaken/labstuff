#!/usr/bin/env python3
"""
Count how many domain boundaries have CTCF motifs within them.
"""
import sys
import os

if len(sys.argv) < 3:
    print("Usage: count_ctcf_in_domains.py motifs.bed domains.bedpe", file=sys.stderr)
    print("  - motifs.bed: CTCF motif BED file (chrom, start, end, strand)", file=sys.stderr)
    print("  - domains.bedpe: Domain calls BEDPE file", file=sys.stderr)
    sys.exit(1)

bed_file = sys.argv[1]
domains_file = sys.argv[2]

print(f"Loading domains from {domains_file}...", file=sys.stderr)

# Load domain boundaries
domain_regions = []
with open(domains_file) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 6:
            continue
        
        # Extract anchor1
        chrom1 = fields[0]
        start1 = int(fields[1])
        end1 = int(fields[2])
        domain_regions.append((chrom1, start1, end1))
        
        # Extract anchor2
        if len(fields) >= 7:
            chrom2 = fields[3]
            start2 = int(fields[4])
            end2 = int(fields[5])
            domain_regions.append((chrom2, start2, end2))

print(f"Loaded {len(domain_regions)} domain boundary regions", file=sys.stderr)

print(f"\nLoading CTCF motifs from {bed_file}...", file=sys.stderr)

motifs = []
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
        motifs.append((chrom, start, end))

print(f"Loaded {len(motifs)} CTCF motifs", file=sys.stderr)

# Count domains with at least one CTCF motif
domains_with_ctcf = 0

for d_chrom, d_start, d_end in domain_regions:
    has_ctcf = False
    for m_chrom, m_start, m_end in motifs:
        if m_chrom == d_chrom:
            # Check if motif overlaps domain
            if not (m_end < d_start or m_start > d_end):
                has_ctcf = True
                break
    
    if has_ctcf:
        domains_with_ctcf += 1

print(f"\nResults:")
print(f"Total domain boundaries: {len(domain_regions)}")
print(f"Domain boundaries with ≥1 CTCF motif: {domains_with_ctcf}")
print(f"Percentage: {100.0 * domains_with_ctcf / len(domain_regions):.2f}%")
