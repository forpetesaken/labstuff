#!/usr/bin/env python3
"""
Conservation scoring per residue for updated_alignment_0611.fas
Splits sequences into vertebrate/invertebrate groups at Ciona (first invertebrate)
"""

import pandas as pd
from pathlib import Path
from collections import Counter

INPUT_FILE = Path("C:/Users/Nat/Downloads/updated_alignment_0611.fas")
OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "conservation_scores.csv"

def read_fasta(filepath):
    """Read FASTA file, return dict of {header: sequence}"""
    fasta = {}
    current_header = None
    current_seq = []
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    fasta[current_header] = ''.join(current_seq)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            fasta[current_header] = ''.join(current_seq)
    
    return fasta

def split_vert_invert(fasta_dict):
    """
    Split sequences into vertebrate and invertebrate groups.
    Vertebrate: all sequences before Ciona_intestinalis
    Invertebrate: from Ciona onwards
    """
    headers = list(fasta_dict.keys())
    split_idx = 0
    
    for i, header in enumerate(headers):
        if 'Ciona' in header and 'intestinalis' in header:
            split_idx = i
            break
    
    vert_headers = headers[:split_idx]
    invert_headers = headers[split_idx:]
    
    print(f"Vertebrate sequences ({len(vert_headers)}): {vert_headers[0]} ... {vert_headers[-1]}")
    print(f"Invertebrate sequences ({len(invert_headers)}): {invert_headers[0]} ... {invert_headers[-1]}")
    
    return {h: fasta_dict[h] for h in vert_headers}, {h: fasta_dict[h] for h in invert_headers}

def class_conservation_gap_aware(col_residues, class_size):
    """
    Calculate conservation % for a column in a class.
    Returns percentage or None if <2 non-gap residues present.
    
    col_residues: list of residues from sequences in this class (may contain '-' for gaps)
    class_size: total number of sequences in this class (all species in group)
    """
    non_gap = [r for r in col_residues if r != '-']
    
    # Require >= 2 non-gap residues, else return None
    if len(non_gap) < 2:
        return None
    
    # Count most common residue
    counter = Counter(non_gap)
    most_common_count = counter.most_common(1)[0][1]
    
    # Denominator = all sequences in class (treating missing as 0)
    conservation_pct = round(100.0 * most_common_count / class_size)
    
    return conservation_pct

def calculate_conservation(fasta_dict, vert_dict, invert_dict):
    """
    Calculate conservation at each alignment column.
    Returns DataFrame with columns: alignment_column, vertebrate%, invertebrate%, all%
    """
    # Align all sequences in order
    headers = list(fasta_dict.keys())
    seq_length = len(next(iter(fasta_dict.values())))
    
    results = []
    
    for col_num in range(seq_length):
        # Extract residues at this column from all sequences
        col_residues_all = [fasta_dict[h][col_num] for h in headers]
        
        # Extract residues for each class
        vert_headers = list(vert_dict.keys())
        invert_headers = list(invert_dict.keys())
        
        col_residues_vert = [fasta_dict[h][col_num] for h in vert_headers]
        col_residues_invert = [fasta_dict[h][col_num] for h in invert_headers]
        
        # Calculate conservations
        cons_vert = class_conservation_gap_aware(col_residues_vert, len(vert_headers))
        cons_invert = class_conservation_gap_aware(col_residues_invert, len(invert_headers))
        cons_all = class_conservation_gap_aware(col_residues_all, len(headers))
        
        results.append({
            'alignment_column': col_num + 1,  # 1-indexed
            'vertebrate_percent': cons_vert,
            'invertebrate_percent': cons_invert,
            'all_percent': cons_all
        })
    
    return pd.DataFrame(results)

def main():
    print(f"Reading FASTA from {INPUT_FILE}...")
    fasta = read_fasta(INPUT_FILE)
    print(f"Loaded {len(fasta)} sequences, sequence length: {len(next(iter(fasta.values())))}")
    
    print("\nSplitting vertebrate/invertebrate...")
    vert_dict, invert_dict = split_vert_invert(fasta)
    
    print("\nCalculating conservation per residue...")
    df = calculate_conservation(fasta, vert_dict, invert_dict)
    
    print(f"\nWriting results to {OUTPUT_FILE}...")
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Done! Conservation scores written to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
