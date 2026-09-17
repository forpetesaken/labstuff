#!/usr/bin/env python3
"""
Generate % conservation table for each amino acid position
starting from human protein position 226 (the "Y").

Works dynamically with any updated alignment.
Outputs both vertebrate and invertebrate versions.
"""

from pathlib import Path
from collections import Counter
import pandas as pd
import sys

def read_fasta_ordered(path):
    """Read FASTA and preserve order."""
    items = []
    sid = None
    seq = []
    for raw in path.read_text(encoding='utf-8').splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if sid is not None:
                items.append((sid, ''.join(seq)))
            sid = line[1:].strip()
            seq = []
        else:
            if sid is None:
                raise ValueError('Invalid FASTA')
            seq.append(line.upper())
    if sid is not None:
        items.append((sid, ''.join(seq)))
    return items

def find_human_sequence(records):
    """Find human sequence (by header containing 'Homo_sapiens' or 'Human')."""
    for h, s in records:
        if 'Homo_sapiens' in h or 'Human' in h:
            return s
    # Fallback: first sequence
    return records[0][1] if records else None

def map_protein_pos_to_alignment_col(human_seq, protein_pos):
    """
    Find the alignment column corresponding to human protein position.
    protein_pos: 1-indexed protein residue number (e.g., 226 for the "Y")
    Returns: 0-indexed alignment column where that residue is located.
    """
    nongap_count = 0
    for col_idx, aa in enumerate(human_seq):
        if aa not in {'-', '.', 'X'}:
            nongap_count += 1
            if nongap_count == protein_pos:
                return col_idx
    raise ValueError(f"Protein position {protein_pos} not found in human sequence")

def truncate_to_human_range(df, start_pos, end_pos):
    """Restrict a conservation table to a human protein position range."""
    return df[(df['protein_position'] >= start_pos) & (df['protein_position'] <= end_pos)].copy()

def calculate_conservation_table(records, start_col_idx, reference_human_seq=None):
    """
    Calculate % conservation from start_col_idx onwards.
    Returns one row per human amino acid position with alignment column,
    human residue, and % conserved.
    """
    human_seq = reference_human_seq if reference_human_seq is not None else find_human_sequence(records)
    if human_seq is None:
        raise ValueError("No human sequence found")
    
    seqs = [s for _, s in records]
    if not seqs:
        raise ValueError("No sequences in group")
    
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Unequal sequence lengths: {sorted(lengths)}")
    
    rows = []
    aln_len = len(seqs[0])
    protein_pos = sum(1 for aa in human_seq[:start_col_idx] if aa not in {'-', '.', 'X'})
    
    for col_idx in range(start_col_idx, aln_len):
        if human_seq[col_idx] in {'-', '.', 'X'}:
            continue

        protein_pos += 1
        col = [s[col_idx] for s in seqs]
        nongap = [aa for aa in col if aa not in {'-', '.', 'X'}]
        
        human_aa = human_seq[col_idx] if col_idx < len(human_seq) else '-'
        
        if len(nongap) == 0:
            pct = float('nan')
        else:
            counts = Counter(nongap)
            pct = 100.0 * max(counts.values()) / len(nongap)
        
        rows.append({
            'protein_position': protein_pos,
            'alignment_column': col_idx + 1,
            'human_residue': human_aa,
            'pct_conserved': pct,
            'nongap_count': len(nongap)
        })
    
    return pd.DataFrame(rows), human_seq

def main():
    in_fasta = Path(r'C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\alignments_out\CTCF\top_hit_CTCF_Xenopus_laevis.fas')
    out_dir = Path(r'C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\alignments_out\CTCF\05')
    
    if not in_fasta.exists():
        print(f"ERROR: FASTA file not found: {in_fasta}")
        sys.exit(1)
    
    records = read_fasta_ordered(in_fasta)
    human_seq = find_human_sequence(records)
    
    # Find alignment column for human protein position 226 (the "Y")
    protein_pos_start = 226
    protein_pos_end = 401
    start_col_idx = map_protein_pos_to_alignment_col(human_seq, protein_pos_start)
    human_aa_at_226 = human_seq[start_col_idx]
    
    print(f"Human protein position {protein_pos_start}: '{human_aa_at_226}'")
    print(f"Alignment column (1-indexed): {start_col_idx + 1}")
    print()
    
    # Find INVERTEBRATES marker
    headers = [h for h, _ in records]
    if 'INVERTEBRATES' not in headers:
        print("WARNING: INVERTEBRATES marker not found; will generate full-set table only")
        invert_idx = None
    else:
        invert_idx = headers.index('INVERTEBRATES')
        vertebrates = records[:invert_idx]
        invertebrates = records[invert_idx + 1:]
    
    # Generate tables
    if invert_idx is not None:
        v_df, v_human = calculate_conservation_table(vertebrates, start_col_idx, reference_human_seq=human_seq)
        i_df, i_human = calculate_conservation_table(invertebrates, start_col_idx, reference_human_seq=human_seq)
        v_df = truncate_to_human_range(v_df, protein_pos_start, protein_pos_end)
        i_df = truncate_to_human_range(i_df, protein_pos_start, protein_pos_end)
        
        v_out = out_dir / 'conservation_table_from_human_pos226_vertebrates.csv'
        i_out = out_dir / 'conservation_table_from_human_pos226_invertebrates.csv'
        
        # Transpose: each position becomes a column
        v_transposed = v_df.set_index('protein_position').T
        i_transposed = i_df.set_index('protein_position').T
        
        v_transposed.to_csv(v_out)
        i_transposed.to_csv(i_out)
        
        print(f"Wrote: {v_out}")
        print(f"Wrote: {i_out}")
        print(f"\nVertebrates: {len(v_df)} positions from alignment column {start_col_idx + 1} onwards")
        print(f"Invertebrates: {len(i_df)} positions from alignment column {start_col_idx + 1} onwards")
        
        # Also generate full (all species ignoring splitter)
        all_records = [r for r in records if r[0] != 'INVERTEBRATES']
        full_df, _ = calculate_conservation_table(all_records, start_col_idx, reference_human_seq=human_seq)
        full_df = truncate_to_human_range(full_df, protein_pos_start, protein_pos_end)
        full_out = out_dir / 'conservation_table_from_human_pos226_all_species.csv'
        
        # Transpose: each position becomes a column, then add split-group conservation rows.
        full_transposed = full_df.set_index('protein_position').T
        full_transposed.loc['pct_conserved_vertebrates'] = v_df.set_index('protein_position')['pct_conserved']
        full_transposed.loc['pct_conserved_invertebrates'] = i_df.set_index('protein_position')['pct_conserved']
        full_transposed.to_csv(full_out)
        
        print(f"Wrote: {full_out}")
        print(
            f"All species (ignoring INVERTEBRATES splitter): {len(full_df)} positions "
            f"from human position {protein_pos_start} to {protein_pos_end}"
        )
    else:
        full_df, _ = calculate_conservation_table(records, start_col_idx, reference_human_seq=human_seq)
        full_df = truncate_to_human_range(full_df, protein_pos_start, protein_pos_end)
        full_out = out_dir / 'conservation_table_from_human_pos226_full.csv'
        
        # Transpose: each position becomes a column
        full_transposed = full_df.set_index('protein_position').T
        full_transposed.to_csv(full_out)
        
        print(f"Wrote: {full_out}")
        print(f"\nFull set: {len(full_df)} positions from alignment column {start_col_idx + 1} onwards")

if __name__ == '__main__':
    main()
