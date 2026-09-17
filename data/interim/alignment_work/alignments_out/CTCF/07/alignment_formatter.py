#!/usr/bin/env python3
"""
Alignment formatter for updated_alignment_0611.fas
Creates formatted Excel with conditional formatting, domain labels, conservation metrics
"""

import pandas as pd
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Border, Side, Alignment
from openpyxl.utils import get_column_letter

INPUT_FASTA = Path("C:/Users/Nat/Downloads/updated_alignment_0611.fas")
CONSERVATION_CSV = Path(__file__).parent / "conservation_scores.csv"
OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "alignment_formatted.xlsx"

# Domain feature definitions (in human CTCF coordinates)
DOMAIN_RANGES = {
    'YXF': (226, 228),
    'ZF1': (266, 288),
    'ZF2': (294, 316),
    'ZF3': (322, 344),
    'ZF4': (350, 372),
    'ZF5': (378, 400),
    'ZF6': (406, 428),
    'ZF7': (434, 456),
    'ZF8': (462, 484),
    'ZF9': (490, 512),
    'ZF10': (518, 540)
}

# Bold positions (all zinc finger residues + YXF)
BOLD_POSITIONS = set()
for start, end in DOMAIN_RANGES.values():
    for pos in range(start, end + 1):
        BOLD_POSITIONS.add(pos)

# Amino acid colors for conditional formatting
AA_COLORS = {
    'A': 'C8E6C9',  # Green
    'R': 'FFCCBC',  # Orange
    'N': 'B3E5FC',  # Cyan
    'D': 'FFCCBC',  # Orange
    'C': 'F8BBD0',  # Pink
    'Q': 'C5E1A5',  # Light Green
    'E': 'FFCCBC',  # Orange
    'G': 'C8E6C9',  # Green
    'H': 'CCCCFF',  # Light Blue
    'I': 'FFE0B2',  # Peach
    'L': 'FFE0B2',  # Peach
    'K': 'FFCCBC',  # Orange
    'M': 'FFE0B2',  # Peach
    'F': 'F0F4C3',  # Light Yellow
    'P': 'D1C4E9',  # Light Purple
    'S': 'B2DFDB',  # Teal
    'T': 'B2DFDB',  # Teal
    'W': 'F0F4C3',  # Light Yellow
    'Y': 'F0F4C3',  # Light Yellow
    'V': 'FFE0B2',  # Peach
    '-': 'EEEEEE',  # Light Gray
}

def read_fasta(filepath):
    """Read FASTA file"""
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
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            fasta[current_header] = ''.join(current_seq)
    
    return fasta

def get_human_mapping(human_seq):
    """Map non-gap positions in human to alignment columns"""
    mapping = {}
    human_pos = 0
    
    for aln_col, residue in enumerate(human_seq):
        if residue != '-':
            human_pos += 1
            mapping[human_pos] = aln_col
    
    return mapping

def color_for_conservation(pct):
    """Get fill color for conservation percentage"""
    if pct is None or pd.isna(pct):
        return 'FFFFFF'
    
    pct = float(pct)
    if pct >= 90:
        return '00B050'  # Dark Green
    elif pct >= 70:
        return 'C6EFCE'  # Light Green
    elif pct >= 50:
        return 'FFEB9C'  # Yellow
    else:
        return 'FFFFFF'  # White

def main():
    print("=== Alignment Formatter ===\n")
    
    # Read data
    print("Reading FASTA...")
    fasta = read_fasta(INPUT_FASTA)
    headers = list(fasta.keys())
    seq_length = len(fasta[headers[0]])
    
    print(f"Loaded {len(headers)} sequences, length {seq_length}")
    
    # Find human
    human_header = None
    for h in headers:
        if 'Human' in h or 'Homo' in h:
            human_header = h
            break
    
    human_seq = fasta[human_header]
    human_mapping = get_human_mapping(human_seq)
    
    print(f"Human mapping: {len(human_mapping)} non-gap residues")
    
    # Read conservation
    print("Reading conservation...")
    cons_df = pd.read_csv(CONSERVATION_CSV)
    
    # Create workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "Alignment"
    
    # Freeze first column (A)
    ws.freeze_panes = 'B1'
    
    # Row 1: Sequence names
    ws['A1'] = 'Sequence'
    for col_idx, header in enumerate(headers, start=1):
        ws.cell(row=1, column=col_idx + 1, value=header[:40])  # Truncate long names
    
    # Row 2-5: Conservation metrics
    metadata_rows = {
        2: 'human_position',
        3: 'vertebrate_percent',
        4: 'invertebrate_percent',
        5: 'all_percent'
    }
    
    for row_num, metric in metadata_rows.items():
        ws.cell(row=row_num, column=1, value=metric)
    
    # Row 6: Domain labels
    ws.cell(row=6, column=1, value='domain')
    
    # Fill metadata rows
    for data_col_idx, cons_col_idx in enumerate(range(1, seq_length + 1), start=2):
        cons_row = cons_df[cons_df['alignment_column'] == cons_col_idx]
        
        if not cons_row.empty:
            cons_row = cons_row.iloc[0]
            
            # Human position
            human_pos = None
            for h_pos, aln_col in human_mapping.items():
                if aln_col == cons_col_idx - 1:
                    human_pos = h_pos
                    break
            
            ws.cell(row=2, column=data_col_idx, value=human_pos)
            
            # Conservation metrics
            ws.cell(row=3, column=data_col_idx, value=cons_row['vertebrate_percent'])
            ws.cell(row=4, column=data_col_idx, value=cons_row['invertebrate_percent'])
            ws.cell(row=5, column=data_col_idx, value=cons_row['all_percent'])
            
            # Color conservation rows
            for row_num in [3, 4, 5]:
                cell = ws.cell(row=row_num, column=data_col_idx)
                if row_num == 3 and not pd.isna(cons_row['vertebrate_percent']):
                    cell.fill = PatternFill(start_color=color_for_conservation(cons_row['vertebrate_percent']), fill_type='solid')
                elif row_num == 4 and not pd.isna(cons_row['invertebrate_percent']):
                    cell.fill = PatternFill(start_color=color_for_conservation(cons_row['invertebrate_percent']), fill_type='solid')
                elif row_num == 5 and not pd.isna(cons_row['all_percent']):
                    cell.fill = PatternFill(start_color=color_for_conservation(cons_row['all_percent']), fill_type='solid')
        
        # Domain labels (row 6)
        human_pos = None
        for h_pos, aln_col in human_mapping.items():
            if aln_col == cons_col_idx - 1:
                human_pos = h_pos
                break
        
        if human_pos:
            for domain_name, (start, end) in DOMAIN_RANGES.items():
                if start <= human_pos <= end:
                    ws.cell(row=6, column=data_col_idx, value=domain_name)
                    break
    
    # Add sequences
    first_seq_row = 7
    for seq_idx, header in enumerate(headers):
        row_num = first_seq_row + seq_idx
        ws.cell(row=row_num, column=1, value=header)
        
        sequence = fasta[header]
        for col_idx, residue in enumerate(sequence, start=1):
            cell = ws.cell(row=row_num, column=col_idx + 1, value=residue)
            
            # Color by amino acid
            cell.fill = PatternFill(start_color=AA_COLORS.get(residue, 'FFFFFF'), fill_type='solid')
            
            # Check if this position should be bold (human mapping + bold positions)
            human_pos = None
            for h_pos, aln_col in human_mapping.items():
                if aln_col == col_idx - 1:
                    human_pos = h_pos
                    break
            
            if human_pos and human_pos in BOLD_POSITIONS:
                cell.font = Font(bold=True, size=11)
            else:
                cell.font = Font(size=10)
    
    last_seq_row = first_seq_row + len(headers) - 1
    
    # Set column widths
    ws.column_dimensions['A'].width = 30
    for col_idx in range(2, seq_length + 2):
        ws.column_dimensions[get_column_letter(col_idx)].width = 2.5
    
    # Set row heights
    for row in ws.iter_rows():
        for cell in row:
            cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=False)
    
    # Save
    print(f"\nSaving to {OUTPUT_FILE}...")
    wb.save(OUTPUT_FILE)
    print(f"Done! Alignment formatted into {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
