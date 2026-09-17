#!/usr/bin/env python3
"""
Plot chromosome-level strand distribution (+ vs -) for all species with 5.5 threshold data.
Creates bar charts showing plus and minus strand counts per chromosome.
"""
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path

def natural_sort_key(chrom):
    """Extract numeric part from chromosome names for natural sorting."""
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

def read_strand_ratios(csv_path):
    """Read strand ratio CSV and return chromosome data."""
    chromosomes = []
    plus_counts = []
    minus_counts = []
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row['Chromosome']
            if chrom == 'TOTAL':
                continue
            chromosomes.append(chrom)
            plus_counts.append(int(row['Plus_Count']))
            minus_counts.append(int(row['Minus_Count']))
    
    return chromosomes, plus_counts, minus_counts

def plot_species_strand_distribution(csv_path, species_name, output_dir):
    """Create a bar chart showing + and - strand distribution per chromosome."""
    chromosomes, plus_counts, minus_counts = read_strand_ratios(csv_path)
    
    if not chromosomes:
        print(f"No data for {species_name}")
        return
    
    # Sort by natural chromosome order
    sorted_data = sorted(zip(chromosomes, plus_counts, minus_counts), 
                        key=lambda x: natural_sort_key(x[0]))
    chromosomes, plus_counts, minus_counts = zip(*sorted_data)
    
    # Calculate percentages
    total_counts = [p + m for p, m in zip(plus_counts, minus_counts)]
    plus_pct = [100 * p / t if t > 0 else 0 for p, t in zip(plus_counts, total_counts)]
    minus_pct = [100 * m / t if t > 0 else 0 for m, t in zip(minus_counts, total_counts)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(chromosomes) * 0.4), 6))
    
    x = np.arange(len(chromosomes))
    width = 0.6
    
    # Plot stacked bars
    bars1 = ax.bar(x, plus_pct, width, label='Plus (+) strand', 
                   color='#45B7D1', alpha=0.8)
    bars2 = ax.bar(x, minus_pct, width, bottom=plus_pct, label='Minus (-) strand', 
                   color='#FFA07A', alpha=0.8)
    
    # Customize plot
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'{species_name} - Strand Distribution by Chromosome (threshold 5.5)', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right')
    ax.set_ylim(0, 100)
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.axhline(y=50, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    # Add total counts to title
    total_plus = sum(plus_counts)
    total_minus = sum(minus_counts)
    ratio = total_minus / total_plus if total_plus > 0 else 0
    ax.text(0.02, 0.98, f'Total: Plus={total_plus:,}, Minus={total_minus:,}, Ratio={ratio:.3f}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    species_slug = species_name.lower().replace(' ', '_')
    output_path = os.path.join(output_dir, f'{species_slug}_chromosome_strands_5.5.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}")

def main():
    base_dir = Path(__file__).parent
    strand_ratio_dir = base_dir / "strand_ratio_analyses"
    output_dir = base_dir / "chromosome_strand_plots"
    
    if not strand_ratio_dir.exists():
        print(f"Error: {strand_ratio_dir} does not exist")
        return
    
    print("Scanning for species with 5.5 threshold data...")
    
    # Find all species directories with strand_ratios_5.5.csv
    species_found = 0
    for species_dir in sorted(strand_ratio_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        
        csv_path = species_dir / "strand_ratios_5.5.csv"
        if not csv_path.exists():
            continue
        
        species_name = species_dir.name.replace('_', ' ').title()
        print(f"\nProcessing: {species_name}")
        
        try:
            plot_species_strand_distribution(csv_path, species_name, output_dir)
            species_found += 1
        except Exception as e:
            print(f"Error processing {species_name}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Generated {species_found} chromosome strand distribution plots")
    print(f"Plots saved in: {output_dir}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
