#!/usr/bin/env python3
"""
Plot chromosome-level CTCF pair orientation distribution for all species.
Creates stacked bar charts showing convergent, divergent, tandem_plus, tandem_minus percentages per chromosome.
Uses 5.5 threshold, 500bp window data.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path
from collections import defaultdict, Counter

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

def read_pair_orientations(tsv_path):
    """Read TSV and count pair orientations per chromosome."""
    chrom_counts = defaultdict(Counter)
    
    with open(tsv_path, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            pair_type = fields[8]
            chrom_counts[chrom][pair_type] += 1
    
    return chrom_counts

def normalize_pair_label(label):
    """Normalize pair type labels."""
    label = label.strip().lower()
    if 'convergent' in label or '><' in label:
        return 'convergent ><'
    elif 'divergent' in label or '<>' in label:
        return 'divergent <>'
    elif 'tandem' in label and ('plus' in label or '>>' in label):
        return 'tandem_plus >>'
    elif 'tandem' in label and ('minus' in label or '<<' in label):
        return 'tandem_minus <<'
    return label

def plot_species_pair_orientations(tsv_path, species_name, output_dir):
    """Create a stacked bar chart showing pair orientation percentages per chromosome."""
    chrom_counts = read_pair_orientations(tsv_path)
    
    if not chrom_counts:
        print(f"No data for {species_name}")
        return
    
    # Sort chromosomes
    chromosomes = sorted(chrom_counts.keys(), key=natural_sort_key)
    
    # Filter out chromosomes with very few pairs (< 10)
    chromosomes = [c for c in chromosomes if sum(chrom_counts[c].values()) >= 10]
    
    if not chromosomes:
        print(f"No chromosomes with sufficient data for {species_name}")
        return
    
    # Define pair types and colors
    pair_types = ['convergent ><', 'divergent <>', 'tandem_plus >>', 'tandem_minus <<']
    colors = {
        'convergent ><': '#FF6B6B',
        'divergent <>': '#4ECDC4',
        'tandem_plus >>': '#45B7D1',
        'tandem_minus <<': '#FFA07A'
    }
    
    # Calculate percentages for each chromosome
    percentages = {pt: [] for pt in pair_types}
    for chrom in chromosomes:
        total = sum(chrom_counts[chrom].values())
        for pt in pair_types:
            count = chrom_counts[chrom][pt]
            pct = 100 * count / total if total > 0 else 0
            percentages[pt].append(pct)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(12, len(chromosomes) * 0.4), 6))
    
    x = np.arange(len(chromosomes))
    width = 0.6
    
    # Plot stacked bars
    bottom = np.zeros(len(chromosomes))
    for pt in pair_types:
        ax.bar(x, percentages[pt], width, label=pt, 
               color=colors[pt], alpha=0.8, bottom=bottom)
        bottom += np.array(percentages[pt])
    
    # Customize plot
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'{species_name} - CTCF Pair Orientations by Chromosome\n(threshold 5.5, window 500bp)', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(chromosomes, rotation=45, ha='right')
    ax.set_ylim(0, 100)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add total counts info
    total_pairs = sum(sum(chrom_counts[c].values()) for c in chromosomes)
    ax.text(0.02, 0.98, f'Total pairs: {total_pairs:,}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    species_slug = species_name.lower().replace(' ', '_')
    output_path = os.path.join(output_dir, f'{species_slug}_pair_orientations_5.5_500bp.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}")

def main():
    base_dir = Path(__file__).parent
    output_dir = base_dir / "chromosome_strand_plots"
    
    # Define species mapping
    species_paths = [
        ("vertebrates/human_motifcalls", "Human"),
        ("vertebrates/alamprey_motifcalls", "Arctic Lamprey"),
        ("vertebrates/chiken_motifcalls", "Chicken"),
        ("vertebrates/frog_motifcalls", "Frog"),
        ("vertebrates/bambooshark_motifcalls", "Bamboo Shark"),
        ("invertebrates/drosophila_motifcalls", "Drosophila"),
        ("invertebrates/cion_intestines_motifcalls", "Ciona"),
        ("invertebrates/sea_urchin_motifcalls", "Sea Urchin"),
        ("invertebrates/arabidopsis_motifcalls", "Arabidopsis"),
        ("invertebrates/yeast_motifcalls", "Yeast"),
        ("invertebrates/lancelet_motifcalls", "Lancelet"),
        ("invertebrates/californiaseahare_motifcalls", "California Sea Hare"),
        ("invertebrates/chineseliverfluke_motifcalls", "Chinese Liver Fluke"),
        ("invertebrates/stonycoral_motifcalls", "Stony Coral"),
        ("invertebrates/tardigrade_motifcalls", "Tardigrade"),
    ]
    
    print("Scanning for species with 5.5 threshold, 500bp window data...")
    
    species_found = 0
    for species_path, species_name in species_paths:
        tsv_path = base_dir / species_path / "tsvs" / "ctcf_pairs_5.5_500bp.tsv"
        
        if not tsv_path.exists():
            print(f"Skipping {species_name}: TSV not found")
            continue
        
        print(f"\nProcessing: {species_name}")
        
        try:
            plot_species_pair_orientations(tsv_path, species_name, output_dir)
            species_found += 1
        except Exception as e:
            print(f"Error processing {species_name}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Generated {species_found} pair orientation plots")
    print(f"Plots saved in: {output_dir}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
