#!/usr/bin/env python3
"""
Analyze strand ratios by chromosome for all species.
Generates CSV files with per-chromosome strand counts and ratios.
"""
import os
import sys
import csv
from collections import defaultdict

def analyze_bed_file(bed_path):
    """Read BED file and count plus/minus strands per chromosome."""
    chrom_counts = defaultdict(lambda: {'plus': 0, 'minus': 0})
    
    if not os.path.isfile(bed_path):
        return None
    
    with open(bed_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            
            chrom = fields[0]
            strand = fields[3]
            
            if strand == '+':
                chrom_counts[chrom]['plus'] += 1
            elif strand == '-':
                chrom_counts[chrom]['minus'] += 1
    
    return chrom_counts

def natural_sort_key(chrom):
    """Extract numeric part from chromosome names for natural sorting."""
    import re
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

def analyze_species(species_name, fimo_dir, output_dir, thresholds=['4.0', '5.5', '6.0', '6.5']):
    """Analyze all thresholds for a species and write CSV output."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {species_name}")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    for threshold in thresholds:
        bed_file = os.path.join(fimo_dir, f"fimo_ctcf_{threshold}.bed")
        
        if not os.path.isfile(bed_file):
            print(f"  Threshold {threshold}: File not found - {bed_file}")
            continue
        
        chrom_counts = analyze_bed_file(bed_file)
        
        if not chrom_counts:
            print(f"  Threshold {threshold}: No data")
            continue
        
        # Sort chromosomes naturally
        sorted_chroms = sorted(chrom_counts.keys(), key=natural_sort_key)
        
        # Calculate totals
        total_plus = sum(chrom_counts[c]['plus'] for c in sorted_chroms)
        total_minus = sum(chrom_counts[c]['minus'] for c in sorted_chroms)
        total_ratio = total_minus / total_plus if total_plus > 0 else 0
        
        print(f"\n  Threshold {threshold}:")
        print(f"    Total: Plus={total_plus}, Minus={total_minus}, Ratio={total_ratio:.3f}")
        
        # Write CSV
        csv_file = os.path.join(output_dir, f"strand_ratios_{threshold}.csv")
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Chromosome', 'Plus_Count', 'Minus_Count', 'Total', 'Ratio_Minus/Plus'])
            
            for chrom in sorted_chroms:
                plus = chrom_counts[chrom]['plus']
                minus = chrom_counts[chrom]['minus']
                total = plus + minus
                ratio = minus / plus if plus > 0 else float('inf')
                writer.writerow([chrom, plus, minus, total, f"{ratio:.3f}"])
            
            # Add totals row
            writer.writerow(['TOTAL', total_plus, total_minus, total_plus + total_minus, f"{total_ratio:.3f}"])
        
        print(f"    CSV saved: {csv_file}")
        
        # Print top biased chromosomes
        biased_chroms = []
        for chrom in sorted_chroms:
            plus = chrom_counts[chrom]['plus']
            minus = chrom_counts[chrom]['minus']
            total = plus + minus
            if total >= 100:  # Only consider chromosomes with at least 100 motifs
                ratio = minus / plus if plus > 0 else float('inf')
                if ratio > 1.5 or ratio < 0.67:
                    biased_chroms.append((chrom, plus, minus, ratio))
        
        if biased_chroms:
            print(f"    Biased chromosomes (ratio > 1.5 or < 0.67, n >= 100):")
            for chrom, plus, minus, ratio in sorted(biased_chroms, key=lambda x: abs(x[3] - 1.0), reverse=True)[:5]:
                print(f"      {chrom}: {plus} plus, {minus} minus, ratio={ratio:.2f}")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define all species with their paths
    species_list = [
        # Vertebrates
        ("Human", os.path.join(base_dir, "vertebrates", "human_motifcalls", "hs_fimo_out")),
        ("Arctic Lamprey", os.path.join(base_dir, "vertebrates", "alamprey_motifcalls", "ArcticLamprey_LJ_fimo_out")),
        ("Chicken", os.path.join(base_dir, "vertebrates", "chiken_motifcalls", "ga", "hs_fimo_out")),
        ("Frog", os.path.join(base_dir, "vertebrates", "frog_motifcalls", "ACfrog_ctcf", "ACfrog", "hs_fimo_out")),
        ("Bamboo Shark", os.path.join(base_dir, "vertebrates", "bambooshark_motifcalls", "BambooShark", "hs_fimo_out")),
        
        # Invertebrates
        ("Drosophila", os.path.join(base_dir, "invertebrates", "drosophila_motifcalls", "DM_fimo_out")),
        ("Ciona", os.path.join(base_dir, "invertebrates", "cion_intestines_motifcalls", "LJ_fimo_out")),
        ("Sea Urchin", os.path.join(base_dir, "invertebrates", "sea_urchin_motifcalls", "SeaUrchin", "hs_fimo_out")),
        ("Arabidopsis", os.path.join(base_dir, "invertebrates", "arabidopsis_motifcalls", "Arabidopsis", "hs_fimo_out")),
        ("Yeast", os.path.join(base_dir, "invertebrates", "yeast_motifcalls", "BakerYeast", "hs_fimo_out")),
        ("Lancelet", os.path.join(base_dir, "invertebrates", "lancelet_motifcalls", "EuropeanLancelet", "hs_fimo_out")),
        ("California Sea Hare", os.path.join(base_dir, "invertebrates", "californiaseahare_motifcalls", "CaliforniaSeaHare", "hs_fimo_out")),
        ("Chinese Liver Fluke", os.path.join(base_dir, "invertebrates", "chineseliverfluke_motifcalls", "ChineseLiverFluke", "hs_fimo_out")),
        ("Stony Coral", os.path.join(base_dir, "invertebrates", "stonycoral_motifcalls", "StonyCoral", "hs_fimo_out")),
        ("Tardigrade", os.path.join(base_dir, "invertebrates", "tardigrade_motifcalls", "Tardigrade", "hs_fimo_out")),
    ]
    
    # Create output directory
    output_base = os.path.join(base_dir, "strand_ratio_analyses")
    os.makedirs(output_base, exist_ok=True)
    
    for species_name, fimo_dir in species_list:
        if not os.path.isdir(fimo_dir):
            print(f"\nSkipping {species_name}: Directory not found - {fimo_dir}")
            continue
        
        # Create species-specific output directory
        species_output = os.path.join(output_base, species_name.replace(" ", "_").lower())
        analyze_species(species_name, fimo_dir, species_output)
    
    print(f"\n{'='*60}")
    print(f"Analysis complete! Results saved in: {output_base}")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
