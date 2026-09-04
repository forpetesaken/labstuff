#!/usr/bin/env python3
"""
Create one-page summaries combining all analyses for each species.
Each species gets a single image with multiple plots arranged in a grid.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

def create_species_summary(species_name, output_dir):
    """Create a one-page summary for a species combining all available plots."""
    base_dir = Path(__file__).parent
    
    # Define plot paths for this species
    # Try different naming patterns for pairs_by_chrom
    pairs_path = base_dir / 'pairsbychromosomes' / f'ctcf_pairs_by_chrom_{species_name}_5.5_500bp.png'
    if not pairs_path.exists():
        pairs_path = base_dir / 'pairsbychromosomes' / f'ctcf_pairs_by_chrom_{species_name}_500bp.png'
    
    plots = {
        'pairs_by_chrom': pairs_path if pairs_path.exists() else None,
        'strand_dist': base_dir / 'chromosome_strand_plots' / f'{species_name.lower()}_chromosome_strands_5.5.png',
        'all_thresholds': None,  # Will search for this
        'threshold_5_5': None,  # Will search for this
    }
    
    # Search for species folder
    for root_folder in ['vertebrates', 'invertebrates']:
        for folder in (base_dir / root_folder).iterdir():
            if not folder.is_dir():
                continue
            
            # Check if this is the species folder (flexible matching)
            folder_clean = folder.name.replace('_motifcalls', '').replace('_', '').lower()
            species_clean = species_name.replace('_', '').lower()
            
            if folder_clean.startswith(species_clean) or species_clean.startswith(folder_clean):
                # Found the species folder
                all_thresh = folder / 'ctcf_all_thresholds_subplot.png'
                thresh_55 = folder / 'ctcf_pairs_5_5_subplots.png'
                
                if all_thresh.exists():
                    plots['all_thresholds'] = all_thresh
                if thresh_55.exists():
                    plots['threshold_5_5'] = thresh_55
                break
    
    # Filter to only existing plots
    existing_plots = {k: v for k, v in plots.items() if v and v.exists()}
    
    if not existing_plots:
        print(f"No plots found for {species_name}")
        return
    
    print(f"Creating summary for {species_name} ({len(existing_plots)} plots)")
    
    # Create figure with grid layout
    n_plots = len(existing_plots)
    if n_plots == 1:
        fig, axes = plt.subplots(1, 1, figsize=(12, 8))
        axes = [axes]
    elif n_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    elif n_plots == 3:
        fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    
    axes = axes.flatten() if n_plots > 1 else axes
    
    # Plot titles
    titles = {
        'pairs_by_chrom': 'Pair Orientations by Chromosome',
        'strand_dist': 'Strand Distribution by Chromosome',
        'all_thresholds': 'All Thresholds Summary',
        'threshold_5_5': 'Threshold 5.5 Windows'
    }
    
    # Load and display each plot
    for idx, (key, path) in enumerate(existing_plots.items()):
        if idx >= len(axes):
            break
        
        img = mpimg.imread(path)
        axes[idx].imshow(img)
        axes[idx].axis('off')
        axes[idx].set_title(titles.get(key, key), fontsize=14, fontweight='bold', pad=10)
    
    # Hide unused subplots
    for idx in range(len(existing_plots), len(axes)):
        axes[idx].axis('off')
    
    # Main title
    species_display = species_name.replace('_', ' ')
    fig.suptitle(f'{species_display} - Analysis Summary', fontsize=18, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f'{species_name}_summary.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_path}")

def main():
    base_dir = Path(__file__).parent
    output_dir = base_dir / 'species_summaries'
    
    # Get all species from pairsbychromosomes folder
    pairs_dir = base_dir / 'pairsbychromosomes'
    if not pairs_dir.exists():
        print(f"Error: {pairs_dir} does not exist")
        return
    
    species_set = set()
    for png in pairs_dir.glob('ctcf_pairs_by_chrom_*_5.5_500bp.png'):
        # Extract species name
        name = png.stem.replace('ctcf_pairs_by_chrom_', '').replace('_5.5_500bp', '')
        species_set.add(name)
    
    print(f"Found {len(species_set)} species\n")
    
    for species in sorted(species_set):
        try:
            create_species_summary(species, output_dir)
        except Exception as e:
            print(f"Error processing {species}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Summaries saved in: {output_dir}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
