# human insulation-boundary analysis

This folder is for CTCF analysis at insulation-score boundaries (human), parallel to your domain-pair workflow.

## Files
- `count_ctcf_at_insulation_boundaries.py`: counts CTCF motifs around boundary midpoints.
- `outputs/`: suggested location for per-boundary TSV outputs.

## Input expectations
- Motifs BED: `hs_fimo_out/fimo_ctcf_5.5.bed` (or another threshold BED)
- Boundaries BED: insulation-boundary BED with at least 3 columns (`chrom start end`)

## Example run
From `ctcf_bedmotif_analysis/vertebrates/human_motifcalls/insulation_boundaries`:

```bash
python count_ctcf_at_insulation_boundaries.py \
  ../hs_fimo_out/fimo_ctcf_5.5.bed \
  ../hs_2000_blocks \
  --window-bp 10000 \
  --out outputs/hs_2000_insulation_ctcf_counts.tsv
```

## Notes
- Boundary position is midpoint: `(start + end) // 2`.
- Window is symmetric around midpoint: `[midpoint - window_bp, midpoint + window_bp]`.
- If your insulation boundaries are in a different file, replace `../hs_2000_blocks` with that path.
