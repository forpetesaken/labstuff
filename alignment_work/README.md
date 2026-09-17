# Alignment Work

This directory now contains the self-contained ConSurf runtime and viewer assets.

- `ConSurf/`: legacy ConSurf source, configurations, input splits, and its path-sensitive output tree.
- `consurf_viewer_site/`: standalone viewer source/site assets.

Maintained alignment and PPI scripts live in `src/alignment_work/`. Reusable alignment inputs and intermediate files live under `data/interim/alignment_work/`; generated or historical results live under `outputs/`.

Run the PPI scripts from the repository root so their project-relative defaults resolve correctly:

```bash
python src/alignment_work/run_ppi_screen.py
```