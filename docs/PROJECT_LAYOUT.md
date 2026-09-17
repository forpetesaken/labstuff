# Project Layout

This repository is transitioning from a flat analysis workspace to an explicit code/input/result lifecycle.

## Classification rule

- **Code** is improved and versioned.
- **Inputs** are preserved and treated as read-only.
- **Results** are regenerated or archived, not hand-renamed.

## Output policy

Use `outputs/latest/` while inspecting a run. When a run matters, archive it as:

```text
outputs/archive/YYYY-MM-DD_HHMM_<name>/
  config.yaml
  run_info.txt
  figures/
  results...
```

Never edit an archived run. Create a new run instead.

## Current migration state

Root-level analysis scripts now live under `src/`. Alignment scripts live under `src/alignment_work/`; alignment inputs/intermediates are under `data/interim/alignment_work/`; historical/generated alignments are under `outputs/archive/alignment_work/`; and the CTCF workflow guide is under `docs/`. The self-contained ConSurf runtime remains under `alignment_work/ConSurf/` because its scripts use path-sensitive legacy conventions.

## Migration policy

The existing `alignment_work/` and `motif_calling/` directories contain active and historical material with different assumptions. They are intentionally not bulk-moved until each pipeline's input/output paths are updated and tested. New scripts should target the lifecycle directories, and migrations should happen one workflow at a time.
