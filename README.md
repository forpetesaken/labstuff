# AIDEN Work

Research and analysis workspace for genome architecture, motif calling, alignments, and ConSurf workflows.

## Project layout

- `src/`: maintained root-level analysis code.
- `configs/`: experiment parameters and run configuration files.
- `data/raw/`: source data; never modify in place.
- `data/interim/`: temporary transformed inputs.
- `data/processed/`: reusable processed datasets.
- `outputs/latest/`: outputs currently under inspection.
- `outputs/archive/`: immutable snapshots of meaningful runs.
- `notebooks/`: exploratory notebooks.
- `docs/`: methods, notes, and figure documentation.
- `scripts/`: one-off utilities and run helpers.

The existing `motif_calling/` tree is preserved as a legacy project area. Alignment work has been reorganized: maintained scripts are in `src/alignment_work/`, intermediate inputs are in `data/interim/alignment_work/`, and historical/generated alignment results are under `outputs/`. The self-contained, path-sensitive ConSurf runtime remains under `alignment_work/ConSurf/`.

The root-level scripts have been moved into `src/`; one-off utilities are in `scripts/`, the CTCF workflow guide is in `docs/`, the large loop-anchor BED is in `data/raw/`, and the generated BEDPE report is in `outputs/latest/`.

## Run lifecycle

1. Edit code under `src/` or the existing legacy code area.
2. Change a config under `configs/`.
3. Write active results to `outputs/latest/`.
4. Archive a meaningful run with `scripts/archive_run.sh`.
5. Commit code and configuration, not large generated data.

Example:

```bash
scripts/archive_run.sh configs/example.yaml outputs/latest peak_calling
```

The archive helper records the timestamp, Git commit, input config, command, and output files. Archived directories should not be edited after creation.
