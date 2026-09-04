# boundary_pairs pie charts

**Script:** ctcf_boundary_pairs.py

**Usage:**
- python ctcf_boundary_pairs.py motifs.bed domains.bedpe [MAX_DIST] [MODE] [LABEL] [OUT_DIR]

**MODE options:**
- all (default): use all motifs within window
- strongest2: use top 2 scored motifs within window

**Notes:**
- MAX_DIST is the boundary window size (default 10000).
- Charts are written into this folder.
