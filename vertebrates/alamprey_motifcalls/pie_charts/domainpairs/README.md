# domainpairs pie charts

**Script:** ctcfbeds_domainpairs.py

**Usage:**
- python ctcfbeds_domainpairs.py motifs.bed domains.bedpe [WINDOW_BP] [MODE] [LABEL] [OUT_DIR]

**MODE options:**
- nearest (default): closest motif to boundary
- interior: most interior motif in window
- discard: skip if multiple motifs in window
- strongest: highest score motif in window

**Notes:**
- WINDOW_BP is the boundary window size (default 10000).
- Charts are written into this folder.
