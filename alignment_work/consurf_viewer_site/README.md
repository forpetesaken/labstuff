# ConSurf Interactive Viewer Site

This folder is a ready-to-publish static website and can be used as its own git repo.
It includes top tabs for:

- RAD21
- STAG1
- STAG2
- CTCF

## Files

- `index.html` - interactive ConSurf viewer (standalone, works without internet once loaded)
- `build_site.py` - builds the site from latest ConSurf `.grades` files
- `update_viewer.ps1` - rebuilds `index.html` from the latest RAD21 ConSurf grade outputs

## Keep Site Updated

Run this from inside this folder whenever new ConSurf outputs are generated:

```powershell
.\update_viewer.ps1
```

This rebuilds tabs/datasets from current files under `ConSurf/output/`.

Then commit and push:

```powershell
git add index.html
git commit -m "Update viewer data"
git push
```

## Gap Handling Notes

ConSurf is run in MSA mode using the provided alignment as-is.

- Full and vertebrate views use the human sequence as the query.
- Invertebrate views use the selected invertebrate reference/query sequence.
- ConSurf reports one conservation score per non-gap residue in the query sequence for that run.
- Gaps in other sequences are treated as missing data, not amino acid substitutions.
- If many non-query sequences are gapped at a query position, that position is still scored, but with less information.
- If the query itself has no residue at a column, that column is not scored at all.

### Human Presence Track

For full and invertebrate views, the site can display a `Human present` track below the plot.

- This indicates whether a scored query position has a human residue aligned at that same column.
- Hovering the track shows the aligned human residue and human residue number when available.
- This is especially useful for invertebrate plots, where conservation can appear over regions that do not actually align to human residues.

## Quick publish options

### Option 1: GitHub Pages

1. Create a new GitHub repository.
2. Upload the contents of this folder (`index.html`).
3. In repository settings, open **Pages**.
4. Set source to **Deploy from a branch** and choose `main` + `/ (root)`.
5. Save. GitHub will provide a public URL.

### Option 2: Netlify (drag and drop)

1. Go to Netlify Drop: https://app.netlify.com/drop
2. Drag this `consurf_viewer_site` folder into the page.
3. Netlify gives you a public share link immediately.

### Option 3: Share as a single file

You can also share just `index.html` directly. People can open it in a browser.
