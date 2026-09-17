#!/usr/bin/env python3
"""Generate 3-panel pie-chart (500 / 1000 / 15000 bp) per threshold from TSV outputs.

Usage:
  python make_threshold_subplots.py <results_dir> [--thresholds 4.0,5.5,6.0,6.5] [--windows 500,1000,15000]

The script looks for TSV files named like: ctcf_pairs_{threshold}_{window}bp.tsv
and writes a PNG named: ctcf_pairs_{threshold}_subplots.png into the same folder (or an output folder if supplied).
"""
import os
import sys
import argparse
import collections
import matplotlib.pyplot as plt


def detect_encoding_and_read_lines(path):
    # Read raw bytes and try to detect encoding using BOM or simple heuristics.
    raw = open(path, 'rb').read()
    # BOM checks
    if raw.startswith(b"\xff\xfe\x00\x00") or raw.startswith(b"\x00\x00\xfe\xff"):
        enc = 'utf-32'
    elif raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        enc = 'utf-16'
    elif raw.startswith(b"\xef\xbb\xbf"):
        enc = 'utf-8-sig'
    else:
        # try utf-8 first, then fallback to utf-16-le/utf-16-be, then latin-1
        for e in ('utf-8', 'utf-16', 'utf-16-le', 'utf-16-be', 'latin-1'):
            try:
                text = raw.decode(e)
                # if decode succeeded and no abundance of null bytes, accept
                if '\x00' in text[:1000]:
                    continue
                enc = e
                break
            except Exception:
                continue
        else:
            enc = 'latin-1'

    try:
        text = raw.decode(enc)
    except Exception:
        text = raw.decode('latin-1', errors='replace')

    # split into lines normalized
    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    return lines


def normalize_pair_label(s):
    import re
    low = s.lower().replace('_', ' ').replace('minus', '-').replace('plus', '+')
    low = re.sub(r'[^a-z0-9+\- ]', '', low)
    # Normalize all possible variants
    if 'converg' in low:
        return 'convergent ><'
    if 'diverg' in low:
        return 'divergent <>'
    if ('tandem' in low and ('+' in low or 'plus' in low)) or ('tandem plus' in low) or ('tandem +' in low):
        return 'tandem >>'
    if ('tandem' in low and ('-' in low or 'minus' in low)) or ('tandem minus' in low) or ('tandem -' in low):
        return 'tandem <<'
    if low in ['tandem', 'tandemplus', 'tandem+']:
        return 'tandem >>'
    if low in ['tandemminus', 'tandem-']:
        return 'tandem <<'
    low = low.strip()
    return low or 'other'


def read_pair_counts(tsv_path):
    """Read TSV and count pair_type (last column). Returns Counter."""
    cnt = collections.Counter()
    if not os.path.isfile(tsv_path):
        return cnt

    lines = detect_encoding_and_read_lines(tsv_path)
    for line in lines:
        if not line or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 1:
            continue
        last = fields[-1].strip()
        # detect header-like last fields by normalizing and matching
        import re
        header_norm = re.sub(r'[^a-z0-9]', '', last.lower())
        if header_norm in ('pairtype', 'pair_type', 'pair'):
            continue

        key = normalize_pair_label(last)
        cnt[key] += 1

    return cnt


def make_subplot_for_threshold(results_dir, threshold, windows, out_dir=None):
    results_dir = os.path.abspath(results_dir)
    if out_dir is None:
        out_dir = results_dir
    os.makedirs(out_dir, exist_ok=True)

    counts_per_window = []
    labels_union = []
    for w in windows:
        # Try tsvs/ subdirectory first, then root directory
        tsv = os.path.join(results_dir, "tsvs", f"ctcf_pairs_{threshold}_{w}bp.tsv")
        if not os.path.isfile(tsv):
            tsv = os.path.join(results_dir, f"ctcf_pairs_{threshold}_{w}bp.tsv")
        cnt = read_pair_counts(tsv)
        counts_per_window.append(cnt)
        for k in cnt.keys():
            if k not in labels_union:
                labels_union.append(k)

    # normalize and enforce a stable label order and color mapping
    preferred_order = ['convergent ><', 'divergent <>', 'tandem >>', 'tandem <<', 'other']
    # determine labels present, preserving preferred order
    labels = [lab for lab in preferred_order if lab in labels_union]
    # append any unexpected labels after the preferred ones
    for lab in labels_union:
        if lab not in labels:
            labels.append(lab)

    # colors (map to preferred_order)
    color_map = {
        'convergent ><': '#FF6B6B',
        'divergent <>':  '#4ECDC4',
        'tandem >>':   '#45B7D1',
        'tandem <<':   '#FFA07A',
        'other':      '#C7CEEA'
    }
    base_colors = [color_map.get(l, '#FFD166') for l in labels]

    fig, axes = plt.subplots(1, len(windows), figsize=(5 * len(windows), 5))
    if len(windows) == 1:
        axes = [axes]

    for ax, w, cnt in zip(axes, windows, counts_per_window):
        sizes = [cnt.get(l, 0) for l in labels]
        total = sum(sizes)
        # if all zeros, draw a placeholder
        if total == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=12)
            ax.set_title(f"{w} bp\n(n=0)")
            ax.axis('off')
            continue

        # trim labels/sizes to non-zero for prettier pies
        disp_labels = [l for l, s in zip(labels, sizes) if s > 0]
        disp_sizes = [s for s in sizes if s > 0]
        disp_colors = [color_map.get(l, '#FFD166') for l in disp_labels]

        wedges, texts, autotexts = ax.pie(
            disp_sizes,
            labels=disp_labels,
            autopct='%1.1f%%',
            startangle=90,
            colors=disp_colors,
            textprops={'fontsize': 10}
        )
        for at in autotexts:
            at.set_color('white')
            at.set_fontsize(9)

        ax.set_title(f"{w} bp (n={total})", fontsize=12)

    plt.suptitle(f"CTCF motif {threshold} — windows: {', '.join(map(str, windows))}", fontsize=14, weight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # make a safe filename for the threshold (avoid '.' in filenames)
    safe_thr = str(threshold).replace('.', '_')
    out_png = os.path.join(out_dir, f"ctcf_pairs_{safe_thr}_subplots.png")
    plt.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_png}")


def main():
    p = argparse.ArgumentParser(description='Make subplot pie charts for thresholds')
    p.add_argument('results_dir', help='Directory containing ctcf_pairs_*.tsv')
    p.add_argument('--thresholds', default='4.0,5.5,6.0,6.5', help='Comma-separated thresholds')
    p.add_argument('--windows', default='500,1000,15000', help='Comma-separated windows (bp)')
    p.add_argument('--out', default=None, help='Optional output directory for PNGs')
    args = p.parse_args()

    thresholds = [t.strip() for t in args.thresholds.split(',') if t.strip()]
    windows = [int(w.strip()) for w in args.windows.split(',') if w.strip()]

    for thr in thresholds:
        make_subplot_for_threshold(args.results_dir, thr, windows, out_dir=args.out)


if __name__ == '__main__':
    main()
