#!/usr/bin/env python3
"""Generate a 4x3 subplot (4 thresholds x 3 windows) showing all pie charts for a species.

Usage:
  python make_all_thresholds_subplot.py <results_dir> [--thresholds 4.0,5.5,6.0,6.5] [--windows 500,1000,15000]
"""
import os
import sys
import argparse
import collections
import matplotlib.pyplot as plt


def detect_encoding_and_read_lines(path):
    raw = open(path, 'rb').read()
    if raw.startswith(b"\xff\xfe\x00\x00") or raw.startswith(b"\x00\x00\xfe\xff"):
        enc = 'utf-32'
    elif raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        enc = 'utf-16'
    elif raw.startswith(b"\xef\xbb\xbf"):
        enc = 'utf-8-sig'
    else:
        for e in ('utf-8', 'utf-16', 'utf-16-le', 'utf-16-be', 'latin-1'):
            try:
                text = raw.decode(e)
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

    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    return lines


def normalize_pair_label(s):
    import re
    low = s.lower().replace('_', ' ').replace('minus', '-').replace('plus', '+')
    low = re.sub(r'[^a-z0-9+\- ]', '', low)
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
        import re
        header_norm = re.sub(r'[^a-z0-9]', '', last.lower())
        if header_norm in ('pairtype', 'pair_type', 'pair'):
            continue

        key = normalize_pair_label(last)
        cnt[key] += 1

    return cnt


def make_all_thresholds_subplot(results_dir, thresholds, windows, out_dir=None):
    results_dir = os.path.abspath(results_dir)
    if out_dir is None:
        out_dir = results_dir
    os.makedirs(out_dir, exist_ok=True)

    color_map = {
        'convergent ><': '#FF6B6B',
        'divergent <>':  '#4ECDC4',
        'tandem >>':   '#45B7D1',
        'tandem <<':   '#FFA07A',
        'other':      '#C7CEEA'
    }

    # Create 4x3 subplot grid (4 thresholds x 3 windows)
    fig, axes = plt.subplots(len(thresholds), len(windows), figsize=(15, 16))
    
    # Collect all possible labels across all data
    all_labels = set()
    for threshold in thresholds:
        for w in windows:
            tsv = os.path.join(results_dir, "tsvs", f"ctcf_pairs_{threshold}_{w}bp.tsv")
            if not os.path.isfile(tsv):
                tsv = os.path.join(results_dir, f"ctcf_pairs_{threshold}_{w}bp.tsv")
            cnt = read_pair_counts(tsv)
            all_labels.update(cnt.keys())
    
    # Order labels
    preferred_order = ['convergent ><', 'divergent <>', 'tandem >>', 'tandem <<', 'other']
    labels = [lab for lab in preferred_order if lab in all_labels]
    for lab in all_labels:
        if lab not in labels:
            labels.append(lab)

    # Plot each threshold x window combination
    for i, threshold in enumerate(thresholds):
        for j, w in enumerate(windows):
            ax = axes[i, j]
            
            # Read data
            tsv = os.path.join(results_dir, "tsvs", f"ctcf_pairs_{threshold}_{w}bp.tsv")
            if not os.path.isfile(tsv):
                tsv = os.path.join(results_dir, f"ctcf_pairs_{threshold}_{w}bp.tsv")
            cnt = read_pair_counts(tsv)
            
            sizes = [cnt.get(l, 0) for l in labels]
            total = sum(sizes)
            
            # If no data, show placeholder
            if total == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=10)
                ax.set_title(f"p={threshold}, {w}bp (n=0)", fontsize=10)
                ax.axis('off')
                continue
            
            # Filter to non-zero sizes for prettier pies
            disp_labels = [l for l, s in zip(labels, sizes) if s > 0]
            disp_sizes = [s for s in sizes if s > 0]
            disp_colors = [color_map.get(l, '#FFD166') for l in disp_labels]
            
            wedges, texts, autotexts = ax.pie(
                disp_sizes,
                labels=disp_labels if j == 0 else None,  # Only show labels in first column
                autopct='%1.1f%%',
                startangle=90,
                colors=disp_colors,
                textprops={'fontsize': 8}
            )
            for at in autotexts:
                at.set_color('white')
                at.set_fontsize(7)
            
            # Set title
            ax.set_title(f"p={threshold}, {w}bp (n={total})", fontsize=10)
    
    # Add overall title
    species_name = os.path.basename(results_dir)
    plt.suptitle(f"CTCF motif pairs: {species_name}", fontsize=16, weight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    out_png = os.path.join(out_dir, f"ctcf_all_thresholds_subplot.png")
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_png}")


def main():
    p = argparse.ArgumentParser(description='Make 4x3 subplot showing all thresholds and windows')
    p.add_argument('results_dir', help='Directory containing tsvs/ subdirectory with ctcf_pairs_*.tsv')
    p.add_argument('--thresholds', default='4.0,5.5,6.0,6.5', help='Comma-separated thresholds')
    p.add_argument('--windows', default='500,1000,15000', help='Comma-separated windows (bp)')
    p.add_argument('--out', default=None, help='Optional output directory for PNG')
    args = p.parse_args()

    thresholds = [t.strip() for t in args.thresholds.split(',') if t.strip()]
    windows = [int(w.strip()) for w in args.windows.split(',') if w.strip()]

    make_all_thresholds_subplot(args.results_dir, thresholds, windows, out_dir=args.out)


if __name__ == '__main__':
    main()
