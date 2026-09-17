#!/usr/bin/env python3
"""
Robust repeat detector for human TSVs.
- Detects file encoding (tries utf-8, utf-16, latin-1)
- Bins distances in 10 bp windows
- Reports per-orientation top bins and example pairs
- Prints concise errors only
"""
import sys
from pathlib import Path
from collections import Counter, defaultdict

BIN_SIZE = 10
THRESHOLD_PERCENT = 20.0
MAX_DISTANCE = 1000
MAX_EXAMPLES = 20

ORIENTATIONS = ['convergent', 'divergent', 'tandem_plus', 'tandem_minus']


def detect_and_read_text(path):
    b = path.read_bytes()
    # quick BOM checks
    if b.startswith(b'\xff\xfe') or b.startswith(b'\xfe\xff'):
        enc = 'utf-16'
    elif b.startswith(b'\xef\xbb\xbf'):
        enc = 'utf-8-sig'
    else:
        # try utf-8, then utf-16, then latin-1
        for trial in ('utf-8', 'utf-16', 'latin-1'):
            try:
                return b.decode(trial), trial
            except Exception:
                continue
        # fallback
        return b.decode('latin-1', errors='replace'), 'latin-1'
    try:
        return b.decode(enc), enc
    except Exception:
        return b.decode('latin-1', errors='replace'), 'latin-1'


def analyze_text(text, filename):
    counters = {o: Counter() for o in ORIENTATIONS}
    examples = {o: defaultdict(list) for o in ORIENTATIONS}
    total_pairs = 0

    lines = text.splitlines()
    if not lines:
        return None
    header = lines[0]
    for i, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 8:
            # maybe a weird line; skip quietly
            continue
        # distance might be field index 7 (0-based)
        # pair_type often at index 8 or last
        try:
            distance = int(parts[7])
        except Exception:
            # try to clean non-digit characters
            clean = ''.join(ch for ch in parts[7] if ch.isdigit())
            if clean == '':
                continue
            try:
                distance = int(clean)
            except Exception:
                continue
        if distance > MAX_DISTANCE:
            continue
        pair_type = parts[8] if len(parts) > 8 else parts[-1]
        pair_type = pair_type.lower()
        found = False
        for o in ORIENTATIONS:
            if o in pair_type:
                bin_start = (distance // BIN_SIZE) * BIN_SIZE
                counters[o][bin_start] += 1
                if len(examples[o][bin_start]) < MAX_EXAMPLES:
                    # build simple descriptor
                    chrom = parts[0]
                    try:
                        m1 = int(parts[1])
                        m2 = int(parts[4])
                        desc = f"{chrom}:{m1}-{m2} ({distance}bp)"
                    except Exception:
                        desc = f"{chrom} (dist={distance})"
                    examples[o][bin_start].append(desc)
                found = True
                total_pairs += 1
                break
        if not found:
            # orientation not recognized; skip
            continue
    results = {'total_pairs': total_pairs, 'by_orientation': {}}
    for o in ORIENTATIONS:
        total = sum(counters[o].values())
        if total == 0:
            results['by_orientation'][o] = {'total': 0, 'top_bins': [], 'flagged': False}
            continue
        top = counters[o].most_common(5)
        top_bins = []
        for bin_start, cnt in top:
            pct = cnt / total * 100
            exs = examples[o].get(bin_start, [])
            top_bins.append({'bin_start': bin_start, 'count': cnt, 'percent_of_orientation': round(pct, 1), 'examples': exs})
        top_bin = top[0]
        pct_top = top_bin[1] / total * 100
        flagged = pct_top > THRESHOLD_PERCENT
        results['by_orientation'][o] = {'total': total, 'top_bins': top_bins, 'flagged': flagged}
    return results


def find_human_tsvs(root):
    root = Path(root)
    candidates = list(root.glob('**/human*/**/tsvs/*.tsv'))
    # also explicit folder names
    candidates += list(root.glob('**/human*/tsvs/*.tsv'))
    # unique
    seen = set()
    out = []
    for p in candidates:
        if str(p) not in seen:
            out.append(p)
            seen.add(str(p))
    return out


def main():
    root = Path('.').resolve()
    tsvs = find_human_tsvs(root)
    if not tsvs:
        print('No human TSVs found.')
        return
    for tsv in tsvs:
        try:
            text, enc = detect_and_read_text(tsv)
        except Exception as e:
            print(f"{tsv}: Error reading file: {e}")
            continue
        result = analyze_text(text, tsv.name)
        if result is None:
            print(f"{tsv.name}: empty or unreadable")
            continue
        print(f"File: {tsv}  (encoding={enc})")
        print(f"  Total pairs (recognized orientations): {result['total_pairs']}")
        for o, info in result['by_orientation'].items():
            if info['total'] == 0:
                print(f"  {o}: no pairs")
                continue
            flag = 'FLAG' if info['flagged'] else ''
            print(f"  {o}: total={info['total']} {flag}")
            for tb in info['top_bins']:
                print(f"    bin {tb['bin_start']}-{tb['bin_start']+BIN_SIZE-1} : {tb['count']} ({tb['percent_of_orientation']}%)")
                if tb['examples']:
                    exs = ', '.join(tb['examples'][:3])
                    print(f"      examples: {exs}")
        print('')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted by user')
        sys.exit(1)
