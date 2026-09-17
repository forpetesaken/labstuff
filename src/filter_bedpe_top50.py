#!/usr/bin/env python3
"""
Filter BEDPE rows to top percentile based on the score column after color.
Input: 5000_blocks.bedpe (or similar)
Output: <input>_topXX.bedpe
"""
import sys
import os

if len(sys.argv) < 2:
    print("Usage: filter_bedpe_top50.py input.bedpe [percentile] [output.bedpe]", file=sys.stderr)
    sys.exit(1)

input_path = sys.argv[1]
percentile = 50.0
output_path = None

if len(sys.argv) >= 3:
    try:
        percentile = float(sys.argv[2])
        output_path = sys.argv[3] if len(sys.argv) >= 4 else None
    except ValueError:
        output_path = sys.argv[2] if len(sys.argv) >= 3 else None

scores = []
rows = []

with open(input_path) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        # score after color is column 12 (0-based index 11)
        if len(fields) <= 11:
            continue
        try:
            score = float(fields[11])
        except ValueError:
            continue
        scores.append(score)
        rows.append((score, line.rstrip('\n')))

if not scores:
    print("No valid rows found.", file=sys.stderr)
    sys.exit(1)

scores.sort()
# percentile threshold (top X%)
scores.sort()
if percentile <= 0 or percentile >= 100:
    print("Percentile must be between 0 and 100.", file=sys.stderr)
    sys.exit(1)

cut_index = int(len(scores) * (percentile / 100.0))
cut_index = max(0, min(cut_index, len(scores) - 1))
threshold = scores[cut_index]

if not output_path:
    base, ext = os.path.splitext(input_path)
    output_path = f"{base}_top{int(percentile)}{ext if ext else '.bedpe'}"

kept = 0
with open(output_path, 'w') as out:
    for score, line in rows:
        if score >= threshold:
            out.write(line + '\n')
            kept += 1

print(f"Total rows: {len(rows)}", file=sys.stderr)
print(f"Threshold (top {percentile}th percentile): {threshold}", file=sys.stderr)
print(f"Kept rows: {kept}", file=sys.stderr)
print(f"Output: {output_path}", file=sys.stderr)
