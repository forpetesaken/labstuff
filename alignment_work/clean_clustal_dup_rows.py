#!/usr/bin/env python3
import sys
from collections import OrderedDict

if len(sys.argv) != 3:
    sys.stderr.write(f"Usage: {sys.argv[0]} input.clustal output.fasta\n")
    sys.exit(1)

inpath = sys.argv[1]
outpath = sys.argv[2]

seqs = OrderedDict()
started = False
block_seen = set()

with open(inpath) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            # blank line: end of block
            block_seen.clear()
            started = True
            continue
        if not started:
            if line.upper().startswith("CLUSTAL"):
                started = True
            continue
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        name, chunk = parts[0], parts[1]
        # skip consensus-like lines
        if set(name) <= set("*:."):
            continue
        # skip duplicate occurrence within the same block
        if name in block_seen:
            continue
        block_seen.add(name)
        if name not in seqs:
            seqs[name] = []
        seqs[name].append(chunk)

# join chunks
for k in list(seqs.keys()):
    seqs[k] = "".join(seqs[k])

# write FASTA
with open(outpath, 'w') as out:
    for name, seq in seqs.items():
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

print(f"Wrote cleaned alignment to {outpath}")
