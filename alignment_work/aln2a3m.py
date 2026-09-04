#!/usr/bin/env python3
import sys
from collections import OrderedDict

def read_fasta_alignment(handle):
    seqs = OrderedDict()
    name = None
    seq_chunks = []
    for line in handle:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(seq_chunks)
            name = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(seq_chunks)
    return seqs

def read_clustal_alignment(handle):
    seqs = OrderedDict()
    started = False
    for line in handle:
        line = line.rstrip("\n")
        if not line.strip():
            # blank line: end of a block
            started = True
            continue
        if not started:
            # skip header like "CLUSTAL W ..."
            if line.upper().startswith("CLUSTAL"):
                started = True
            continue

        # skip consensus lines (those made of spaces, *, :, .)
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        name, chunk = parts[0], parts[1]
        if set(name) <= set("*:."):
            # likely a consensus line
            continue

        if name not in seqs:
            seqs[name] = []
        seqs[name].append(chunk)

    # join all chunks
    for k in list(seqs.keys()):
        seqs[k] = "".join(seqs[k])
    return seqs

def read_alignment(path):
    with open(path) as f:
        first = f.readline()
        if not first:
            raise ValueError("Empty alignment file")
        # rewind
        f.seek(0)
        if first.upper().startswith("CLUSTAL"):
            return read_clustal_alignment(f)
        else:
            return read_fasta_alignment(f)

def alignment_to_a3m(seqs):
    """
    seqs: OrderedDict{name: aligned_seq}
    Use first sequence as master. Remove columns where master has a gap.
    Return new OrderedDict{name: a3m_seq}
    """
    if not seqs:
        return seqs

    names = list(seqs.keys())
    master_name = names[0]
    master = seqs[master_name]
    L = len(master)

    # sanity check: all seqs same length
    for n, s in seqs.items():
        if len(s) != L:
            raise ValueError(f"Sequence {n} has length {len(s)} != master length {L}")

    keep_positions = [i for i, c in enumerate(master) if c != "-"]

    out = OrderedDict()
    for n in names:
        s = seqs[n]
        new_seq = "".join(s[i] for i in keep_positions)
        out[n] = new_seq
    return out

def write_a3m(seqs, handle=sys.stdout, line_width=80):
    for name, seq in seqs.items():
        handle.write(f">{name}\n")
        for i in range(0, len(seq), line_width):
            handle.write(seq[i:i+line_width] + "\n")

def main():
    if len(sys.argv) != 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} alignment.aln > output.a3m\n")
        sys.exit(1)
    aln_path = sys.argv[1]
    seqs = read_alignment(aln_path)
    a3m_seqs = alignment_to_a3m(seqs)
    write_a3m(a3m_seqs)

if __name__ == "__main__":
    main()
