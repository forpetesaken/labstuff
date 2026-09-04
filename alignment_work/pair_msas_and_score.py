#!/usr/bin/env python3
"""Pair two unpaired MSAs by species and compute a first-pass coevolution proxy score.

This script is intended for quick triage of candidate PPIs from existing single-protein MSAs.
It is not a replacement for RF2-PPI / AF-Multimer confidence modeling.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


FASTARecord = Tuple[str, str]


def read_fasta(path: Path) -> List[FASTARecord]:
    records: List[FASTARecord] = []
    header = None
    seq_chunks: List[str] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks).upper()))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        records.append((header, "".join(seq_chunks).upper()))

    return records


def _canonical_species_key(header: str) -> str:
    # Prefer UniProt-style OS=Species name OX=TaxId when available.
    m = re.search(r"\bOS=([^=]+?)\s+OX=", header)
    if m:
        species = m.group(1).strip()
    else:
        species = header.split()[0].strip()

    species = species.replace(" ", "_").lower()
    species = re.sub(r"[^a-z0-9_]+", "_", species)
    species = re.sub(r"_(ncbi|dnazoo)$", "", species)
    species = re.sub(r"_[0-9]+$", "", species)
    species = re.sub(r"_+", "_", species).strip("_")
    return species


def _best_record_per_species(records: Iterable[FASTARecord]) -> Dict[str, FASTARecord]:
    best: Dict[str, FASTARecord] = {}
    for header, seq in records:
        key = _canonical_species_key(header)
        # Keep the least gappy sequence when duplicates exist for one species.
        nongap = sum(1 for c in seq if c != "-")
        if key not in best:
            best[key] = (header, seq)
            continue
        prev_nongap = sum(1 for c in best[key][1] if c != "-")
        if nongap > prev_nongap:
            best[key] = (header, seq)
    return best


def _find_human_key(keys: Iterable[str]) -> str | None:
    key_set = set(keys)
    for candidate in ("homo_sapiens", "human"):
        if candidate in key_set:
            return candidate
    return None


def _ungapped_identity(seq_a: str, seq_b: str) -> float | None:
    matches = 0
    denom = 0
    for aa, bb in zip(seq_a, seq_b):
        if aa == "-" or bb == "-":
            continue
        denom += 1
        if aa == bb:
            matches += 1
    if denom == 0:
        return None
    return matches / denom


def _pearson(values_x: List[float], values_y: List[float]) -> float | None:
    n = len(values_x)
    if n < 3:
        return None
    mean_x = sum(values_x) / n
    mean_y = sum(values_y) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(values_x, values_y))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in values_x))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in values_y))
    den = den_x * den_y
    if den == 0:
        return None
    return num / den


def _shannon_entropy(seq: str) -> float:
    nongap = [aa for aa in seq if aa != "-"]
    if not nongap:
        return 0.0
    counts = Counter(nongap)
    total = sum(counts.values())
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy


def pair_and_score(msa_a: Path, msa_b: Path, out_paired: Path, out_report: Path) -> None:
    recs_a = read_fasta(msa_a)
    recs_b = read_fasta(msa_b)

    species_a = _best_record_per_species(recs_a)
    species_b = _best_record_per_species(recs_b)

    shared = sorted(set(species_a) & set(species_b))
    if not shared:
        raise ValueError("No shared species found between the two MSAs.")

    human_key = _find_human_key(shared)

    out_paired.parent.mkdir(parents=True, exist_ok=True)
    out_report.parent.mkdir(parents=True, exist_ok=True)

    with out_paired.open("w", encoding="utf-8") as handle:
        for sp in shared:
            _, seq_a = species_a[sp]
            _, seq_b = species_b[sp]
            handle.write(f">{sp}\n")
            handle.write(f"{seq_a}:{seq_b}\n")

    # Coevolution proxy: correlation of per-species divergence from human in the two proteins.
    ids_a: List[float] = []
    ids_b: List[float] = []

    if human_key is not None:
        human_a = species_a[human_key][1]
        human_b = species_b[human_key][1]
        for sp in shared:
            seq_a = species_a[sp][1]
            seq_b = species_b[sp][1]
            ida = _ungapped_identity(seq_a, human_a)
            idb = _ungapped_identity(seq_b, human_b)
            if ida is None or idb is None:
                continue
            ids_a.append(1.0 - ida)
            ids_b.append(1.0 - idb)

    corr = _pearson(ids_a, ids_b) if ids_a and ids_b else None

    # Secondary sanity metrics.
    human_entropy_a = _shannon_entropy(species_a[human_key][1]) if human_key else float("nan")
    human_entropy_b = _shannon_entropy(species_b[human_key][1]) if human_key else float("nan")

    with out_report.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"msa_a\t{msa_a}\n")
        handle.write(f"msa_b\t{msa_b}\n")
        handle.write(f"unique_species_a\t{len(species_a)}\n")
        handle.write(f"unique_species_b\t{len(species_b)}\n")
        handle.write(f"shared_species\t{len(shared)}\n")
        handle.write(f"human_key\t{human_key if human_key else 'NA'}\n")
        handle.write(f"paired_depth\t{len(shared)}\n")
        handle.write(f"divergence_corr_to_human\t{corr if corr is not None else 'NA'}\n")
        handle.write(f"human_seq_entropy_a\t{human_entropy_a}\n")
        handle.write(f"human_seq_entropy_b\t{human_entropy_b}\n")



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pair two non-paired MSAs by species and compute a quick coevolution proxy."
    )
    parser.add_argument("--msa-a", required=True, type=Path, help="Path to protein A MSA FASTA")
    parser.add_argument("--msa-b", required=True, type=Path, help="Path to protein B MSA FASTA")
    parser.add_argument(
        "--out-paired",
        required=True,
        type=Path,
        help="Output paired MSA FASTA path (concatenated as seqA:seqB)",
    )
    parser.add_argument(
        "--out-report",
        required=True,
        type=Path,
        help="Output TSV report with pairing and score summary",
    )
    args = parser.parse_args()

    pair_and_score(args.msa_a, args.msa_b, args.out_paired, args.out_report)


if __name__ == "__main__":
    main()
