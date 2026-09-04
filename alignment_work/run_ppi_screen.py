#!/usr/bin/env python3
"""Run a lightweight PPI screen from non-paired MSAs.

Workflow:
1) Resolve one input MSA per protein.
2) Build paired MSAs for all protein pairs.
3) Score each pair using pair_msas_and_score.py.
4) Write ranked summary tables.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from pair_msas_and_score import pair_and_score


@dataclass
class PairResult:
    protein_a: str
    protein_b: str
    msa_a: Path
    msa_b: Path
    shared_species: int
    corr: Optional[float]
    tier: str
    screen_score: float
    paired_msa_path: Path
    report_path: Path


def _tokenize_path(path: Path) -> str:
    return str(path).replace("\\", "/").lower()


def _candidate_score(path: Path, protein: str) -> int:
    s = _tokenize_path(path)
    score = 0

    # Prefer direct protein folders and curated snapshots.
    if f"/output/{protein.lower()}/" in s:
        score += 5
    if "0708" in s:
        score += 4
    if "updated" in s:
        score += 3
    if "alignment" in s:
        score += 2

    # Avoid sub-collections used for ConSurf visualization buckets.
    if "consurf_full" in s:
        score -= 3
    if "consurf_vertebrates" in s:
        score -= 3
    if "consurf_invertebrates" in s:
        score -= 3
    if "consurf_run" in s:
        score -= 2

    # Generic files named Human_*_msa_file are often intermediate snapshots.
    if "human_" in s and "_msa_file" in s:
        score -= 2

    # Prefer shorter paths when tied (usually less-derived files).
    score -= len(path.parts) // 8
    return score


def _fasta_record_count(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def resolve_msa_for_protein(input_root: Path, protein: str) -> Path:
    aliases = {
        "PDS5A": ["PDS5A", "PSD5A"],
    }
    tokens = aliases.get(protein.upper(), [protein])

    patterns = [
        "**/*{token}*.fas",
        "**/*{token}*.fasta",
        "**/*{token}*.fa",
        "**/*{token}*.aln",
        "**/*{token}*.a3m",
    ]

    candidates: List[Path] = []
    for token in tokens:
        for pat in patterns:
            candidates.extend(input_root.glob(pat.format(token=token)))

    candidates = [p for p in candidates if p.is_file()]
    if not candidates:
        raise FileNotFoundError(f"No MSA files found for protein {protein} under {input_root}")

    # Prefer multi-species MSAs and avoid one-sequence files.
    deep_candidates = [p for p in candidates if _fasta_record_count(p) >= 10]
    if deep_candidates:
        candidates = deep_candidates

    candidates.sort(key=lambda p: _candidate_score(p, protein), reverse=True)
    return candidates[0]


def _read_report(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        next(handle)  # skip header
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            k, v = line.split("\t", 1)
            data[k] = v
    return data


def _screen_score(shared_species: int, corr: Optional[float]) -> float:
    corr_value = corr if corr is not None else 0.0
    return corr_value * math.log2(shared_species + 1)


def _tier(shared_species: int, corr: Optional[float]) -> str:
    if corr is None:
        return "D"
    if shared_species >= 60 and corr >= 0.90:
        return "A"
    if shared_species >= 50 and corr >= 0.85:
        return "B"
    if shared_species >= 30 and corr >= 0.75:
        return "C"
    return "D"


def run_screen(
    input_root: Path,
    proteins: Iterable[str],
    out_dir: Path,
) -> List[PairResult]:
    proteins = [p.strip().upper() for p in proteins if p.strip()]
    if len(proteins) < 2:
        raise ValueError("Need at least two proteins to run all-pairs screen.")

    resolved: Dict[str, Path] = {}
    for protein in proteins:
        resolved[protein] = resolve_msa_for_protein(input_root, protein)

    out_dir.mkdir(parents=True, exist_ok=True)
    pairs_dir = out_dir / "pairs"
    pairs_dir.mkdir(parents=True, exist_ok=True)

    results: List[PairResult] = []

    for a, b in itertools.combinations(sorted(proteins), 2):
        pair_name = f"{a}__{b}"
        out_paired = pairs_dir / f"{pair_name}.paired.fas"
        out_report = pairs_dir / f"{pair_name}.report.tsv"

        pair_and_score(resolved[a], resolved[b], out_paired, out_report)
        rep = _read_report(out_report)

        shared = int(rep["shared_species"])
        corr = None if rep["divergence_corr_to_human"] == "NA" else float(rep["divergence_corr_to_human"])
        score = _screen_score(shared, corr)
        tier = _tier(shared, corr)

        results.append(
            PairResult(
                protein_a=a,
                protein_b=b,
                msa_a=resolved[a],
                msa_b=resolved[b],
                shared_species=shared,
                corr=corr,
                tier=tier,
                screen_score=score,
                paired_msa_path=out_paired,
                report_path=out_report,
            )
        )

    results.sort(
        key=lambda r: (
            -r.screen_score,
            -r.shared_species,
            -(r.corr if r.corr is not None else -1.0),
        )
    )

    return results


def write_outputs(results: List[PairResult], out_dir: Path) -> None:
    summary = out_dir / "ppi_screen_summary.tsv"
    with summary.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "pair",
                "tier",
                "screen_score",
                "shared_species",
                "divergence_corr_to_human",
                "msa_a",
                "msa_b",
                "paired_msa",
                "pair_report",
            ]
        )
        for r in results:
            pair = f"{r.protein_a}__{r.protein_b}"
            w.writerow(
                [
                    pair,
                    r.tier,
                    f"{r.screen_score:.6f}",
                    r.shared_species,
                    "NA" if r.corr is None else f"{r.corr:.6f}",
                    r.msa_a,
                    r.msa_b,
                    r.paired_msa_path,
                    r.report_path,
                ]
            )

    top = out_dir / "ppi_screen_top_hits.tsv"
    with top.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair", "tier", "screen_score", "shared_species", "divergence_corr_to_human"])
        for r in results[:20]:
            pair = f"{r.protein_a}__{r.protein_b}"
            w.writerow(
                [
                    pair,
                    r.tier,
                    f"{r.screen_score:.6f}",
                    r.shared_species,
                    "NA" if r.corr is None else f"{r.corr:.6f}",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all-pairs PPI screen from non-paired MSAs")
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("Code/alignment_work/ConSurf/output"),
        help="Root folder containing single-protein MSAs",
    )
    parser.add_argument(
        "--proteins",
        nargs="+",
        default=["NIPBL", "RAD21", "SMC1", "SMC3", "PDS5A", "PDS5B", "WAPL"],
        help="Protein names to include in the all-pairs run",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("protein-protein-interactions/ppi_screen_run"),
        help="Output folder for paired MSAs and score tables",
    )
    args = parser.parse_args()

    results = run_screen(args.input_root, args.proteins, args.out_dir)
    write_outputs(results, args.out_dir)

    print(f"Wrote {len(results)} pair results to {args.out_dir}")


if __name__ == "__main__":
    main()
