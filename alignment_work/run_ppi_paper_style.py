#!/usr/bin/env python3
"""Paper-structured PPI run from non-paired MSAs.

This script follows the staged logic of Zhang et al. (2025) on local data:
1) Build paired MSAs by species (pMSA)
2) DCA-like prefilter (APC-corrected mutual information proxy)
3) Interaction scoring from co-divergence and pMSA depth
4) High-confidence shortlist

Note: This is a practical local surrogate for RF2-PPI/AF2 stages when model
weights/GPU infrastructure are unavailable.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from pair_msas_and_score import pair_and_score, read_fasta


AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT = {aa: i for i, aa in enumerate(AA20)}
GAP_INT = 20


@dataclass
class PairStage:
    pair: str
    protein_a: str
    protein_b: str
    msa_a: Path
    msa_b: Path
    paired_msa: Path
    pair_report: Path
    shared_species: int
    divergence_corr: Optional[float]
    dca_proxy_apc_max: float
    dca_proxy_apc_q95: float
    dca_prefilter_pass: bool
    stage2_score: float
    high_confidence_80: bool
    high_confidence_90: bool


def canonical_species_key(header: str) -> str:
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


def best_record_per_species(records: Iterable[Tuple[str, str]]) -> Dict[str, Tuple[str, str]]:
    best: Dict[str, Tuple[str, str]] = {}
    for header, seq in records:
        key = canonical_species_key(header)
        nongap = sum(1 for c in seq if c != "-")
        if key not in best:
            best[key] = (header, seq)
            continue
        prev_nongap = sum(1 for c in best[key][1] if c != "-")
        if nongap > prev_nongap:
            best[key] = (header, seq)
    return best


def read_report(path: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            k, v = line.split("\t", 1)
            data[k] = v
    return data


def fasta_record_count(path: Path) -> int:
    c = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                c += 1
    return c


def tokenize_path(path: Path) -> str:
    return str(path).replace("\\", "/").lower()


def candidate_score(path: Path, protein: str) -> int:
    s = tokenize_path(path)
    score = 0

    if f"/output/{protein.lower()}/" in s:
        score += 5
    if "0708" in s:
        score += 4
    if "updated" in s:
        score += 3
    if "alignment" in s:
        score += 2

    if "consurf_full" in s:
        score -= 3
    if "consurf_vertebrates" in s:
        score -= 3
    if "consurf_invertebrates" in s:
        score -= 3
    if "consurf_run" in s:
        score -= 2
    if "human_" in s and "_msa_file" in s:
        score -= 2

    score -= len(path.parts) // 8
    return score


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
        raise FileNotFoundError(f"No MSA files found for {protein} under {input_root}")

    deep = [p for p in candidates if fasta_record_count(p) >= 10]
    if deep:
        candidates = deep

    candidates.sort(key=lambda p: candidate_score(p, protein), reverse=True)
    return candidates[0]


def encode_alignment_columnwise(seqs: Sequence[str]) -> np.ndarray:
    n = len(seqs)
    l = len(seqs[0]) if n else 0
    arr = np.full((n, l), GAP_INT, dtype=np.int16)
    for i, seq in enumerate(seqs):
        for j, aa in enumerate(seq):
            arr[i, j] = AA_TO_INT.get(aa, GAP_INT)
    return arr


def shannon_entropy_from_counts(counts: np.ndarray) -> float:
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log2(p)).sum())


def select_informative_columns(arr: np.ndarray, max_cols: int, max_gap_frac: float, min_entropy: float) -> np.ndarray:
    n, l = arr.shape
    if l == 0:
        return np.array([], dtype=np.int32)

    keep: List[Tuple[int, float]] = []
    for col in range(l):
        v = arr[:, col]
        nongap = v[v != GAP_INT]
        gap_frac = 1.0 - (len(nongap) / n)
        if gap_frac > max_gap_frac:
            continue
        counts = np.bincount(nongap, minlength=20)[:20]
        ent = shannon_entropy_from_counts(counts)
        if ent < min_entropy:
            continue
        keep.append((col, ent))

    if not keep:
        return np.array([], dtype=np.int32)

    keep.sort(key=lambda t: t[1], reverse=True)
    selected = [c for c, _ in keep[:max_cols]]
    return np.array(selected, dtype=np.int32)


def mi_pair(col_a: np.ndarray, col_b: np.ndarray, min_obs: int) -> float:
    mask = (col_a != GAP_INT) & (col_b != GAP_INT)
    if int(mask.sum()) < min_obs:
        return 0.0

    xa = col_a[mask].astype(np.int32)
    xb = col_b[mask].astype(np.int32)

    joint = np.bincount(xa * 20 + xb, minlength=400).reshape(20, 20).astype(np.float64)
    total = joint.sum()
    if total == 0:
        return 0.0

    pxy = joint / total
    px = pxy.sum(axis=1, keepdims=True)
    py = pxy.sum(axis=0, keepdims=True)

    nz = pxy > 0
    denom = px @ py
    return float((pxy[nz] * np.log(pxy[nz] / denom[nz])).sum())


def interprotein_apc_mi_score(
    seqs_a: Sequence[str],
    seqs_b: Sequence[str],
    max_cols_per_protein: int,
    max_gap_frac: float,
    min_entropy: float,
    min_obs: int,
) -> Tuple[float, float]:
    if not seqs_a or not seqs_b:
        return 0.0, 0.0

    arr_a = encode_alignment_columnwise(seqs_a)
    arr_b = encode_alignment_columnwise(seqs_b)

    cols_a = select_informative_columns(arr_a, max_cols=max_cols_per_protein, max_gap_frac=max_gap_frac, min_entropy=min_entropy)
    cols_b = select_informative_columns(arr_b, max_cols=max_cols_per_protein, max_gap_frac=max_gap_frac, min_entropy=min_entropy)

    if len(cols_a) == 0 or len(cols_b) == 0:
        return 0.0, 0.0

    mi = np.zeros((len(cols_a), len(cols_b)), dtype=np.float64)
    for i, ca in enumerate(cols_a):
        va = arr_a[:, ca]
        for j, cb in enumerate(cols_b):
            vb = arr_b[:, cb]
            mi[i, j] = mi_pair(va, vb, min_obs=min_obs)

    row_mean = mi.mean(axis=1, keepdims=True)
    col_mean = mi.mean(axis=0, keepdims=True)
    overall = float(mi.mean())
    if overall <= 1e-12:
        apc = mi
    else:
        apc = mi - (row_mean @ col_mean) / overall

    max_apc = float(np.max(apc))
    q95_apc = float(np.quantile(apc, 0.95))
    return max_apc, q95_apc


def collect_species_paired_sequences(msa_a: Path, msa_b: Path) -> Tuple[List[str], List[str], int]:
    a = best_record_per_species(read_fasta(msa_a))
    b = best_record_per_species(read_fasta(msa_b))
    shared = sorted(set(a) & set(b))

    seqs_a: List[str] = []
    seqs_b: List[str] = []
    for sp in shared:
        seqs_a.append(a[sp][1])
        seqs_b.append(b[sp][1])
    return seqs_a, seqs_b, len(shared)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run paper-structured local PPI screen")
    parser.add_argument("--input-root", type=Path, default=Path("Code/alignment_work/ConSurf/output"))
    parser.add_argument(
        "--proteins",
        nargs="+",
        default=["NIPBL", "RAD21", "SMC1", "SMC3", "PDS5A", "PDS5B", "WAPL"],
    )
    parser.add_argument("--out-dir", type=Path, default=Path("protein-protein-interactions/paper_style_run"))
    parser.add_argument("--dca-keep-fraction", type=float, default=0.4)
    parser.add_argument("--max-cols-per-protein", type=int, default=80)
    parser.add_argument("--max-gap-frac", type=float, default=0.6)
    parser.add_argument("--min-entropy", type=float, default=0.2)
    parser.add_argument("--min-obs", type=int, default=20)
    parser.add_argument("--min-shared-species", type=int, default=20)
    args = parser.parse_args()

    proteins = [p.strip().upper() for p in args.proteins if p.strip()]
    if len(proteins) < 2:
        raise ValueError("Need at least 2 proteins")

    resolved: Dict[str, Path] = {p: resolve_msa_for_protein(args.input_root, p) for p in proteins}

    out_dir = args.out_dir
    pairs_dir = out_dir / "pairs"
    pairs_dir.mkdir(parents=True, exist_ok=True)

    stages: List[PairStage] = []

    for a, b in itertools.combinations(sorted(proteins), 2):
        pair = f"{a}__{b}"
        msa_a = resolved[a]
        msa_b = resolved[b]

        out_paired = pairs_dir / f"{pair}.paired.fas"
        out_report = pairs_dir / f"{pair}.report.tsv"
        pair_and_score(msa_a, msa_b, out_paired, out_report)

        rep = read_report(out_report)
        shared = int(rep["shared_species"])
        corr = None if rep["divergence_corr_to_human"] == "NA" else float(rep["divergence_corr_to_human"])

        seqs_a, seqs_b, shared_check = collect_species_paired_sequences(msa_a, msa_b)
        if shared_check != shared:
            shared = min(shared, shared_check)

        if shared >= args.min_shared_species:
            dca_max, dca_q95 = interprotein_apc_mi_score(
                seqs_a,
                seqs_b,
                max_cols_per_protein=args.max_cols_per_protein,
                max_gap_frac=args.max_gap_frac,
                min_entropy=args.min_entropy,
                min_obs=args.min_obs,
            )
        else:
            dca_max, dca_q95 = 0.0, 0.0

        corr_val = corr if corr is not None else 0.0
        stage2_score = corr_val * math.log2(shared + 1)

        stages.append(
            PairStage(
                pair=pair,
                protein_a=a,
                protein_b=b,
                msa_a=msa_a,
                msa_b=msa_b,
                paired_msa=out_paired,
                pair_report=out_report,
                shared_species=shared,
                divergence_corr=corr,
                dca_proxy_apc_max=dca_max,
                dca_proxy_apc_q95=dca_q95,
                dca_prefilter_pass=False,
                stage2_score=stage2_score,
                high_confidence_80=False,
                high_confidence_90=False,
            )
        )

    # Stage 1 DCA-like prefilter: keep top fraction by APC max.
    sorted_dca = sorted(stages, key=lambda s: s.dca_proxy_apc_max, reverse=True)
    keep_n = max(1, int(round(len(sorted_dca) * args.dca_keep_fraction)))
    keep_pairs = {s.pair for s in sorted_dca[:keep_n]}

    for s in stages:
        s.dca_prefilter_pass = s.pair in keep_pairs

        # Stage 2/3 confidence tiers, loosely mirroring paper's stricter cutoffs.
        if s.dca_prefilter_pass and s.shared_species >= 55 and (s.divergence_corr or 0.0) >= 0.88:
            s.high_confidence_80 = True
        if s.dca_prefilter_pass and s.shared_species >= 60 and (s.divergence_corr or 0.0) >= 0.92:
            s.high_confidence_90 = True

    stages.sort(key=lambda s: (not s.dca_prefilter_pass, -(s.stage2_score), -s.dca_proxy_apc_max))

    summary_path = out_dir / "paper_style_stage_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "pair",
                "shared_species",
                "divergence_corr_to_human",
                "dca_proxy_apc_max",
                "dca_proxy_apc_q95",
                "dca_prefilter_pass",
                "stage2_score",
                "high_confidence_80",
                "high_confidence_90",
                "msa_a",
                "msa_b",
                "paired_msa",
                "pair_report",
            ]
        )
        for s in stages:
            w.writerow(
                [
                    s.pair,
                    s.shared_species,
                    "NA" if s.divergence_corr is None else f"{s.divergence_corr:.6f}",
                    f"{s.dca_proxy_apc_max:.6f}",
                    f"{s.dca_proxy_apc_q95:.6f}",
                    int(s.dca_prefilter_pass),
                    f"{s.stage2_score:.6f}",
                    int(s.high_confidence_80),
                    int(s.high_confidence_90),
                    s.msa_a,
                    s.msa_b,
                    s.paired_msa,
                    s.pair_report,
                ]
            )

    pred80_path = out_dir / "paper_style_predicted_ppi_80.tsv"
    with pred80_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair", "shared_species", "divergence_corr_to_human", "dca_proxy_apc_max", "stage2_score"])
        for s in stages:
            if s.high_confidence_80:
                w.writerow(
                    [
                        s.pair,
                        s.shared_species,
                        "NA" if s.divergence_corr is None else f"{s.divergence_corr:.6f}",
                        f"{s.dca_proxy_apc_max:.6f}",
                        f"{s.stage2_score:.6f}",
                    ]
                )

    pred90_path = out_dir / "paper_style_predicted_ppi_90.tsv"
    with pred90_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair", "shared_species", "divergence_corr_to_human", "dca_proxy_apc_max", "stage2_score"])
        for s in stages:
            if s.high_confidence_90:
                w.writerow(
                    [
                        s.pair,
                        s.shared_species,
                        "NA" if s.divergence_corr is None else f"{s.divergence_corr:.6f}",
                        f"{s.dca_proxy_apc_max:.6f}",
                        f"{s.stage2_score:.6f}",
                    ]
                )

    manifest_path = out_dir / "paper_style_run_manifest.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "resolved_msa"])
        for p in sorted(resolved):
            w.writerow([p, resolved[p]])

    print(f"Wrote staged outputs to {out_dir}")
    print(f"Total pairs: {len(stages)}")
    print(f"DCA-prefilter kept: {sum(1 for s in stages if s.dca_prefilter_pass)}")
    print(f"Predicted PPI (80-style): {sum(1 for s in stages if s.high_confidence_80)}")
    print(f"Predicted PPI (90-style): {sum(1 for s in stages if s.high_confidence_90)}")


if __name__ == "__main__":
    main()
