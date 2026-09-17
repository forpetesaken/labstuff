from __future__ import annotations

import argparse
import gzip
import io
from pathlib import Path
from typing import Iterator, Tuple
from urllib.parse import urlparse
from urllib.request import urlopen


def _parse_fasta_handle(handle: io.TextIOBase) -> Iterator[Tuple[str, str]]:
    rec_id = None
    seq_parts = []

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if rec_id is not None:
                yield rec_id, "".join(seq_parts).upper()
            rec_id = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)

    if rec_id is not None:
        yield rec_id, "".join(seq_parts).upper()


def _is_url(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in {"http", "https"}


def fasta_records(source: str) -> Iterator[Tuple[str, str]]:
    """Yield (record_id, sequence) pairs from local FASTA or URL (optionally .gz)."""
    if _is_url(source):
        with urlopen(source) as response:
            if source.lower().endswith(".gz"):
                with gzip.GzipFile(fileobj=response) as gz_handle:
                    with io.TextIOWrapper(gz_handle, encoding="utf-8", newline="") as text_handle:
                        yield from _parse_fasta_handle(text_handle)
            else:
                with io.TextIOWrapper(response, encoding="utf-8", newline="") as text_handle:
                    yield from _parse_fasta_handle(text_handle)
        return

    path = Path(source)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from _parse_fasta_handle(handle)


def find_gaps(seq: str, min_gap: int) -> Iterator[Tuple[int, int]]:
    """Yield 0-based half-open intervals for N-runs of length >= min_gap."""
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == "N":
            start = i
            while i < n and seq[i] == "N":
                i += 1
            if i - start >= min_gap:
                yield start, i
        else:
            i += 1


def write_gap_bed(genome_fasta: str, output_bed: Path, min_gap: int) -> int:
    count = 0
    with output_bed.open("w", encoding="ascii", newline="") as out:
        for chrom, seq in fasta_records(genome_fasta):
            for start, end in find_gaps(seq, min_gap=min_gap):
                out.write(f"{chrom}\t{start}\t{end}\n")
                count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract genome gap intervals (runs of N) into BED format."
    )
    parser.add_argument(
        "genome_fasta",
        help="Input genome FASTA path or http/https URL (.fa/.fasta/.fna, optionally .gz)",
    )
    parser.add_argument("-o", "--output", type=Path, default=Path("gaps.bed"), help="Output BED path")
    parser.add_argument("--min-gap", type=int, default=10, help="Minimum N-run length to report")
    args = parser.parse_args()

    if args.min_gap < 1:
        raise ValueError("--min-gap must be >= 1")

    n_gaps = write_gap_bed(args.genome_fasta, args.output, args.min_gap)
    print(f"Wrote {n_gaps} gaps to {args.output}")


if __name__ == "__main__":
    main()
