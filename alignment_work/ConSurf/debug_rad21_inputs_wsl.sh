#!/usr/bin/env bash
set -euo pipefail

IN_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/input_splits/RAD21"

echo "VERT_HEADERS"
grep -c '^>' "$IN_BASE/rad21_vertebrates_sanitized.fasta"

echo "INVERT_HEADERS"
grep -c '^>' "$IN_BASE/rad21_invertebrates_sanitized.fasta"

echo "VERT_FIRST"
grep '^>' "$IN_BASE/rad21_vertebrates_sanitized.fasta" | head -n 8

echo "INVERT_FIRST"
grep '^>' "$IN_BASE/rad21_invertebrates_sanitized.fasta" | head -n 8
