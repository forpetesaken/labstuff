#!/usr/bin/env bash
set -euo pipefail

CONF="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1"
CONSURF_BIN="/home/atforpetesaken/consurf_ws/consurf"

BASE="/home/atforpetesaken/consurf_runs/rad21_stag1_style_20260707"
mkdir -p "$BASE/work/full" "$BASE/work/vertebrates" "$BASE/work/invertebrates"

IN_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/input_splits/RAD21"
OUT_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/output/RAD21"
mkdir -p "$OUT_BASE/rad21_consurf_full" "$OUT_BASE/rad21_consurf_vertebrates" "$OUT_BASE/rad21_consurf_invertebrates"

export CONSURFCONF="$CONF"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/rad21_full_reduced50.fasta" \
  -SEQ_NAME Human_RAD21 \
  -Out_Dir "$OUT_BASE/rad21_consurf_full" \
  -w "$BASE/work/full"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/rad21_vertebrates_sanitized.fasta" \
  -SEQ_NAME Human_RAD21 \
  -Out_Dir "$OUT_BASE/rad21_consurf_vertebrates" \
  -w "$BASE/work/vertebrates"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/rad21_invertebrates_sanitized.fasta" \
  -SEQ_NAME Ciona_intestinalis \
  -Out_Dir "$OUT_BASE/rad21_consurf_invertebrates" \
  -w "$BASE/work/invertebrates"

echo "RAD21 STAG1-style ConSurf runs complete"
