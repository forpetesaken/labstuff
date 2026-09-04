#!/usr/bin/env bash
set -euo pipefail

CONF="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1"
CONSURF_BIN="/home/atforpetesaken/consurf_ws/consurf"

BASE="/home/atforpetesaken/consurf_runs/pds5a_0708"
mkdir -p "$BASE/work/full" "$BASE/work/vertebrates" "$BASE/work/invertebrates"
mkdir -p "$BASE/out/full" "$BASE/out/vertebrates" "$BASE/out/invertebrates"

IN_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/input_splits/PDS5A"
OUT_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/output/PDS5A"
mkdir -p "$OUT_BASE/pds5a_consurf_full" "$OUT_BASE/pds5a_consurf_vertebrates" "$OUT_BASE/pds5a_consurf_invertebrates"

export CONSURFCONF="$CONF"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/pds5a_full_0708_sanitized.fasta" \
  -SEQ_NAME Human_PDS5A \
  -Out_Dir "$BASE/out/full" \
  -w "$BASE/work/full"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/pds5a_vertebrates_sanitized.fasta" \
  -SEQ_NAME Human_PDS5A \
  -Out_Dir "$BASE/out/vertebrates" \
  -w "$BASE/work/vertebrates"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/pds5a_invertebrates_sanitized.fasta" \
  -SEQ_NAME Ciona_intestinalis \
  -Out_Dir "$BASE/out/invertebrates" \
  -w "$BASE/work/invertebrates"

cp -f "$BASE/out/full"/* "$OUT_BASE/pds5a_consurf_full/"
cp -f "$BASE/out/vertebrates"/* "$OUT_BASE/pds5a_consurf_vertebrates/"
cp -f "$BASE/out/invertebrates"/* "$OUT_BASE/pds5a_consurf_invertebrates/"

echo "PDS5A 0708 ConSurf runs complete"
