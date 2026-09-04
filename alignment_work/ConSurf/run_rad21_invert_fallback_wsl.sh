#!/usr/bin/env bash
set -euo pipefail

CONF="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1"
CONSURF_BIN="/home/atforpetesaken/consurf_ws/consurf"
IN_BASE="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/input_splits/RAD21"

export CONSURFCONF="$CONF"

perl "$CONSURF_BIN" \
  -MSA "$IN_BASE/rad21_invertebrates_sanitized.fasta" \
  -SEQ_NAME Ciona_intestinalis \
  -Out_Dir "/mnt/c/Users/Nat/Downloads/AIDEN" \
  -w /home/atforpetesaken/consurf_runs/rad21_reclass_20260707/work/invertebrates_fallback

echo "RAD21 invertebrate fallback run complete"
