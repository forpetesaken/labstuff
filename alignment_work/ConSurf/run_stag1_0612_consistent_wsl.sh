#!/usr/bin/env bash
set -euo pipefail

BASE=/home/atforpetesaken/consurf_runs/stag1_0612_consistent_20260707
mkdir -p "$BASE/input" "$BASE/out_full" "$BASE/out_vert" "$BASE/out_inv" "$BASE/work_full" "$BASE/work_vert" "$BASE/work_inv"

cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1 "$BASE/consurfrc.local"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/STAG1/stag1_full_reduced50.fasta "$BASE/input/full.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/STAG1/stag1_vertebrates_sanitized.fasta "$BASE/input/vert.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/STAG1/stag1_invertebrates_sanitized.fasta "$BASE/input/inv.fasta"

export CONSURFCONF="$BASE/consurfrc.local"

perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/full.fasta" \
  -SEQ_NAME Human_STAG1 \
  -Out_Dir "$BASE/out_full" \
  -w "$BASE/work_full"

perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/vert.fasta" \
  -SEQ_NAME Human_STAG1 \
  -Out_Dir "$BASE/out_vert" \
  -w "$BASE/work_vert"

perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/inv.fasta" \
  -SEQ_NAME Branchiostoma_lanceolatum \
  -Out_Dir "$BASE/out_inv" \
  -w "$BASE/work_inv"

echo "STAG1 0612-consistent ConSurf runs complete"
