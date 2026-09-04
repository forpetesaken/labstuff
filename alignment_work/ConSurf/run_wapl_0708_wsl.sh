#!/usr/bin/env bash
set -euo pipefail

BASE=/home/atforpetesaken/consurf_runs/wapl_0708_20260707
mkdir -p "$BASE/input" "$BASE/out_full" "$BASE/out_vert" "$BASE/out_inv" "$BASE/work_full" "$BASE/work_vert" "$BASE/work_inv"

cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1 "$BASE/consurfrc.local"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_full_all.fasta "$BASE/input/full.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_vertebrates_0708_sanitized.fasta "$BASE/input/vert.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_invertebrates_0708_sanitized.fasta "$BASE/input/inv.fasta"

export CONSURFCONF="$BASE/consurfrc.local"
perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/full.fasta" \
  -SEQ_NAME Human_WAPL \
  -Out_Dir "$BASE/out_full" \
  -w "$BASE/work_full"

perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/vert.fasta" \
  -SEQ_NAME Human_WAPL \
  -Out_Dir "$BASE/out_vert" \
  -w "$BASE/work_vert"

perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/inv.fasta" \
  -SEQ_NAME Ciona_intestinalis \
  -Out_Dir "$BASE/out_inv" \
  -w "$BASE/work_inv"

echo "WAPL 0708 ConSurf runs complete"
