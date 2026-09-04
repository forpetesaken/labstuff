#!/usr/bin/env bash
set -euo pipefail

BASE=/home/atforpetesaken/consurf_runs/wapl_20260707
mkdir -p "$BASE/input" "$BASE/out_full" "$BASE/out_vert" "$BASE/out_inv" "$BASE/work_full" "$BASE/work_vert" "$BASE/work_inv"

cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1 "$BASE/consurfrc.local"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_full_reduced50.fasta "$BASE/input/full.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_vertebrates_sanitized.fasta "$BASE/input/vert.fasta"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/input_splits/WAPL/wapl_invertebrates_sanitized.fasta "$BASE/input/inv.fasta"

# Relax fragment threshold for small invertebrate subsets.
sed 's/^MINIMUM_FRAGMENTS_FOR_MSA=.*/MINIMUM_FRAGMENTS_FOR_MSA=4/' "$BASE/consurfrc.local" > "$BASE/consurfrc.run"

export CONSURFCONF="$BASE/consurfrc.run"

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

echo "WAPL ConSurf runs complete"
