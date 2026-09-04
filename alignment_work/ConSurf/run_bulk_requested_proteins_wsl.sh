#!/usr/bin/env bash
set -euo pipefail

SITE=/mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work
BASE=/home/atforpetesaken/consurf_runs/bulk_requested_20260707
CONSURF_BIN=/home/atforpetesaken/consurf_ws/consurf

mkdir -p "$BASE"
cp "$SITE/ConSurf/consurfrc.local.stag1" "$BASE/consurfrc.local"
sed -i 's/^MINIMUM_FRAGMENTS_FOR_MSA=.*/MINIMUM_FRAGMENTS_FOR_MSA=4/' "$BASE/consurfrc.local"
export CONSURFCONF="$BASE/consurfrc.local"

proteins=(PDS5A PDS5B ESCO1 ESCO2 NIPBL WAPL SMC1 SMC3)

for P in "${proteins[@]}"; do
  lower=$(echo "$P" | tr 'A-Z' 'a-z')
  human="Human_${P}"

  in_dir="$SITE/ConSurf/input_splits/$P"
  stage_in="$BASE/input/$P"
  stage_out="$BASE/out/$P"
  stage_work="$BASE/work/$P"

  mkdir -p "$stage_in" "$stage_out/full" "$stage_out/vertebrates" "$stage_out/invertebrates" "$stage_work/full" "$stage_work/vertebrates" "$stage_work/invertebrates"

  cp "$in_dir/${lower}_full_reduced50.fasta" "$stage_in/full.fasta"
  cp "$in_dir/${lower}_vertebrates_sanitized.fasta" "$stage_in/vertebrates.fasta"
  cp "$in_dir/${lower}_invertebrates_sanitized.fasta" "$stage_in/invertebrates.fasta"

  perl "$CONSURF_BIN" -MSA "$stage_in/full.fasta" -SEQ_NAME "$human" -Out_Dir "$stage_out/full" -w "$stage_work/full"
  perl "$CONSURF_BIN" -MSA "$stage_in/vertebrates.fasta" -SEQ_NAME "$human" -Out_Dir "$stage_out/vertebrates" -w "$stage_work/vertebrates"
  perl "$CONSURF_BIN" -MSA "$stage_in/invertebrates.fasta" -SEQ_NAME Ciona_intestinalis -Out_Dir "$stage_out/invertebrates" -w "$stage_work/invertebrates"

  dest_root="$SITE/ConSurf/output/$P"
  mkdir -p "$dest_root/${lower}_consurf_full" "$dest_root/${lower}_consurf_vertebrates" "$dest_root/${lower}_consurf_invertebrates"
  cp -f "$stage_out/full"/* "$dest_root/${lower}_consurf_full/"
  cp -f "$stage_out/vertebrates"/* "$dest_root/${lower}_consurf_vertebrates/"
  cp -f "$stage_out/invertebrates"/* "$dest_root/${lower}_consurf_invertebrates/"

  echo "Completed $P"
done

echo "All requested proteins complete"
