#!/usr/bin/env bash
set -euo pipefail

RUN=/home/atforpetesaken/consurf_runs/rad21_20260707
IN="$RUN/input"
OUT="$RUN/output"
WORK="$RUN/work"

mkdir -p "$IN" "$OUT" "$WORK"
cp "/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/output/RAD21/updated_RAD21alignment_0625.fas" "$IN/rad21_raw.fasta"

perl -ne 'BEGIN{$i=0;%seen=();} chomp; s/\r$//; if(/^>/){$i++; if($i==1){$id="Human_RAD21";} else {s/^>//; ($id)=split(/\s+/, $_); $id =~ s/[^A-Za-z0-9_.|:-]/_/g; $id = "SEQ_$i" if $id eq "";} my $base=$id; my $k=1; while(exists $seen{$id}){$k++; $id = $base."_".$k;} $seen{$id}=1; print ">$id\n";} else {s/\s+//g; print "$_\n";}' "$IN/rad21_raw.fasta" > "$IN/rad21_sanitized.fasta"

export CONSURFCONF="/mnt/c/Users/Nat/Downloads/AIDEN Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1"
perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$IN/rad21_sanitized.fasta" \
  -SEQ_NAME Human_RAD21 \
  -Out_Dir "$OUT" \
  -w "$WORK"

echo "RAD21 ConSurf run complete"
echo "Output: $OUT"
