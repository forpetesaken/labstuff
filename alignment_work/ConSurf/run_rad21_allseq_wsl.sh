#!/usr/bin/env bash
set -euo pipefail

BASE=/home/atforpetesaken/consurf_runs/rad21_allseq_20260707
mkdir -p "$BASE/input" "$BASE/out_full" "$BASE/work_full"

cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/consurfrc.local.stag1 "$BASE/consurfrc.local"
cp /mnt/c/Users/Nat/Downloads/AIDEN\ Lab/Code/alignment_work/ConSurf/output/RAD21/updated_RAD21alignment_0625.fas "$BASE/input/rad21_raw.fasta"

python3 - <<'PY'
from pathlib import Path
inp = Path('/home/atforpetesaken/consurf_runs/rad21_allseq_20260707/input/rad21_raw.fasta')
out = Path('/home/atforpetesaken/consurf_runs/rad21_allseq_20260707/input/rad21_full_all.fasta')
records=[]; header=None; buf=[]
for raw in inp.read_text(encoding='utf-8').splitlines():
    line=raw.strip()
    if not line:
        continue
    if line.startswith('>'):
        if header is not None:
            records.append((header,''.join(buf).upper()))
        header=line[1:].strip(); buf=[]
    else:
        buf.append(line)
if header is not None:
    records.append((header,''.join(buf).upper()))
used=set(); counter=1
with out.open('w', encoding='utf-8', newline='\n') as f:
    for h,s in records:
        if 'Homo sapiens' in h or 'Homo_sapiens' in h or 'Human' in h:
            sid='Human_RAD21'
        elif 'Ciona_intestinalis' in h:
            sid='Ciona_intestinalis'
        else:
            while True:
                sid=f'ALL_{counter:03d}'
                counter+=1
                if sid not in used:
                    break
        base=sid; i=2
        while sid in used:
            sid=f'{base}_{i}'; i+=1
        used.add(sid)
        f.write(f'>{sid}\n')
        for j in range(0, len(s), 80):
            f.write(s[j:j+80]+'\n')
PY

export CONSURFCONF="$BASE/consurfrc.local"
perl /home/atforpetesaken/consurf_ws/consurf \
  -MSA "$BASE/input/rad21_full_all.fasta" \
  -SEQ_NAME Human_RAD21 \
  -Out_Dir "$BASE/out_full" \
  -w "$BASE/work_full"

echo "RAD21 all-sequence full ConSurf run complete"
