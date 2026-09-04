from pathlib import Path

p = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\input_splits\STAG1\stag1_invertebrates_sanitized.fasta")
records = []
header = None
buf = []
for raw in p.read_text(encoding="utf-8").splitlines():
    line = raw.strip()
    if not line:
        continue
    if line.startswith(">"):
        if header is not None:
            records.append((header, "".join(buf)))
        header = line[1:].strip()
        buf = []
    else:
        buf.append(line)
if header is not None:
    records.append((header, "".join(buf)))

records[0] = ("Ciona_intestinalis", records[0][1])

with p.open("w", encoding="utf-8", newline="\n") as f:
    for h, s in records:
        f.write(f">{h}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i+80] + "\n")

print("updated first header to Ciona_intestinalis")
