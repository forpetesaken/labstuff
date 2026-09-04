from pathlib import Path

IN_FASTA = Path(r"C:\Users\Nat\Downloads\updated_alignment_0611.fas")
OUT_DIR = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\input_splits")


def read_fasta(path: Path):
    records = []
    header = None
    seq_parts = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip()
            seq_parts = []
        else:
            if header is None:
                raise ValueError("Invalid FASTA")
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8", newline="\n") as out:
        for h, s in records:
            out.write(f">{h}\n")
            for i in range(0, len(s), 80):
                out.write(s[i : i + 80] + "\n")


def sanitize_split(records, split_name):
    out = []
    used = set()
    counter = 1
    for h, s in records:
        if "Homo_sapiens" in h or "Human" in h:
            sid = "Human_CTCF"
        elif "Ciona_intestinalis" in h:
            sid = "Ciona_intestinalis"
        else:
            while True:
                sid = f"{split_name}_{counter:03d}"
                counter += 1
                if sid not in used:
                    break
        if sid in used:
            i = 2
            base = sid
            while f"{base}_{i}" in used:
                i += 1
            sid = f"{base}_{i}"
        used.add(sid)
        out.append((sid, s))
    return out


def main():
    records = read_fasta(IN_FASTA)
    if not records:
        raise ValueError("No records found")

    split_idx = None
    for i, (h, _) in enumerate(records):
        if "Ciona_intestinalis" in h:
            split_idx = i
            break
    if split_idx is None:
        raise ValueError("Could not find Ciona_intestinalis split point")

    vertebrates = records[:split_idx]
    invertebrates = records[split_idx:]
    if not vertebrates or not invertebrates:
        raise ValueError("Empty split produced")

    vertebrates_s = sanitize_split(vertebrates, "VERT")
    invertebrates_s = sanitize_split(invertebrates, "INVERT")

    if not any(h == "Human_CTCF" for h, _ in vertebrates_s):
        raise ValueError("Human_CTCF not present in vertebrate split")
    if not any(h == "Ciona_intestinalis" for h, _ in invertebrates_s):
        raise ValueError("Ciona_intestinalis not present in invertebrate split")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v_path = OUT_DIR / "ctcf_vertebrates_sanitized.fasta"
    i_path = OUT_DIR / "ctcf_invertebrates_sanitized.fasta"
    write_fasta(v_path, vertebrates_s)
    write_fasta(i_path, invertebrates_s)

    summary = OUT_DIR / "ctcf_split_summary.txt"
    summary.write_text(
        "\n".join(
            [
                f"Total: {len(records)}",
                f"Vertebrates: {len(vertebrates_s)}",
                f"Invertebrates: {len(invertebrates_s)}",
                f"Vertebrate query: Human_CTCF",
                f"Invertebrate query: Ciona_intestinalis",
                f"Vertebrate FASTA: {v_path}",
                f"Invertebrate FASTA: {i_path}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {v_path}")
    print(f"Wrote {i_path}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
