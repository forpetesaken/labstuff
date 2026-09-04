from pathlib import Path

IN_FASTA = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output\NIPBL\NIPBL_0708.fas")
OUT_DIR = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\input_splits\NIPBL")
INVERTEBRATE_MARKERS = [
    "Branchiostoma",
    "Ciona",
    "Ascidian",
    "Tunicate",
    "Strongylocentrotus",
    "Acorn_worm",
    "Sea_Cucumber",
    "Common_earthworm",
    "Leech",
    "Bone_eating_worm",
    "Scale_worm",
    "Mint_sauce_worm",
    "squid",
    "octopus",
    "chiton",
    "limpet",
    "Drosophila",
    "Aedes",
    "mosquito",
    "dragonfly",
    "grasshopper",
    "cockroach",
    "silkworm",
    "butterfly",
    "Beetle",
    "Camponotus",
    "springtail",
    "shrimp",
    "tick",
    "Latrodectus",
    "Caenorhabditis",
    "Blood_Fluke",
    "Bloodfluke_planorb",
    "Anemone",
    "Hydra",
    "sponge",
    "Schizosaccharomyces",
    "Arabidopsis",
    "Rice",
    "Sequence_1",
]


def read_fasta(path: Path):
    records = []
    header = None
    buf = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf).upper()))
            header = line[1:].strip()
            buf = []
        else:
            buf.append(line)
    if header is not None:
        records.append((header, "".join(buf).upper()))
    return records


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for h, s in records:
            f.write(f">{h}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")


def sanitize(records, prefix, query_key, query_name):
    out = []
    used = set()
    counter = 1
    for h, s in records:
        if query_key.lower() in h.lower():
            sid = query_name
        elif "Ciona_intestinalis" in h:
            sid = "Ciona_intestinalis"
        else:
            while True:
                sid = f"{prefix}_{counter:03d}"
                counter += 1
                if sid not in used:
                    break
        base = sid
        i = 2
        while sid in used:
            sid = f"{base}_{i}"
            i += 1
        used.add(sid)
        out.append((sid, s))
    return out


def pick_evenly(records, target, must_include_keys):
    must = [(h, s) for h, s in records if any(k.lower() in h.lower() for k in must_include_keys)]
    rest = [(h, s) for h, s in records if not any(k.lower() in h.lower() for k in must_include_keys)]
    n_rest = target - len(must)
    if n_rest <= 0:
        return must
    if len(rest) <= n_rest:
        return must + rest
    step = len(rest) / n_rest
    sampled = [rest[int(i * step)] for i in range(n_rest)]
    return must + sampled


def main():
    records = read_fasta(IN_FASTA)
    vertebrates = []
    invertebrates = []
    for h, s in records:
        if any(marker.lower() in h.lower() for marker in INVERTEBRATE_MARKERS):
            invertebrates.append((h, s))
        else:
            vertebrates.append((h, s))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    full_s = sanitize(records, "ALL", "human", "Human_NIPBL")
    write_fasta(OUT_DIR / "nipbl_full_all80.fasta", full_s)

    vert_s = sanitize(vertebrates, "VERT", "human", "Human_NIPBL")
    write_fasta(OUT_DIR / "nipbl_vertebrates_0708_sanitized.fasta", vert_s)

    invert_s = sanitize(invertebrates, "INVERT", "ciona", "Ciona_intestinalis")
    write_fasta(OUT_DIR / "nipbl_invertebrates_0708_sanitized.fasta", invert_s)

    summary = [
        f"Protein: NIPBL",
        f"Total: {len(records)}",
        f"Vertebrates: {len(vertebrates)}",
        f"Invertebrates: {len(invertebrates)}",
        f"Full all: {len(full_s)}",
        f"Vertebrate query: Human_NIPBL",
        f"Invertebrate query: Ciona_intestinalis",
    ]
    (OUT_DIR / "nipbl_0708_split_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("Done")


if __name__ == "__main__":
    main()
