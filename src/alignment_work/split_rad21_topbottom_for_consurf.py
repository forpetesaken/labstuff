from pathlib import Path

IN_FASTA = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output\RAD21\updated_RAD21alignment_0625.fas")
OUT_DIR = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\input_splits\RAD21")

# First invertebrate boundary should start at the lower block in this ordered alignment.
BOUNDARY_MARKERS = [
    "Branchiostoma_lanceolatum",
    "Branchiostoma_lanceolatum_NCBI",
    "Lancelet_",
    "Ciona_intestinalis",
    "Ascidian",
    "Strongylocentrotus",
    "Sea_Cucumber",
    "Common_earthworm",
    "Leech",
    "Drosophila",
    "Aedes",
    "Schizosaccharomyces",
    "Arabidopsis",
    "Rice",
    "Sequence_1",
]

# If any obvious invertebrate appears above the block boundary, move it to invertebrates.
OUTLIER_INVERTEBRATE_MARKERS = [
    "Branchiostoma",
    "Lancelet",
    "Ciona",
    "Ascidian",
    "Strongylocentrotus",
    "Sea_Cucumber",
    "earthworm",
    "Leech",
    "squid",
    "octopus",
    "chiton",
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
    "Blood_Fluke",
    "Bloodfluke",
    "Anemone",
    "coral",
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
                f.write(s[i : i + 80] + "\n")


def find_boundary_index(records):
    for idx, (h, _) in enumerate(records):
        h_lower = h.lower()
        for marker in BOUNDARY_MARKERS:
            m = marker.lower()
            if m == "branchiostoma_lanceolatum":
                if h_lower == m or h_lower == "branchiostoma_lanceolatum_ncbi":
                    return idx
            elif m.endswith("_"):
                if h_lower.startswith(m):
                    return idx
            else:
                if m in h_lower:
                    return idx
    raise ValueError("Could not detect invertebrate boundary in RAD21 alignment")


def is_invertebrate_outlier(header: str):
    h = header.lower()
    return any(marker.lower() in h for marker in OUTLIER_INVERTEBRATE_MARKERS)


def sanitize(records, prefix, query_key, query_name):
    out = []
    used = set()
    counter = 1
    for h, s in records:
        hl = h.lower()
        if query_key.lower() in hl:
            sid = query_name
        elif "ciona_intestinalis" in hl:
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


def main():
    records = read_fasta(IN_FASTA)
    boundary_idx = find_boundary_index(records)

    vertebrates = []
    invertebrates = []
    moved_outliers = []

    for idx, (h, s) in enumerate(records):
        if idx >= boundary_idx:
            invertebrates.append((h, s))
        elif is_invertebrate_outlier(h):
            invertebrates.append((h, s))
            moved_outliers.append(h)
        else:
            vertebrates.append((h, s))

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    vert_s = sanitize(vertebrates, "VERT", "human", "Human_RAD21")
    invert_s = sanitize(invertebrates, "INVERT", "ciona", "Ciona_intestinalis")

    write_fasta(OUT_DIR / "rad21_vertebrates_sanitized.fasta", vert_s)
    write_fasta(OUT_DIR / "rad21_invertebrates_sanitized.fasta", invert_s)

    summary = [
        f"Total: {len(records)}",
        f"Boundary index (1-based): {boundary_idx + 1}",
        f"Vertebrates: {len(vertebrates)}",
        f"Invertebrates: {len(invertebrates)}",
        "Vertebrate query: Human_RAD21",
        "Invertebrate query: Ciona_intestinalis",
    ]

    if moved_outliers:
        summary.append(f"Moved invertebrate outliers above boundary: {len(moved_outliers)}")
        summary.extend([f"  - {name}" for name in moved_outliers])

    (OUT_DIR / "rad21_split_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("Done")


if __name__ == "__main__":
    main()
