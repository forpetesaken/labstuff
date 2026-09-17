from pathlib import Path

IN_FASTA = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\output\PDS5B\pds5b__0708.fas")
OUT_DIR = Path(r"C:\Users\Nat\Downloads\AIDEN Lab\Code\alignment_work\ConSurf\input_splits\PDS5B")

# Alignment is expected to be top=vertebrates, bottom=invertebrates.
# Boundary appears where the invertebrate block starts.
BOUNDARY_MARKERS = [
    "Branchiostoma",
    "Lancelet",
]

# Require a substantive top vertebrate block before accepting a boundary marker.
MIN_PREFIX_RECORDS = 10

# If obvious invertebrate taxa occur above boundary, move them down.
OUTLIER_INVERTEBRATE_MARKERS = [
    "Branchiostoma",
    "Lancelet",
    "Ciona",
    "Ascidian",
    "Tunicate",
    "Leech",
    "Bone_eating_worm",
    "Mint_sauce_worm",
    "dragonfly",
    "grasshopper",
    "coral",
    "Rice",
    "Drosophila",
    "Aedes",
    "mosquito",
    "Strongylocentrotus",
    "earthworm",
    "squid",
    "octopus",
    "chiton",
    "shrimp",
    "tick",
    "Latrodectus",
    "Blood_Fluke",
    "Hydra",
    "sponge",
    "Schizosaccharomyces",
    "Arabidopsis",
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
    fallback = None
    for idx, (h, _) in enumerate(records):
        hl = h.lower()
        if any(marker.lower() in hl for marker in BOUNDARY_MARKERS):
            if fallback is None:
                fallback = idx
            if idx >= MIN_PREFIX_RECORDS:
                return idx
    if fallback is not None:
        return fallback
    raise ValueError("Could not detect PDS5B invertebrate boundary")


def is_invertebrate_outlier(header: str):
    hl = header.lower()
    return any(marker.lower() in hl for marker in OUTLIER_INVERTEBRATE_MARKERS)


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

    full_s = sanitize(records, "ALL", "human", "Human_PDS5B")
    vert_s = sanitize(vertebrates, "VERT", "human", "Human_PDS5B")
    invert_s = sanitize(invertebrates, "INVERT", "ciona", "Ciona_intestinalis")

    write_fasta(OUT_DIR / "pds5b_full_0708_sanitized.fasta", full_s)
    write_fasta(OUT_DIR / "pds5b_vertebrates_sanitized.fasta", vert_s)
    write_fasta(OUT_DIR / "pds5b_invertebrates_sanitized.fasta", invert_s)

    summary = [
        f"Total: {len(records)}",
        f"Boundary index (1-based): {boundary_idx + 1}",
        f"Vertebrates: {len(vertebrates)}",
        f"Invertebrates: {len(invertebrates)}",
        f"Full all: {len(full_s)}",
        "Vertebrate query: Human_PDS5B",
        "Invertebrate query: Ciona_intestinalis",
    ]
    if moved_outliers:
        summary.append(f"Moved invertebrate outliers above boundary: {len(moved_outliers)}")
        summary.extend([f"  - {name}" for name in moved_outliers])

    (OUT_DIR / "pds5b_split_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("Done")


if __name__ == "__main__":
    main()
