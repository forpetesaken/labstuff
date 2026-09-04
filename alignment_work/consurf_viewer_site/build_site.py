import argparse
import json
import re
from pathlib import Path
from urllib.request import urlopen

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"


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
                records.append((header, "".join(buf)))
            header = line[1:].strip()
            buf = []
        else:
            if header is None:
                raise ValueError(f"Invalid FASTA: {path}")
            buf.append(line)
    if header is not None:
        records.append((header, "".join(buf)))
    return records


def parse_grades(grades_path: Path):
    rows = []
    with grades_path.open("r", encoding="utf-8") as handle:
        in_table = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not in_table:
                if line.strip().startswith("POS"):
                    in_table = True
                continue

            if not line.strip():
                continue
            if line.lstrip().startswith("*"):
                break

            parts = line.split()
            if len(parts) < 4:
                continue
            if not parts[0].isdigit():
                continue

            pos = int(parts[0])
            aa = parts[1]
            score = float(parts[2])
            color_raw = parts[3]
            grade = int(color_raw.replace("*", ""))
            low_conf = "*" in color_raw
            rows.append(
                {
                    "pos": pos,
                    "aa": aa,
                    "score": score,
                    "grade": grade,
                    "low_conf": low_conf,
                }
            )

    if not rows:
        raise ValueError(f"No rows parsed from {grades_path}")

    rows.sort(key=lambda r: r["pos"])
    return rows


def query_name_from_grades_path(grades_path: Path):
    return grades_path.stem.removesuffix("_consurf")


def extract_species_label(header: str):
  os_match = re.search(r"\bOS=([^=]+?)(?:\sOX=|\sGN=|\sPE=|\sSV=|$)", header)
  if os_match:
    return os_match.group(1).strip()
  return header.split()[0].strip()


def build_display_name_lookup(source_records):
    lookup = {}
    for header, seq in source_records:
        key = seq.upper()
        lookup.setdefault(key, []).append(extract_species_label(header))
    return lookup


def parse_msa_for_viewer(msa_path: Path, query_name: str, display_source_path: Path | None = None):
    records = read_fasta(msa_path)
    if not records:
        return None

    display_lookup = None
    if display_source_path and display_source_path.exists():
        display_lookup = build_display_name_lookup(read_fasta(display_source_path))

    query_seq = None
    viewer_records = []
    for header, seq in records:
        if header == query_name:
            query_seq = seq
        display_name = header
        if display_lookup:
            names = display_lookup.get(seq.upper())
            if names:
                display_name = names.pop(0)
        viewer_records.append({"name": display_name, "seq": seq})

    if query_seq is None:
        return None

    pos_to_col = {}
    query_pos = 0
    for col_idx, aa in enumerate(query_seq):
        if aa != "-":
            query_pos += 1
            pos_to_col[query_pos] = col_idx

    return {
        "query_name": query_name,
        "aligned_length": len(query_seq),
        "pos_to_col": pos_to_col,
      "records": viewer_records,
    }


def make_plotly_script(embed_plotly: bool):
    if not embed_plotly:
        return f'<script src="{PLOTLY_CDN}"></script>'

    with urlopen(PLOTLY_CDN, timeout=30) as response:
        plotly_js = response.read().decode("utf-8")
    return f"<script>{plotly_js}</script>"


def compute_human_mapping_from_source(source_path: Path, source_query_key: str):
    records = read_fasta(source_path)
    human_seq = None
    query_seq = None
    for header, seq in records:
        lower = header.lower()
        if human_seq is None and ("homo_sapiens" in lower or "human" in lower):
            human_seq = seq
        if query_seq is None and source_query_key.lower() in lower:
            query_seq = seq

    if human_seq is None or query_seq is None:
        return None
    if len(human_seq) != len(query_seq):
        return None

    human_index = 0
    mapping = []
    for human_char, query_char in zip(human_seq, query_seq):
        if human_char != "-":
            human_index += 1
        if query_char != "-":
            if human_char != "-":
                mapping.append(
                    {
                        "present": 1,
                        "aa": human_char,
                        "count": human_index,
                    }
                )
            else:
                mapping.append(
                    {
                        "present": None,
                        "aa": None,
                        "count": None,
                    }
                )
    return mapping


def gather_data(project_root: Path):
    cfg = {
        "RAD21": {
            "Full": project_root / "ConSurf/output/RAD21/rad21_consurf_full/Human_RAD21_consurf.grades",
            "Vertebrates": project_root / "ConSurf/output/RAD21/rad21_consurf_vertebrates/Human_RAD21_consurf.grades",
            "Invertebrates": project_root / "ConSurf/output/RAD21/rad21_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
            "_source": project_root / "ConSurf/output/RAD21/updated_RAD21alignment_0625.fas",
            "_source_query": "Ciona_intestinalis",
            "_display_source": project_root / "ConSurf/output/RAD21/updated_RAD21alignment_0625.fas",
        },
        "STAG1": {
            "Full": project_root / "ConSurf/output/STAG1/stag1_consurf_full/Human_STAG1_consurf.grades",
            "Vertebrates": project_root / "ConSurf/output/STAG1/stag1_consurf_vertebrates/Human_STAG1_consurf.grades",
            "Invertebrates": project_root / "ConSurf/output/STAG1/stag1_consurf_invertebrates/Branchiostoma_lanceolatum_consurf.grades",
            "_source": Path(r"C:/Users/Nat/Downloads/Bioinformatics/updated_STAG1alignment_0612.fas"),
            "_source_query": "Branchiostoma_lanceolatum",
            "_display_source": Path(r"C:/Users/Nat/Downloads/Bioinformatics/updated_STAG1alignment_0612.fas"),
        },
        "STAG2": {
            "Full": project_root / "ConSurf/output/STAG2/stag2_consurf_full/Human_STAG2_consurf.grades",
            "Vertebrates": project_root / "ConSurf/output/STAG2/stag2_consurf_vertebrates/Human_STAG2_consurf.grades",
            "Invertebrates": project_root / "ConSurf/output/STAG2/stag2_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
            "_source": project_root / "alignments_out/STAG2/01_STAG2_aligned.fasta",
            "_source_query": "Ciona_intestinalis",
            "_display_source": project_root / "alignments_out/STAG2/01_STAG2_aligned.fasta",
        },
        "CTCF": {
            "Full": project_root / "ConSurf/output/CTCF/ctcf_consurf_run/Human_CTCF_consurf.grades",
            "Vertebrates": project_root / "ConSurf/output/CTCF/ctcf_consurf_vertebrates/Human_CTCF_consurf.grades",
            "Invertebrates": project_root / "ConSurf/output/CTCF/ctcf_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
            "_source": project_root / "alignments_out/CTCF/01_CTCF_aligned.fasta",
            "_source_query": "Ciona_intestinalis",
            "_display_source": Path(r"C:/Users/Nat/Downloads/Bioinformatics/updated_alignment_0611.fas"),
        },
        "WAPL": {
            "Full": project_root / "ConSurf/output/WAPL/wapl_consurf_full/Human_WAPL_consurf.grades",
            "Vertebrates": project_root / "ConSurf/output/WAPL/wapl_consurf_vertebrates/Human_WAPL_consurf.grades",
            "Invertebrates": project_root / "ConSurf/output/WAPL/wapl_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
            "_source": project_root / "ConSurf/output/WAPL/WAPLalignment_0708.fas",
            "_source_query": "Ciona_intestinalis",
            "_display_source": project_root / "alignments_out/WAPL/01_WAPL_aligned.fasta",
        },
        "PDS5A": {
          "Full": project_root / "ConSurf/output/PDS5A/pds5a_consurf_full/Human_PDS5A_consurf.grades",
          "Vertebrates": project_root / "ConSurf/output/PDS5A/pds5a_consurf_vertebrates/Human_PDS5A_consurf.grades",
          "Invertebrates": project_root / "ConSurf/output/PDS5A/pds5a_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
          "_source": project_root / "ConSurf/output/PDS5A/psd5A_0708.fas",
          "_source_query": "Ciona_intestinalis",
          "_display_source": project_root / "ConSurf/output/PDS5A/psd5A_0708.fas",
        },
        "PDS5B": {
          "Full": project_root / "ConSurf/output/PDS5B/pds5b_consurf_full/Human_PDS5B_consurf.grades",
          "Vertebrates": project_root / "ConSurf/output/PDS5B/pds5b_consurf_vertebrates/Human_PDS5B_consurf.grades",
          "Invertebrates": project_root / "ConSurf/output/PDS5B/pds5b_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
          "_source": project_root / "ConSurf/output/PDS5B/pds5b__0708.fas",
          "_source_query": "Ciona_intestinalis",
          "_display_source": project_root / "ConSurf/output/PDS5B/pds5b__0708.fas",
        },
        "SMC1": {
          "Full": project_root / "ConSurf/output/SMC1/smc1_consurf_full/Human_SMC1_consurf.grades",
          "Vertebrates": project_root / "ConSurf/output/SMC1/smc1_consurf_vertebrates/Human_SMC1_consurf.grades",
          "Invertebrates": project_root / "ConSurf/output/SMC1/smc1_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
          "_source": project_root / "ConSurf/output/SMC1/SMC1_0708.fas",
          "_source_query": "Ciona_intestinalis",
          "_display_source": project_root / "alignments_out/SMC1/01_SMC1_aligned.fasta",
        },
        "SMC3": {
          "Full": project_root / "ConSurf/output/SMC3/smc3_consurf_full/Human_SMC3_consurf.grades",
          "Vertebrates": project_root / "ConSurf/output/SMC3/smc3_consurf_vertebrates/Human_SMC3_consurf.grades",
          "Invertebrates": project_root / "ConSurf/output/SMC3/smc3_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
          "_source": project_root / "ConSurf/output/SMC3/smc3_0708.fas",
          "_source_query": "Ciona_intestinalis",
          "_display_source": project_root / "ConSurf/output/SMC3/smc3_0708.fas",
        },
        "NIPBL": {
          "Full": project_root / "ConSurf/output/NIPBL/nipbl_consurf_full/Human_NIPBL_consurf.grades",
          "Vertebrates": project_root / "ConSurf/output/NIPBL/nipbl_consurf_vertebrates/Human_NIPBL_consurf.grades",
          "Invertebrates": project_root / "ConSurf/output/NIPBL/nipbl_consurf_invertebrates/Ciona_intestinalis_consurf.grades",
          "_source": project_root / "ConSurf/output/NIPBL/NIPBL_0708.fas",
          "_source_query": "Ciona_intestinalis",
          "_display_source": project_root / "ConSurf/output/NIPBL/NIPBL_0708.fas",
        },
    }

    payload = {}
    for protein, datasets in cfg.items():
        payload[protein] = {}
        source_path = datasets.get("_source")
        source_query_key = datasets.get("_source_query")
        display_source_path = datasets.get("_display_source")
        invertebrate_mapping = None
        if source_path:
            if not source_path.exists():
                raise FileNotFoundError(f"Missing expected source alignment: {source_path}")
            if source_query_key:
                invertebrate_mapping = compute_human_mapping_from_source(source_path, source_query_key)

        for dataset_name, path in datasets.items():
            if dataset_name.startswith("_"):
                continue
            if not path.exists():
                raise FileNotFoundError(f"Missing expected grades file: {path}")
            rows = parse_grades(path)
            query_name = query_name_from_grades_path(path)
            msa_path = path.with_name(f"{query_name}_msa_file.fas")
            msa_data = None
            if msa_path.exists():
              msa_data = parse_msa_for_viewer(msa_path, query_name, display_source_path)
            human_presence = None
            if dataset_name == "Full":
                human_presence = [
                    {"present": 1, "aa": row["aa"], "count": row["pos"]}
                    for row in rows
                ]
            elif dataset_name == "Invertebrates" and invertebrate_mapping is not None:
                if len(invertebrate_mapping) == len(rows):
                    human_presence = invertebrate_mapping

            payload[protein][dataset_name] = {
                "rows": rows,
                "human_presence": human_presence,
                "msa": msa_data,
              "n_sequences": len(msa_data["records"]) if msa_data else None,
              "msa_source": msa_path.name if msa_data else None,
            }
    return payload


def build_html(payload, plotly_script_tag: str):
    datasets_json = json.dumps(payload)
    overview_json = json.dumps(
        {
      "CTCF": {
        "summary": "CTCF is a chromatin architectural factor that positions loop boundaries and recruits cohesin to specific DNA elements.",
        "sections": [
          "KTYQR motif 39-43 lies in the N-terminal regulatory region associated with chromatin-linked CTCF activity.",
          "YDF motif 226-228 is positioned near the transition into the DNA-binding core and is useful as a conserved anchor feature.",
          "Lysine-rich linker segment 214-224 (KKTKKTKKSKL) provides a basic connector between N-terminal and zinc-finger regions.",
          "ZF array 266-577 mediates sequence-specific DNA binding that underlies boundary/insulator function.",
        ],
        "citations": [
          {"label": "UniProt P49711", "url": "https://www.uniprot.org/uniprotkb/P49711/entry"},
          {"label": "PMID 8649389", "url": "https://pubmed.ncbi.nlm.nih.gov/8649389/"},
          {"label": "PMID 17827499", "url": "https://pubmed.ncbi.nlm.nih.gov/17827499/"},
          {"label": "PMID 18550811", "url": "https://pubmed.ncbi.nlm.nih.gov/18550811/"},
        ],
      },
      "RAD21": {
        "summary": "RAD21 is the kleisin bridge between SMC1A and SMC3 and its cleavage and interaction interfaces are central to cohesion release and loader coupling.",
        "sections": [
          "Separase cleavage sites around R172 and R450 are required for anaphase chromatid separation.",
          "Caspase-sensitive segment 276-280 (including D279) controls apoptotic RAD21 processing.",
          "Region 154-171 contacts NIPBL; region 287-403 contributes WAPL/PDS5B/STAG1 interactions.",
        ],
        "citations": [
          {"label": "UniProt O60216", "url": "https://www.uniprot.org/uniprotkb/O60216/entry"},
          {"label": "PMID 11509732", "url": "https://pubmed.ncbi.nlm.nih.gov/11509732/"},
          {"label": "PMID 12417729", "url": "https://pubmed.ncbi.nlm.nih.gov/12417729/"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
      "STAG1": {
        "summary": "STAG1 is a cohesin SA subunit whose SCD/HEAT architecture supports cohesin assembly and chromatin interactions.",
        "sections": [
          "SCD domain 296-381 is a conserved STAG family module in the cohesin-associated scaffold.",
          "HEAT-repeat arm (approximately 86-915) forms a structural platform for cohesin-partner interactions.",
          "Regulatory phosphosites around S1062 and S1093 mark a C-terminal control region linked to cohesin cycle timing.",
        ],
        "citations": [
          {"label": "UniProt Q8WVM7", "url": "https://www.uniprot.org/uniprotkb/Q8WVM7/entry"},
          {"label": "PMID 11076961", "url": "https://pubmed.ncbi.nlm.nih.gov/11076961/"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
      "STAG2": {
        "summary": "STAG2 is the alternative SA subunit in many cohesin complexes and uses SCD/HEAT features to stabilize cohesin architecture.",
        "sections": [
          "Residue N327 in the SCD region is function-sensitive for RAD21/cohesin interaction.",
          "SCD domain 293-378 is a conserved STAG module implicated in cohesin architecture.",
          "C-terminal regulatory phosphosites around S1058, S1061 and T1112 mark a control region during cohesin cycling.",
        ],
        "citations": [
          {"label": "UniProt Q8N3U4", "url": "https://www.uniprot.org/uniprotkb/Q8N3U4/entry"},
          {"label": "PMID 11076961", "url": "https://pubmed.ncbi.nlm.nih.gov/11076961/"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
      "WAPL": {
        "summary": "WAPL promotes cohesin release from DNA, and short FGF-containing motifs are key determinants of interactions in the release pathway.",
        "sections": [
          "FGF motifs at 73-75, 429-431 and 453-455 are annotated interaction hotspots.",
          "Mutagenesis of the 429-431 and 453-455 FGF motifs reduces or abolishes PDS5B-dependent interactions.",
          "These motifs map to the WAPL-PDS5 regulatory interface that controls cohesin unloading timing.",
        ],
        "citations": [
          {"label": "UniProt Q7Z5K2", "url": "https://www.uniprot.org/uniprotkb/Q7Z5K2/entry"},
          {"label": "PMID 19696148", "url": "https://pubmed.ncbi.nlm.nih.gov/19696148/"},
        ],
      },
      "PDS5A": {
        "summary": "PDS5A is a cohesin regulatory factor that binds cohesin and WAPL to tune residence time and release dynamics on chromatin.",
        "sections": [
          "External residue-level evidence (UniProt Q29RF7) includes PTM sites S1097 and K1146 in the C-terminal regulatory region.",
          "A compact PTM cluster is annotated at S1195, T1208 and K1211.",
          "Additional externally curated PTM residues include K1290 and S1305.",
        ],
        "citations": [
          {"label": "UniProt Q29RF7", "url": "https://www.uniprot.org/uniprotkb/Q29RF7/entry"},
        ],
      },
      "PDS5B": {
        "summary": "PDS5B is a cohesin regulatory subunit that partners with cohesin and WAPL to influence cohesin retention and release from chromatin.",
        "sections": [
          "PDS5B contains HEAT-repeat architecture that supports broad protein interaction surfaces.",
          "Comparing full, vertebrate, and invertebrate ConSurf runs helps localize regions under strong evolutionary constraint.",
          "This tab reports residue-level conservation from the 0708 PDS5B alignment split into vertebrate and invertebrate sets.",
        ],
        "citations": [
          {"label": "UniProt Q9NTI5", "url": "https://www.uniprot.org/uniprotkb/Q9NTI5/entry"},
        ],
      },
      "NIPBL": {
        "summary": "NIPBL is the core cohesin loader and uses specific motifs and HEAT-repeat regions to engage cohesin and chromatin.",
        "sections": [
          "PxVxL motif around 996-1009 contributes partner recognition in chromatin regulatory contexts.",
          "V1003 and L1005 within this motif are mutagenesis-sensitive residues that weaken reported binding interactions.",
          "HEAT-repeat region is directly implicated in contacts with cohesin subunits during loading complexes.",
        ],
        "citations": [
          {"label": "UniProt Q6KC79", "url": "https://www.uniprot.org/uniprotkb/Q6KC79/entry"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
      "SMC1": {
        "summary": "SMC1A is a cohesin ATPase subunit whose ATP-binding and regulatory phosphosites are linked to checkpoint and cohesion control.",
        "sections": [
          "ATP-binding motif 32-39 marks the N-terminal head region central to ATPase cycling.",
          "Hinge domain 515-629 forms the SMC1A-SMC3 dimerization interface critical for ring architecture.",
          "Checkpoint-regulated phosphosites S957 and S966 are mutagenesis-sensitive in DNA damage response assays.",
        ],
        "citations": [
          {"label": "UniProt Q14683", "url": "https://www.uniprot.org/uniprotkb/Q14683/entry"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
      "SMC3": {
        "summary": "SMC3 partners with SMC1A in cohesin ATPase cycling; specific lysines and head/hinge regions are key to stable cohesion and loop control.",
        "sections": [
          "ATP-binding motif 32-39 marks the catalytic head region used in ATP-driven conformational transitions.",
          "Acetylation-sensitive K105 and K106 are classic cohesion-control residues with strong functional evidence.",
          "Hinge region supports SMC1A-SMC3 dimer architecture and thus constrains ring mechanics during extrusion.",
        ],
        "citations": [
          {"label": "UniProt Q9UQE7", "url": "https://www.uniprot.org/uniprotkb/Q9UQE7/entry"},
          {"label": "PMID 32409525", "url": "https://pubmed.ncbi.nlm.nih.gov/32409525/"},
        ],
      },
        }
    )
    return f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>ConSurf Interactive Viewer</title>
  {plotly_script_tag}
  <style>
    body {{
      font-family: Segoe UI, Arial, sans-serif;
      margin: 0;
      background: #f7f9fc;
      color: #1f2937;
    }}
    .wrap {{
      width: min(98vw, 2200px);
      max-width: 2200px;
      margin: 12px auto;
      padding: 0 12px;
    }}
    .tabs {{
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      margin-bottom: 12px;
    }}
    .overview {{
      border: 1px solid #cbd5e1;
      border-radius: 10px;
      background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%);
      padding: 10px 12px;
      margin-bottom: 10px;
      color: #0f172a;
      font-size: 13px;
      line-height: 1.45;
    }}
    .overview p {{
      margin: 0 0 7px 0;
    }}
    .overview-callout {{
      border-left: 4px solid #2563eb;
      background: #eff6ff;
      padding: 7px 9px;
      border-radius: 6px;
      margin: 0 0 7px 0;
    }}
    .overview-citations {{
      font-size: 12px;
      color: #334155;
      margin: 0;
    }}
    .overview a {{
      color: #1d4ed8;
      text-decoration: none;
    }}
    .overview a:hover {{
      text-decoration: underline;
    }}
    .tab {{
      border: 1px solid #cbd5e1;
      background: #fff;
      border-radius: 8px;
      padding: 7px 12px;
      cursor: pointer;
      font-weight: 600;
      color: #334155;
    }}
    .tab.active {{
      background: #1d4ed8;
      border-color: #1d4ed8;
      color: white;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 10px;
      margin-bottom: 10px;
    }}
    select, input, button {{
      font-size: 14px;
      padding: 6px 10px;
      border: 1px solid #cbd5e1;
      border-radius: 6px;
      background: #fff;
    }}
    button {{
      cursor: pointer;
      background: #eef2ff;
      border-color: #c7d2fe;
    }}
    button:hover {{
      background: #e0e7ff;
    }}
    .hint {{
      font-size: 13px;
      color: #475569;
    }}
    .status {{
      font-size: 12px;
      color: #334155;
      min-height: 18px;
      margin-bottom: 8px;
    }}
    .dataset-meta {{
      font-size: 12px;
      color: #475569;
      margin-bottom: 8px;
    }}
    .panel {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) 280px;
      gap: 12px;
    }}
    #plot {{
      width: 100%;
      height: 82vh;
      min-height: 740px;
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      background: #fff;
    }}
    .highlights {{
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      background: #fff;
      padding: 10px;
      max-height: 82vh;
      min-height: 740px;
      overflow: auto;
    }}
    .highlights h3 {{
      margin: 2px 0 8px 0;
      font-size: 15px;
    }}
    .highlights ul {{
      list-style: none;
      padding: 0;
      margin: 0;
      display: flex;
      flex-direction: column;
      gap: 6px;
    }}
    .highlights li {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 8px;
      border: 1px solid #e2e8f0;
      border-radius: 6px;
      padding: 6px 8px;
      font-size: 12px;
      background: #f8fafc;
      cursor: pointer;
    }}
    .highlights li.selected {{
      border-color: #2563eb;
      background: #dbeafe;
    }}
    .chip {{
      width: 10px;
      height: 10px;
      border-radius: 999px;
      display: inline-block;
      margin-right: 6px;
      border: 1px solid #94a3b8;
    }}
    .empty {{
      color: #64748b;
      font-size: 12px;
    }}
    .msa-overlay {{
      position: fixed;
      inset: 0;
      background: rgba(15, 23, 42, 0.62);
      display: none;
      align-items: center;
      justify-content: center;
      z-index: 9999;
      padding: 20px;
    }}
    .msa-overlay.open {{
      display: flex;
    }}
    .msa-modal {{
      width: min(92vw, 1600px);
      height: min(88vh, 980px);
      background: linear-gradient(180deg, #fff 0%, #f8fafc 100%);
      border: 1px solid #cbd5e1;
      border-radius: 18px;
      box-shadow: 0 24px 80px rgba(15, 23, 42, 0.28);
      display: grid;
      grid-template-rows: auto auto 1fr;
      overflow: hidden;
    }}
    .msa-head {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 12px;
      padding: 14px 18px 10px 18px;
      border-bottom: 1px solid #e2e8f0;
      background: rgba(255,255,255,0.85);
      backdrop-filter: blur(8px);
    }}
    .msa-title {{
      font-size: 18px;
      font-weight: 700;
      color: #0f172a;
    }}
    .msa-close {{
      border-radius: 999px;
      padding: 8px 12px;
      background: #e2e8f0;
      border-color: #cbd5e1;
    }}
    .msa-meta {{
      font-size: 12px;
      color: #475569;
      padding: 10px 18px;
      border-bottom: 1px solid #e2e8f0;
      background: rgba(248,250,252,0.9);
    }}
    .msa-scroll {{
      overflow: auto;
      padding: 14px 18px 18px 18px;
      background:
        radial-gradient(circle at top left, rgba(191,219,254,0.18), transparent 28%),
        radial-gradient(circle at bottom right, rgba(254,202,202,0.16), transparent 24%),
        #ffffff;
    }}
    .msa-chunk {{
      margin-bottom: 10px;
    }}
    .msa-chunk-head {{
      font-family: Consolas, "Courier New", monospace;
      font-size: 11px;
      color: #475569;
      margin-bottom: 4px;
    }}
    .msa-row {{
      display: grid;
      grid-template-columns: 240px 1fr;
      gap: 12px;
      align-items: start;
      margin-bottom: 0;
      font-family: Consolas, "Courier New", monospace;
      font-size: 12px;
      line-height: 1;
    }}
    .msa-name {{
      position: sticky;
      left: 0;
      background: rgba(255, 255, 255, 0.96);
      color: #0f172a;
      border: 1px solid #e2e8f0;
      padding: 2px 6px;
      border-radius: 6px;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }}
    .msa-seq {{
      display: flex;
      flex-wrap: wrap;
      gap: 0;
      align-items: flex-start;
      align-content: flex-start;
    }}
    .aa-box {{
      width: 16px;
      height: 16px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      border-radius: 0;
      font-size: 11px;
      font-weight: 800;
      line-height: 1;
      box-shadow: inset 0 0 0 1px rgba(15, 23, 42, 0.04);
      color: #0f172a;
      background: #e2e8f0;
    }}
    .aa-gap {{ background: #ffffff; color: #cbd5e1; }}
    .aa-hydrophobic {{ background: #fde68a; color: #92400e; }}
    .aa-polar {{ background: #bbf7d0; color: #166534; }}
    .aa-positive {{ background: #bfdbfe; color: #1d4ed8; }}
    .aa-negative {{ background: #fecaca; color: #991b1b; }}
    .aa-gly {{ background: #ddd6fe; color: #6d28d9; }}
    .aa-pro {{ background: #fbcfe8; color: #be185d; }}
    .aa-cys {{ background: #a5f3fc; color: #155e75; }}
    .aa-aromatic {{ background: #fed7aa; color: #9a3412; }}
    @media (max-width: 1000px) {{
      .panel {{
        grid-template-columns: 1fr;
      }}
      .highlights {{
        max-height: none;
      }}
      .msa-modal {{
        width: 96vw;
        height: 92vh;
      }}
      .msa-row {{
        grid-template-columns: 1fr;
      }}
      .msa-name {{
        position: static;
      }}
    }}
  </style>
</head>
<body>
  <div class=\"wrap\">
    <h2>ConSurf Interactive Viewer</h2>
    <div class=\"tabs\" id=\"protein-tabs\"></div>
    <div class=\"overview\" id=\"protein-overview\"></div>

    <div class=\"toolbar\">
      <label for=\"dataset\">Dataset:</label>
      <select id=\"dataset\"></select>
      <input id="hl-start" type="number" min="1" step="1" placeholder="Start residue" style="width:130px;" />
      <input id="hl-end" type="number" min="1" step="1" placeholder="End residue" style="width:130px;" />
      <button id="hl-add" type="button">Add by residues</button>
      <button id="hl-template" type="button">Load templates</button>
      <div class=\"hint\">Hover points for exact residue number, amino acid, score, and grade.</div>
    </div>


    <div class="toolbar">
      <input id="seq-query" type="text" placeholder="Find motif in current reference sequence (e.g. FGF)" />
      <button id="seq-find" type="button">Find sequence</button>
    </div>
    <div class=\"toolbar\">
      <input id=\"hl-label\" type=\"text\" placeholder=\"Label\" value=\"Region\" />
      <input id=\"hl-color\" type=\"color\" value=\"#f59e0b\" />
      <button id=\"hl-pick\" type=\"button\">Pick from plot</button>
      <button id=\"hl-clear\" type=\"button\">Clear all</button>
    </div>

    <div class=\"status\" id=\"status\"></div>
    <div class=\"dataset-meta\" id=\"dataset-meta\"></div>

    <div class=\"panel\">
      <div id=\"plot\"></div>
      <div class=\"highlights\">
        <h3>Highlights</h3>
        <ul id=\"highlight-list\"></ul>
        <div class="empty" style="margin-top:12px;">Click a highlight to open the colored MSA slice viewer.</div>
      </div>
    </div>
  </div>

  <div class="msa-overlay" id="msa-overlay">
    <div class="msa-modal">
      <div class="msa-head">
        <div class="msa-title">MSA Slice Viewer</div>
        <button class="msa-close" id="msa-close" type="button">Close</button>
      </div>
      <div class="msa-meta" id="msa-meta">Select or create a highlight to view the aligned region.</div>
      <div class="msa-scroll" id="msa-view"></div>
    </div>
  </div>

  <script>
    const proteinDatasets = {datasets_json};

    const proteinTabs = document.getElementById('protein-tabs');
    const proteinOverviewEl = document.getElementById('protein-overview');
    const datasetSelect = document.getElementById('dataset');
    const plotEl = document.getElementById('plot');
    const listEl = document.getElementById('highlight-list');
    const statusEl = document.getElementById('status');
    const datasetMetaEl = document.getElementById('dataset-meta');
    const msaOverlayEl = document.getElementById('msa-overlay');
    const msaMetaEl = document.getElementById('msa-meta');
    const msaViewEl = document.getElementById('msa-view');
    const msaCloseBtn = document.getElementById('msa-close');
    const pickBtn = document.getElementById('hl-pick');
    const clearBtn = document.getElementById('hl-clear');
    const addBtn = document.getElementById('hl-add');
    const templateBtn = document.getElementById('hl-template');
    const findBtn = document.getElementById('seq-find');
    const labelInput = document.getElementById('hl-label');
    const colorInput = document.getElementById('hl-color');
    const startInput = document.getElementById('hl-start');
    const endInput = document.getElementById('hl-end');
    const seqQueryInput = document.getElementById('seq-query');

    const proteins = Object.keys(proteinDatasets);
    const proteinOverviews = {overview_json};
    let currentProtein = proteins[0];
    let currentDataset = Object.keys(proteinDatasets[currentProtein])[0];
    let selectedHighlightId = null;

    const highlights = {{}};
    proteins.forEach((protein) => {{
      highlights[protein] = {{}};
      Object.keys(proteinDatasets[protein]).forEach((dataset) => {{
        highlights[protein][dataset] = [];
      }});
    }});

    let pickMode = false;
    let pickStart = null;

    function setStatus(message) {{
      statusEl.textContent = message;
    }}

    function makeHighlightId() {{
      return `hl_${{Date.now()}}_${{Math.floor(Math.random() * 10000)}}`;
    }}

    function normalizeRange(a, b) {{
      return [Math.min(a, b), Math.max(a, b)];
    }}

    function removeSearchHighlights() {{
      highlights[currentProtein][currentDataset] = highlights[currentProtein][currentDataset].filter((row) => !row.isSearch);
    }}

    function getCurrentHighlights() {{
      return highlights[currentProtein][currentDataset];
    }}

    function getSelectedHighlight() {{
      return getCurrentHighlights().find((row) => row.id === selectedHighlightId) || null;
    }}

    function setSelectedHighlight(id) {{
      selectedHighlightId = id;
      renderHighlightList();
      openMsaModal();
    }}

    function closeMsaModal() {{
      msaOverlayEl.classList.remove('open');
    }}

    function openMsaModal() {{
      renderMsaSlice();
      msaOverlayEl.classList.add('open');
    }}

    function escapeHtml(value) {{
      return String(value)
        .replaceAll('&', '&amp;')
        .replaceAll('<', '&lt;')
        .replaceAll('>', '&gt;')
        .replaceAll('"', '&quot;')
        .replaceAll("'", '&#39;');
    }}

    function aaClass(aa) {{
      if (aa === '-') return 'aa-gap';
      if ('AILMV'.includes(aa)) return 'aa-hydrophobic';
      if ('STNQH'.includes(aa)) return 'aa-polar';
      if ('KR'.includes(aa)) return 'aa-positive';
      if ('DE'.includes(aa)) return 'aa-negative';
      if (aa === 'G') return 'aa-gly';
      if (aa === 'P') return 'aa-pro';
      if (aa === 'C') return 'aa-cys';
      if ('FWY'.includes(aa)) return 'aa-aromatic';
      return '';
    }}

    function renderColoredSeq(seq) {{
      return seq
        .split('')
        .map((aa) => `<span class="aa-box ${{aaClass(aa)}}">${{escapeHtml(aa)}}</span>`)
        .join('');
    }}

    function renderCitations(citations) {{
      if (!Array.isArray(citations) || !citations.length) return 'No citations available.';
      return citations
        .map((c) => `<a href="${{escapeHtml(c.url)}}" target="_blank" rel="noopener noreferrer">${{escapeHtml(c.label)}}</a>`)
        .join('; ');
    }}

    function renderTabs() {{
      proteinTabs.innerHTML = '';
      proteins.forEach((protein) => {{
        const btn = document.createElement('button');
        btn.type = 'button';
        btn.className = 'tab' + (protein === currentProtein ? ' active' : '');
        btn.textContent = protein;
        btn.addEventListener('click', () => {{
          currentProtein = protein;
          currentDataset = Object.keys(proteinDatasets[currentProtein])[0];
          renderAll();
        }});
        proteinTabs.appendChild(btn);
      }});
    }}

    function renderOverview() {{
      const ov = proteinOverviews[currentProtein];
      if (!ov) {{
        proteinOverviewEl.innerHTML = '<p><strong>Overview:</strong> No overview available.</p>';
        return;
      }}

      const sectionParts = Array.isArray(ov.sections) ? ov.sections : [];
      const sectionsHtml = sectionParts.length
        ? '<ul style="margin: 6px 0 0 18px; padding: 0;">' +
          sectionParts.map((text) => `<li>${{escapeHtml(text)}}</li>`).join('') +
          '</ul>'
        : 'No residue-level external annotations are currently listed.';

      proteinOverviewEl.innerHTML =
        `<p><strong>Overview:</strong> ${{escapeHtml(ov.summary || '')}}</p>` +
        `<div class="overview-callout"><strong>DNA extrusion-relevant residues/regions (external evidence):</strong> ${{sectionsHtml}}</div>` +
        `<p class="overview-citations"><strong>Citations:</strong> ${{renderCitations(ov.citations)}}</p>`;
    }}

    function renderDatasetSelect() {{
      datasetSelect.innerHTML = '';
      Object.keys(proteinDatasets[currentProtein]).forEach((dataset) => {{
        const opt = document.createElement('option');
        opt.value = dataset;
        opt.textContent = dataset;
        datasetSelect.appendChild(opt);
      }});
      datasetSelect.value = currentDataset;
    }}

    function renderHighlightList() {{
      const rows = getCurrentHighlights();
      listEl.innerHTML = '';
      if (!rows.length) {{
        const empty = document.createElement('div');
        empty.className = 'empty';
        empty.textContent = 'No highlights yet.';
        listEl.appendChild(empty);
        selectedHighlightId = null;
        closeMsaModal();
        return;
      }}

      if (!rows.some((row) => row.id === selectedHighlightId)) {{
        selectedHighlightId = rows[rows.length - 1].id;
      }}

      rows.forEach((row) => {{
        const li = document.createElement('li');
        if (row.id === selectedHighlightId) {{
          li.classList.add('selected');
        }}
        li.addEventListener('click', () => setSelectedHighlight(row.id));
        const left = document.createElement('div');

        const chip = document.createElement('span');
        chip.className = 'chip';
        chip.style.background = row.color;
        left.appendChild(chip);

        const text = document.createElement('span');
        text.textContent = `${{row.label}}: ${{row.start}}-${{row.end}}`;
        left.appendChild(text);

        const del = document.createElement('button');
        del.type = 'button';
        del.textContent = 'Delete';
        del.addEventListener('click', (ev) => {{
          ev.stopPropagation();
          highlights[currentProtein][currentDataset] = highlights[currentProtein][currentDataset].filter((r) => r.id !== row.id);
          renderPlot();
          setStatus(`Removed highlight ${{row.start}}-${{row.end}}`);
        }});

        li.appendChild(left);
        li.appendChild(del);
        listEl.appendChild(li);
      }});
    }}

    function buildHighlightShapes() {{
      return highlights[currentProtein][currentDataset].map((row) => ({{
        type: 'rect',
        xref: 'x',
        yref: 'paper',
        x0: row.start - 0.5,
        x1: row.end + 0.5,
        y0: 0,
        y1: 1,
        fillcolor: row.color,
        opacity: 0.18,
        line: {{ width: 1, color: row.color }},
        layer: 'below'
      }}));
    }}

    function addHighlight(start, end, label, color) {{
      const [s, e] = normalizeRange(start, end);
      highlights[currentProtein][currentDataset].push({{
        id: makeHighlightId(),
        start: s,
        end: e,
        label: label || `Region ${{highlights[currentProtein][currentDataset].length + 1}}`,
        color: color || '#f59e0b',
        isSearch: false
      }});
      selectedHighlightId = highlights[currentProtein][currentDataset][highlights[currentProtein][currentDataset].length - 1].id;
      renderPlot();
      setStatus(`Added highlight ${{s}}-${{e}}`);
    }}

    function addHighlightFromInputs() {{
      const datasetRows = proteinDatasets[currentProtein][currentDataset]?.rows || [];
      if (!datasetRows.length) {{
        setStatus('No residues available for this dataset.');
        return;
      }}

      const minResidue = datasetRows[0].pos;
      const maxResidue = datasetRows[datasetRows.length - 1].pos;
      const startVal = Number(startInput.value);
      const endVal = Number(endInput.value);

      if (!Number.isInteger(startVal) || !Number.isInteger(endVal)) {{
        setStatus('Enter integer start and end residue counts.');
        return;
      }}
      if (startVal < minResidue || startVal > maxResidue || endVal < minResidue || endVal > maxResidue) {{
        setStatus(`Residues must be between ${{minResidue}} and ${{maxResidue}} for this dataset.`);
        return;
      }}

      addHighlight(startVal, endVal, labelInput.value.trim(), colorInput.value);
    }}

    function extractTemplateResiduesFromOverview() {{
      const ov = proteinOverviews[currentProtein] || {{}};
      const sections = Array.isArray(ov.sections) ? ov.sections : [];
      const templates = [];

      sections.forEach((text) => {{
        const line = String(text || '');

        const rangeRegex = /\\b([0-9]{{1,5}}) *- *([0-9]{{1,5}})\\b/g;
        let rangeMatch;
        while ((rangeMatch = rangeRegex.exec(line)) !== null) {{
          templates.push({{
            start: Number(rangeMatch[1]),
            end: Number(rangeMatch[2]),
            label: `Template ${{rangeMatch[1]}}-${{rangeMatch[2]}}`,
          }});
        }}

        const residueRegex = /\\b([A-Z])([0-9]{{1,5}})\\b/g;
        let residueMatch;
        while ((residueMatch = residueRegex.exec(line)) !== null) {{
          templates.push({{
            start: Number(residueMatch[2]),
            end: Number(residueMatch[2]),
            label: `Template ${{residueMatch[1]}}${{residueMatch[2]}}`,
          }});
        }}
      }});

      const unique = [];
      const seen = new Set();
      templates.forEach((row) => {{
        const [s, e] = normalizeRange(row.start, row.end);
        const key = `${{s}}-${{e}}`;
        if (seen.has(key)) return;
        seen.add(key);
        unique.push({{ ...row, start: s, end: e }});
      }});
      return unique;
    }}

    function preloadTemplateHighlights() {{
      const datasetRows = proteinDatasets[currentProtein][currentDataset]?.rows || [];
      if (!datasetRows.length) {{
        setStatus('No residues available for this dataset.');
        return;
      }}

      const minResidue = datasetRows[0].pos;
      const maxResidue = datasetRows[datasetRows.length - 1].pos;
      const templates = extractTemplateResiduesFromOverview();
      if (!templates.length) {{
        setStatus('No residue templates found in the overview sections.');
        return;
      }}

      const kept = [];
      let skipped = 0;
      templates.forEach((tpl) => {{
        if (tpl.start < minResidue || tpl.end > maxResidue) {{
          skipped += 1;
          return;
        }}
        kept.push(tpl);
      }});

      if (!kept.length) {{
        setStatus(`No template residues fall within this dataset range (${{minResidue}}-${{maxResidue}}).`);
        return;
      }}

      highlights[currentProtein][currentDataset] = highlights[currentProtein][currentDataset].filter((row) => !row.isTemplate);
      kept.forEach((tpl) => {{
        highlights[currentProtein][currentDataset].push({{
          id: makeHighlightId(),
          start: tpl.start,
          end: tpl.end,
          label: tpl.label,
          color: '#14b8a6',
          isSearch: false,
          isTemplate: true,
        }});
      }});
      selectedHighlightId = highlights[currentProtein][currentDataset][highlights[currentProtein][currentDataset].length - 1].id;
      renderPlot();
      if (skipped > 0) {{
        setStatus(`Loaded ${{kept.length}} template highlight(s); skipped ${{skipped}} outside dataset range.`);
      }} else {{
        setStatus(`Loaded ${{kept.length}} template highlight(s) from overview residues.`);
      }}
    }}

    function addSearchHighlights(query) {{
      const datasetObj = proteinDatasets[currentProtein][currentDataset];
      const rows = datasetObj.rows;
      const motif = query.trim().toUpperCase();
      if (!motif) {{
        setStatus('Enter a sequence motif to search.');
        return;
      }}

      const sequence = rows.map((r) => r.aa).join('');
      const matches = [];
      let start = 0;
      while (true) {{
        const idx = sequence.indexOf(motif, start);
        if (idx === -1) break;
        matches.push([rows[idx].pos, rows[idx + motif.length - 1].pos]);
        start = idx + 1;
      }}

      removeSearchHighlights();
      if (!matches.length) {{
        renderPlot();
        setStatus(`No matches found for ${{motif}}.`);
        return;
      }}

      matches.forEach(([s, e], i) => {{
        highlights[currentProtein][currentDataset].push({{
          id: makeHighlightId(),
          start: s,
          end: e,
          label: `${{motif}} #${{i + 1}}`,
          color: '#ef4444',
          isSearch: true
        }});
      }});
      if (highlights[currentProtein][currentDataset].length) {{
        selectedHighlightId = highlights[currentProtein][currentDataset][highlights[currentProtein][currentDataset].length - matches.length].id;
      }}
      renderPlot();
      setStatus(`Found ${{matches.length}} match(es) for ${{motif}}.`);
    }}

    function renderMsaSlice() {{
      const datasetObj = proteinDatasets[currentProtein][currentDataset];
      const msa = datasetObj.msa;
      const hl = getSelectedHighlight();
      if (!hl) {{
        msaMetaEl.textContent = 'Select or create a highlight to view the aligned region.';
        msaViewEl.innerHTML = '';
        return;
      }}
      if (!msa || !msa.pos_to_col) {{
        msaMetaEl.textContent = `MSA slice unavailable for ${{hl.label}}.`;
        msaViewEl.innerHTML = '';
        return;
      }}

      const startCol = msa.pos_to_col[String(hl.start)] ?? msa.pos_to_col[hl.start];
      const endCol = msa.pos_to_col[String(hl.end)] ?? msa.pos_to_col[hl.end];
      if (startCol === undefined || endCol === undefined) {{
        msaMetaEl.textContent = `Could not map ${{hl.start}}-${{hl.end}} onto the aligned MSA.`;
        msaViewEl.innerHTML = '';
        return;
      }}

      const flank = 5;
      const fromCol = Math.max(0, Math.min(startCol, endCol) - flank);
      const toCol = Math.min(msa.aligned_length - 1, Math.max(startCol, endCol) + flank);
      msaMetaEl.textContent = `${{hl.label}} | query residues ${{hl.start}}-${{hl.end}} | aligned columns ${{fromCol + 1}}-${{toCol + 1}}`;
      const wrapCols = 80;
      const chunks = [];
      for (let chunkStart = fromCol; chunkStart <= toCol; chunkStart += wrapCols) {{
        const chunkEnd = Math.min(toCol, chunkStart + wrapCols - 1);
        const rowsHtml = msa.records.map((r) => `
          <div class="msa-row">
            <div class="msa-name">${{escapeHtml(r.name)}}</div>
            <div class="msa-seq">${{renderColoredSeq(r.seq.slice(chunkStart, chunkEnd + 1))}}</div>
          </div>
        `).join('');
        chunks.push(`
          <div class="msa-chunk">
            <div class="msa-chunk-head">Columns ${{chunkStart + 1}}-${{chunkEnd + 1}}</div>
            ${{rowsHtml}}
          </div>
        `);
      }}
      msaViewEl.innerHTML = chunks.join('');
    }}

    function renderPlot() {{
      const datasetObj = proteinDatasets[currentProtein][currentDataset];
      const data = datasetObj.rows;
      const x = data.map((r) => r.pos);
      const score = data.map((r) => r.score);
      const grade = data.map((r) => r.grade);
      const custom = data.map((r) => [r.aa, r.grade, r.low_conf ? 'yes' : 'no', r.score]);
      const humanPresence = datasetObj.human_presence || null;
      const nSequences = Number.isFinite(datasetObj.n_sequences) ? datasetObj.n_sequences : null;
      const msaSource = datasetObj.msa_source || null;
      const showHumanTrack =
        (currentDataset === 'Full' || currentDataset === 'Invertebrates') &&
        Array.isArray(humanPresence) &&
        humanPresence.length === x.length;
      const humanPresenceY = showHumanTrack ? humanPresence.map((r) => r.present) : null;
      const humanPresenceCustom = showHumanTrack
        ? humanPresence.map((r) => [r.aa, r.count])
        : null;

      const scoreTrace = {{
        x: x,
        y: score,
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Normalized score',
        line: {{ color: '#1d4ed8', width: 1.8 }},
        marker: {{ size: 5, color: '#1d4ed8' }},
        customdata: custom,
        hovertemplate:
          'Residue: %{{x}}<br>' +
          'AA: %{{customdata[0]}}<br>' +
          'Score: %{{y:.4f}}<br>' +
          'Grade: %{{customdata[1]}}<br>' +
          'Low confidence: %{{customdata[2]}}<extra></extra>',
        yaxis: 'y'
      }};

      const gradeTrace = {{
        x: x,
        y: grade,
        type: 'scatter',
        mode: 'markers',
        name: 'ConSurf grade',
        marker: {{
          size: 7,
          color: grade,
          colorscale: 'Viridis',
          cmin: 1,
          cmax: 9,
          opacity: 0.85,
          colorbar: {{ title: 'Grade' }}
        }},
        customdata: custom,
        hovertemplate:
          'Residue: %{{x}}<br>' +
          'AA: %{{customdata[0]}}<br>' +
          'Grade: %{{y}}<br>' +
          'Score: %{{customdata[3]:.4f}}<br>' +
          'Low confidence: %{{customdata[2]}}<extra></extra>',
        yaxis: 'y2'
      }};

      const humanPresenceTrace = {{
        x: x,
        y: humanPresenceY,
        type: 'scatter',
        mode: 'lines',
        name: 'Human sequence present',
        line: {{ color: '#059669', width: 6, shape: 'hv' }},
        customdata: humanPresenceCustom,
        hovertemplate:
          'Residue: %{{x}}<br>' +
          'Human sequence present here: yes<br>' +
          'Human residue: %{{customdata[0]}}<br>' +
          'Human residue number: %{{customdata[1]}}<extra></extra>',
        yaxis: 'y3'
      }};

      const traces = showHumanTrack ? [scoreTrace, gradeTrace, humanPresenceTrace] : [scoreTrace, gradeTrace];
      const titleSuffix = nSequences ? ` (N=${{nSequences}})` : '';

      const layout = {{
        title: `${{currentProtein}} - ${{currentDataset}}${{titleSuffix}}`,
        margin: {{ l: 65, r: 65, t: 55, b: 55 }},
        hovermode: 'x unified',
        xaxis: {{ title: 'Residue position' }},
        yaxis: {{
          title: 'Normalized score',
          zeroline: true,
          zerolinecolor: '#94a3b8',
          domain: showHumanTrack ? [0.46, 1.0] : [0.35, 1.0]
        }},
        yaxis2: {{
          title: 'ConSurf grade',
          domain: showHumanTrack ? [0.16, 0.32] : [0.0, 0.24],
          range: [0.5, 9.5],
          tickvals: [1,2,3,4,5,6,7,8,9]
        }},
        yaxis3: {{
          title: 'Human present',
          domain: showHumanTrack ? [0.0, 0.07] : [0.0, 0.0],
          range: [0, 1.2],
          tickvals: [1],
          ticktext: ['yes'],
          visible: showHumanTrack,
          showgrid: false,
          zeroline: false
        }},
        shapes: buildHighlightShapes(),
        legend: {{ orientation: 'h' }}
      }};

      Plotly.newPlot(plotEl, traces, layout, {{ responsive: true }});
      if (nSequences && msaSource) {{
        datasetMetaEl.textContent = `N = ${{nSequences}} sequences | source: ${{msaSource}}`;
      }} else if (nSequences) {{
        datasetMetaEl.textContent = `N = ${{nSequences}} sequences`;
      }} else {{
        datasetMetaEl.textContent = 'N unavailable for this dataset.';
      }}

      if (typeof plotEl.removeAllListeners === 'function') {{
        plotEl.removeAllListeners('plotly_click');
      }}
      plotEl.on('plotly_click', (ev) => {{
        if (!pickMode) return;
        const xVal = ev?.points?.[0]?.x;
        if (!Number.isFinite(xVal)) return;
        const residue = Math.round(xVal);
        if (pickStart === null) {{
          pickStart = residue;
          setStatus(`Pick mode: start set to ${{residue}}. Click an end residue.`);
        }} else {{
          addHighlight(pickStart, residue, labelInput.value.trim(), colorInput.value);
          pickStart = null;
          pickMode = false;
          pickBtn.textContent = 'Pick from plot';
        }}
      }});

      renderHighlightList();
      if (!pickMode) {{
        if ((currentDataset === 'Full' || currentDataset === 'Invertebrates') && !showHumanTrack) {{
          setStatus('Click Pick from plot, then click two residues to create a highlight. Human-presence track unavailable for this dataset.');
        }} else {{
          setStatus('Click Pick from plot, then click two residues to create a highlight.');
        }}
      }}
    }}

    pickBtn.addEventListener('click', () => {{
      pickMode = !pickMode;
      pickStart = null;
      if (pickMode) {{
        setStatus('Pick mode: click first residue, then second residue.');
        pickBtn.textContent = 'Cancel pick';
      }} else {{
        setStatus('Pick mode cancelled.');
        pickBtn.textContent = 'Pick from plot';
      }}
    }});

    clearBtn.addEventListener('click', () => {{
      highlights[currentProtein][currentDataset] = [];
      selectedHighlightId = null;
      closeMsaModal();
      renderPlot();
      setStatus('Cleared highlights for current selection.');
    }});

    addBtn.addEventListener('click', addHighlightFromInputs);
    templateBtn.addEventListener('click', preloadTemplateHighlights);
    [startInput, endInput].forEach((inputEl) => {{
      inputEl.addEventListener('keydown', (ev) => {{
        if (ev.key === 'Enter') {{
          addHighlightFromInputs();
        }}
      }});
    }});

    msaCloseBtn.addEventListener('click', closeMsaModal);
    msaOverlayEl.addEventListener('click', (ev) => {{
      if (ev.target === msaOverlayEl) closeMsaModal();
    }});

    findBtn.addEventListener('click', () => {{
      addSearchHighlights(seqQueryInput.value);
    }});

    seqQueryInput.addEventListener('keydown', (ev) => {{
      if (ev.key === 'Enter') {{
        addSearchHighlights(seqQueryInput.value);
      }}
    }});

    datasetSelect.addEventListener('change', () => {{
      currentDataset = datasetSelect.value;
      renderPlot();
    }});

    function renderAll() {{
      selectedHighlightId = null;
      renderTabs();
      renderOverview();
      renderDatasetSelect();
      renderPlot();
    }}

    renderAll();
  </script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Build multi-protein ConSurf viewer site")
    parser.add_argument("--output", type=Path, default=Path("index.html"), help="Output HTML path")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Path to alignment_work project root (default: parent of this script dir)",
    )
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="Use Plotly from CDN instead of embedding JS",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_root = args.project_root.resolve() if args.project_root else script_dir.parent.resolve()

    payload = gather_data(project_root)
    plotly_script_tag = make_plotly_script(embed_plotly=not args.cdn)
    html = build_html(payload, plotly_script_tag)

    out_path = args.output if args.output.is_absolute() else script_dir / args.output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    mode = "cdn" if args.cdn else "standalone"
    print(f"Wrote {mode} site: {out_path}")


if __name__ == "__main__":
    main()
