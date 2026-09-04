import argparse
import json
from pathlib import Path
from urllib.request import urlopen


PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"


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
            color = int(color_raw.replace("*", ""))
            low_conf = "*" in color_raw
            rows.append(
                {
                    "pos": pos,
                    "aa": aa,
                    "score": score,
                    "grade": color,
                    "low_conf": low_conf,
                }
            )

    if not rows:
        raise ValueError(f"No rows parsed from {grades_path}")

    rows.sort(key=lambda r: r["pos"])
    return rows


def make_plotly_script(embed_plotly: bool):
  if not embed_plotly:
    return f'<script src="{PLOTLY_CDN}"></script>'

  with urlopen(PLOTLY_CDN, timeout=30) as response:
    plotly_js = response.read().decode("utf-8")
  return f"<script>{plotly_js}</script>"


def build_html(datasets, plotly_script_tag: str):
    payload = json.dumps(datasets)
    return f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>RAD21 Interactive ConSurf Viewer</title>
  {plotly_script_tag}
  <style>
    body {{
      font-family: Segoe UI, Arial, sans-serif;
      margin: 0;
      background: #f7f9fc;
      color: #1f2937;
    }}
    .wrap {{
      max-width: 1200px;
      margin: 20px auto;
      padding: 0 16px;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
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
    .controls {{
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }}
    .status {{
      font-size: 12px;
      color: #334155;
      min-height: 18px;
    }}
    .panel {{
      display: grid;
      grid-template-columns: 1fr 320px;
      gap: 12px;
    }}
    .hint {{
      font-size: 13px;
      color: #475569;
    }}
    #plot {{
      width: 100%;
      height: 720px;
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      background: #fff;
    }}
    .highlights {{
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      background: #fff;
      padding: 10px;
      max-height: 720px;
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
    @media (max-width: 980px) {{
      .panel {{
        grid-template-columns: 1fr;
      }}
      .highlights {{
        max-height: none;
      }}
    }}
  </style>
</head>
<body>
  <div class=\"wrap\">
    <h2>RAD21 Interactive ConSurf Viewer</h2>
    <div class=\"toolbar\">
      <label for=\"dataset\">Dataset:</label>
      <select id=\"dataset\"></select>
      <div class=\"hint\">Hover points to see exact residue number, amino acid, normalized score, and grade.</div>
    </div>
    <div class="toolbar controls">
      <input id="hl-label" type="text" placeholder="Label" value="Region" />
      <input id="hl-start" type="number" placeholder="Start" />
      <input id="hl-end" type="number" placeholder="End" />
      <input id="hl-color" type="color" value="#f59e0b" />
      <button id="hl-add" type="button">Add highlight</button>
      <button id="hl-pick" type="button">Pick from plot</button>
      <button id="hl-clear" type="button">Clear all</button>
    </div>
    <div class="status" id="status"></div>
    <div class="panel">
      <div id="plot"></div>
      <div class="highlights">
        <h3>Highlights</h3>
        <ul id="highlight-list"></ul>
      </div>
    </div>
  </div>

  <script>
    const datasets = {payload};
    const select = document.getElementById('dataset');
    const plotEl = document.getElementById('plot');
    const listEl = document.getElementById('highlight-list');
    const statusEl = document.getElementById('status');
    const addBtn = document.getElementById('hl-add');
    const pickBtn = document.getElementById('hl-pick');
    const clearBtn = document.getElementById('hl-clear');
    const labelInput = document.getElementById('hl-label');
    const startInput = document.getElementById('hl-start');
    const endInput = document.getElementById('hl-end');
    const colorInput = document.getElementById('hl-color');

    const highlightsByDataset = {{}};
    Object.keys(datasets).forEach((name) => {{
      highlightsByDataset[name] = [];
    }});

    let pickMode = false;
    let pickStart = null;

    Object.keys(datasets).forEach((name) => {{
      const opt = document.createElement('option');
      opt.value = name;
      opt.textContent = name;
      select.appendChild(opt);
    }});

    function setStatus(message) {{
      statusEl.textContent = message;
    }}

    function normalizeRange(a, b) {{
      return [Math.min(a, b), Math.max(a, b)];
    }}

    function makeHighlightId() {{
      return `hl_${{Date.now()}}_${{Math.floor(Math.random() * 10000)}}`;
    }}

    function getCurrentDataset() {{
      return select.value;
    }}

    function renderHighlightList(datasetName) {{
      const rows = highlightsByDataset[datasetName];
      listEl.innerHTML = '';
      if (!rows.length) {{
        const empty = document.createElement('div');
        empty.className = 'empty';
        empty.textContent = 'No highlights yet.';
        listEl.appendChild(empty);
        return;
      }}

      rows.forEach((row) => {{
        const li = document.createElement('li');
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
        del.addEventListener('click', () => {{
          highlightsByDataset[datasetName] = highlightsByDataset[datasetName].filter((r) => r.id !== row.id);
          render(datasetName);
          setStatus(`Removed highlight ${{row.start}}-${{row.end}}`);
        }});

        li.appendChild(left);
        li.appendChild(del);
        listEl.appendChild(li);
      }});
    }}

    function buildHighlightShapes(datasetName) {{
      const rows = highlightsByDataset[datasetName];
      return rows.map((row) => ({{
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

    function addHighlight(datasetName, start, end, label, color) {{
      const [s, e] = normalizeRange(start, end);
      highlightsByDataset[datasetName].push({{
        id: makeHighlightId(),
        start: s,
        end: e,
        label: label || `Region ${{highlightsByDataset[datasetName].length + 1}}`,
        color: color || '#f59e0b'
      }});
      render(datasetName);
      setStatus(`Added highlight ${{s}}-${{e}}`);
    }}

    function render(name) {{
      const data = datasets[name];
      const x = data.map(r => r.pos);
      const score = data.map(r => r.score);
      const grade = data.map(r => r.grade);

      const custom = data.map((r) => [r.aa, r.grade, r.low_conf ? 'yes' : 'no', r.score]);

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

      const layout = {{
        title: name,
        margin: {{ l: 65, r: 65, t: 55, b: 55 }},
        hovermode: 'x unified',
        xaxis: {{ title: 'Residue position' }},
        yaxis: {{
          title: 'Normalized score',
          zeroline: true,
          zerolinecolor: '#94a3b8',
          domain: [0.35, 1.0]
        }},
        yaxis2: {{
          title: 'ConSurf grade',
          domain: [0.0, 0.24],
          range: [0.5, 9.5],
          tickvals: [1,2,3,4,5,6,7,8,9]
        }},
        shapes: buildHighlightShapes(name),
        legend: {{ orientation: 'h' }}
      }};

      Plotly.newPlot(plotEl, [scoreTrace, gradeTrace], layout, {{responsive: true}});

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
          const dataset = getCurrentDataset();
          addHighlight(dataset, pickStart, residue, labelInput.value.trim(), colorInput.value);
          pickStart = null;
          pickMode = false;
          pickBtn.textContent = 'Pick from plot';
        }}
      }});

      renderHighlightList(name);

      if (!pickMode) {{
        setStatus('Add a highlight using Start/End or click Pick from plot.');
      }}
    }}

    addBtn.addEventListener('click', () => {{
      const dataset = getCurrentDataset();
      const start = Number(startInput.value);
      const end = Number(endInput.value);
      if (!Number.isFinite(start) || !Number.isFinite(end)) {{
        setStatus('Enter numeric Start and End values.');
        return;
      }}
      addHighlight(dataset, Math.round(start), Math.round(end), labelInput.value.trim(), colorInput.value);
    }});

    pickBtn.addEventListener('click', () => {{
      pickMode = !pickMode;
      pickStart = null;
      if (pickMode) {{
        setStatus('Pick mode: click first residue, then click second residue to create a highlight.');
        pickBtn.textContent = 'Cancel pick';
      }} else {{
        setStatus('Pick mode cancelled.');
        pickBtn.textContent = 'Pick from plot';
      }}
    }});

    clearBtn.addEventListener('click', () => {{
      const dataset = getCurrentDataset();
      highlightsByDataset[dataset] = [];
      render(dataset);
      setStatus('Cleared highlights for current dataset.');
    }});

    select.addEventListener('change', () => render(select.value));
    render(select.value);
  </script>
</body>
</html>
"""


def parse_spec(text):
    if "=" not in text:
        raise ValueError("Each input must be in Label=path format")
    label, path = text.split("=", 1)
    return label.strip(), Path(path.strip())


def main():
    parser = argparse.ArgumentParser(description="Build interactive HTML viewer for ConSurf grades")
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Input spec in Label=path_to_grades format. Repeat for multiple datasets.",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output HTML path")
    parser.add_argument(
      "--standalone",
      action="store_true",
      help="Embed Plotly JS in the HTML for offline/shareable viewing.",
    )
    args = parser.parse_args()

    datasets = {}
    for spec in args.input:
        label, path = parse_spec(spec)
        datasets[label] = parse_grades(path)

    plotly_script_tag = make_plotly_script(embed_plotly=args.standalone)
    html = build_html(datasets, plotly_script_tag)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    mode = "standalone" if args.standalone else "cdn"
    print(f"Wrote interactive viewer ({mode}): {args.output}")


if __name__ == "__main__":
    main()
