"""
helixmind_dash.py
Helix Mind - Dash app main file (refreshed)

"""
import sys
import os
from pathlib import Path
import io
import base64
import tempfile
import subprocess
from datetime import datetime
import re

# make sure project root is importable so `helixmind` package resolves
ROOT = Path(os.path.dirname(os.path.abspath(__file__)))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dotenv import load_dotenv
load_dotenv()

# Dash + Flask
from dash import Dash, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_table
from flask import send_from_directory, abort

# Bio / data libs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
Entrez.email = os.getenv("NCBI_EMAIL", "example@example.com")  # good practice

import pandas as pd
import plotly.express as px
import requests

# constants & folders
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
DATA_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

LOCAL_DB_DEFAULT = "localdb/queries_db"
MAX_SEQ_DISPLAY = 10  # avoid huge interactive runs
ONLINE_BLAST_MAX = 1  # limit online preview to 1 query to be nice to NCBI

# attempt to import internal helper modules, if present
try:
    from modules.fasta_parser import parse_fasta_file as helix_parse_fasta
except Exception:
    helix_parse_fasta = None

try:
    # kegg_module: prefer search_kegg, fetch_kegg_details, fetch_kegg_enzyme_info
    from modules.kegg_module import search_kegg as helix_search_kegg
    from modules.kegg_module import fetch_kegg_details as helix_fetch_kegg_details
    from modules.kegg_module import fetch_kegg_enzyme_info as helix_fetch_kegg_enzyme_info
except Exception:
    helix_search_kegg = None
    helix_fetch_kegg_details = None
    helix_fetch_kegg_enzyme_info = None

try:
    from modules.pubmed_module import fetch_pubmed_articles as helix_fetch_pubmed
except Exception:
    helix_fetch_pubmed = None

try:
    from modules.blast_module import run_blast as helix_run_blast
except Exception:
    helix_run_blast = None

try:
    from modules.structure_viewer import render_structure as helix_render_structure  # not used in this JS-based viewer
except Exception:
    helix_render_structure = None

# Helper: heatmap generator (returns a Plotly figure)
def generate_heatmap(df: pd.DataFrame, x_label="Features", y_label="Sequences", title="Heatmap"):
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
    fig = px.imshow(
        df.values,
        labels=dict(x=x_label, y=y_label, color="Value"),
        x=df.columns.tolist(),
        y=df.index.tolist(),
        color_continuous_scale=px.colors.sequential.Aggrnyl,
        title=title,
        aspect="auto",
    )
    fig.update_layout(margin=dict(l=40, r=40, t=40, b=40))
    return fig

# Initialize Dash app
external_scripts = ["https://3Dmol.org/build/3Dmol-min.js"]
app = Dash(__name__, suppress_callback_exceptions=True,
           external_stylesheets=[dbc.themes.BOOTSTRAP],
           external_scripts=external_scripts)
server = app.server
app.title = "Helix Mind Bioinformatics Toolkit"

# Serve downloads from results directory (safe)
@server.route("/download/<path:filename>")
def download_file(filename):
    full = RESULTS_DIR / filename
    try:
        resolved = full.resolve()
    except Exception:
        abort(404)
    if not str(resolved).startswith(str(RESULTS_DIR.resolve())):
        abort(404)
    if not full.exists():
        abort(404)
    return send_from_directory(RESULTS_DIR, filename, as_attachment=True)

# --- Small helpers ---
def gc_content(seq_str: str) -> float:
    s = seq_str.upper()
    if not s:
        return 0.0
    return round((s.count("G") + s.count("C")) / len(s) * 100, 2)

def is_ec_number(text: str) -> bool:
    # Simple EC matcher like EC:1.1.1.1 or 1.1.1.1
    return bool(re.match(r"^(EC:)?\d+(\.\d+){1,3}$", text.strip(), flags=re.IGNORECASE))

# App layout - tabs for all features
app.layout = dbc.Container([
    html.H1("Helix Mind", className="my-3 text-center"),

    dcc.Tabs(id="tabs", value="tab-fasta", children=[
        dcc.Tab(label="FASTA Parser", value="tab-fasta"),
        dcc.Tab(label="Heatmap", value="tab-heatmap"),
        dcc.Tab(label="BLAST", value="tab-blast"),
        dcc.Tab(label="KEGG", value="tab-kegg"),
        dcc.Tab(label="PubMed", value="tab-pubmed"),
        dcc.Tab(label="Structure Viewer", value="tab-structure"),
        dcc.Tab(label="ChatGPT Assistant", value="tab-chatgpt"),
    ]),
    html.Div(id="tab-content", className="mt-4"),

    html.Hr(),  # separator line
    html.Footer([
        "© 2025 ",
        html.A("@biostackaryan", href="https://github.com/biostackaryan", target="_blank"),
        " - HelixMind"
    ], style={"textAlign": "center", "marginTop": "20px", "padding": "10px"})
], fluid=True)

# Render tab contents
@app.callback(Output("tab-content", "children"), Input("tabs", "value"))
def render_tab(tab):
    if tab == "tab-fasta":
        return dbc.Row([
            dbc.Col([
                html.H4("FASTA Parser"),
                dcc.Upload(
                    id="upload-fasta",
                    children=html.Div(["Drag and drop or ", html.A("select a FASTA file")]),
                    style={"width": "100%", "height": "70px", "lineHeight": "70px",
                           "borderWidth": "1px", "borderStyle": "dashed",
                           "borderRadius": "5px", "textAlign": "center"},
                    multiple=False
                ),
                html.Br(),
                dbc.Button("Parse FASTA", id="btn-parse-fasta", color="primary"),
                html.Div(id="fasta-parse-output", className="mt-3")
            ], width=10)
        ])

    if tab == "tab-heatmap":
        return dbc.Row([
            dbc.Col([
                html.H4("Generate Heatmap"),
                dcc.Upload(
                    id="upload-heatmap-fasta",
                    children=html.Div(["Upload FASTA for heatmap (GC% & length)"]),
                    style={"width": "100%", "height": "70px", "lineHeight": "70px",
                           "borderWidth": "1px", "borderStyle": "dashed",
                           "borderRadius": "5px", "textAlign": "center"},
                    multiple=False
                ),
                html.Br(),
                dbc.Button("Generate Heatmap", id="btn-gen-heatmap", color="success"),
                html.Div(id="heatmap-status", className="mt-2"),
                dcc.Graph(id="heatmap-figure", style={"display": "none"})
            ], width=10)
        ])

    if tab == "tab-blast":
        return dbc.Row([
            dbc.Col([
                html.H4("BLAST Search (local or NCBI online)"),
                dcc.Upload(
                    id="upload-blast",
                    children=html.Div(["Upload FASTA to BLAST"]),
                    style={"width": "100%", "height": "70px", "lineHeight": "70px",
                           "borderWidth": "1px", "borderStyle": "dashed",
                           "borderRadius": "5px", "textAlign": "center"},
                    multiple=False
                ),
                html.Br(),
                dbc.Row([
                    dbc.Col(dbc.Label("Mode"), width=2),
                    dbc.Col(dcc.Dropdown(id="blast-mode", options=[
                        {"label": "Local BLAST (requires local BLAST DB & blastn)", "value": "local"},
                        {"label": "Online BLAST (NCBI qblast; preview only)", "value": "online"},
                    ], value="local"), width=6),
                ], align="center"),
                html.Br(),
                dbc.Row([
                    dbc.Col(dbc.Label("Database / DB path"), width=2),
                    dbc.Col(dbc.Input(id="blast-db", value=LOCAL_DB_DEFAULT, placeholder="localdb/queries_db or nt"), width=6),
                ], align="center"),
                html.Br(),
                dbc.Row([
                    dbc.Col(dbc.Label("E-value"), width=2),
                    dbc.Col(dbc.Input(id="blast-evalue", value="0.001"), width=6),
                ], align="center"),
                html.Br(),
                dbc.Checklist(options=[{"label": "Show raw BLAST output", "value": "raw"}], value=[], id="blast-show-raw"),
                html.Br(),
                dbc.Button("Run BLAST (limited preview)", id="btn-run-blast", color="danger"),
                html.Div(id="blast-status", className="mt-2"),
            ], width=5),
            dbc.Col([
                html.H5("Parsed BLAST Results"),
                html.Div(id="blast-table-area"),
                html.Hr(),
                html.H5("Raw BLAST Output"),
                html.Div(id="blast-raw-area", style={"whiteSpace": "pre-wrap", "maxHeight": "420px", "overflow": "auto", "background": "#f8f9fa", "padding": "0.5rem"})
            ], width=7)
        ])

    if tab == "tab-kegg":
        return dbc.Row([
            dbc.Col([
                html.H4("KEGG Pathway / EC Search"),
                dcc.Input(id="kegg-query", placeholder="Pathway name/ID (path:hsa00010) or EC (EC:1.1.1.1)", style={"width": "100%"}),
                html.Br(), html.Br(),
                dbc.Button("Search KEGG", id="btn-search-kegg", color="primary"),
                html.Div(id="kegg-results-area", className="mt-3")
            ], width=8)
        ])

    if tab == "tab-pubmed":
        return dbc.Row([
            dbc.Col([
                html.H4("PubMed Search"),
                dcc.Input(id="pubmed-query", placeholder="Search term (gene, keyword, author...)", style={"width": "100%"}),
                html.Br(), html.Br(),
                dbc.Button("Search PubMed", id="btn-search-pubmed", color="primary"),
                html.Div(id="pubmed-results-area", className="mt-3")
            ], width=8)
        ])

    if tab == "tab-structure":
        return dbc.Row([
            dbc.Col([
                html.H4("3D Structure Viewer (PDB or PubChem CID)"),
                dcc.Input(id="structure-query", placeholder="pdb:1CRN or cid:2244", style={"width": "100%"}),
                html.Br(), html.Br(),
                dbc.Button("Load Structure", id="btn-load-structure", color="info"),
                html.Div(id="structure-output", className="mt-3", style={"whiteSpace": "pre-wrap"}),
                html.Div(id="viewer-div", style={"width": "100%", "height": "500px", "background": "#000"})
            ], width=10)
        ])

    if tab == "tab-chatgpt":
        return dbc.Row([
            dbc.Col([
                html.H4("ChatGPT Bioinformatics Assistant (stub if no module)"),
                dcc.Textarea(id="chat-prompt", placeholder="Ask something bioinformatics-related...", style={"width": "100%", "height": "150px"}),
                html.Br(),
                dbc.Button("Ask", id="btn-ask-chat", color="primary"),
                html.Div(id="chat-response", className="mt-3", style={"whiteSpace": "pre-wrap"})
            ], width=8)
        ])

    return html.Div("Select a tab")

# ---------- FASTA parser callback ----------
@app.callback(
    Output("fasta-parse-output", "children"),
    Input("btn-parse-fasta", "n_clicks"),
    State("upload-fasta", "contents"),
    State("upload-fasta", "filename"),
    prevent_initial_call=True
)
def parse_fasta(n_clicks, contents, filename):
    if not contents:
        return dbc.Alert("Please upload a FASTA file.", color="warning")

    # decode
    try:
        content_type, content_string = contents.split(",", 1)
        decoded = base64.b64decode(content_string)
        fasta_text = decoded.decode("utf-8")
    except Exception as e:
        return dbc.Alert(f"Failed to decode uploaded file: {e}", color="danger")

    records = []
    # Prefer user module parser if available, but ensure we still compute stats
    if helix_parse_fasta:
        try:
            parsed = helix_parse_fasta(io.StringIO(fasta_text))
            # normalize to list of SeqRecord or dict of id->sequence
            if isinstance(parsed, list):
                # If they're SeqRecord-like
                for x in parsed:
                    if isinstance(x, SeqRecord):
                        records.append(x)
                    else:
                        try:
                            # expect tuple (id, seq) or obj with .id and .seq
                            rec_id = getattr(x, "id", None) or str(x[0])
                            rec_seq = str(getattr(x, "seq", None) or x[1])
                            records.append(SeqRecord(Seq(rec_seq), id=rec_id, description=""))
                        except Exception:
                            pass
            elif isinstance(parsed, dict):
                for k, v in parsed.items():
                    records.append(SeqRecord(Seq(str(v)), id=str(k), description=""))
            elif isinstance(parsed, tuple) and len(parsed) >= 1:
                # assume first element is iterable of records
                for r in parsed[0]:
                    if isinstance(r, SeqRecord):
                        records.append(r)
                    else:
                        try:
                            rec_id = getattr(r, "id", None) or str(r[0])
                            rec_seq = str(getattr(r, "seq", None) or r[1])
                            records.append(SeqRecord(Seq(rec_seq), id=rec_id, description=""))
                        except Exception:
                            pass
        except Exception:
            records = []

    # Fallback: simple Biopython parse
    if not records:
        try:
            handle = io.StringIO(fasta_text)
            records = list(SeqIO.parse(handle, "fasta"))
        except Exception as e:
            return dbc.Alert(f"Error parsing FASTA: {e}", color="danger")

    if not records:
        return dbc.Alert("No sequences found in the uploaded FASTA.", color="warning")

    # compute stats
    infos = []
    for r in records:
        seq_str = str(r.seq)
        infos.append({
            "ID": r.id,
            "Length_bp": len(seq_str),
            "GC_percent": gc_content(seq_str)
        })
    df = pd.DataFrame(infos).sort_values("Length_bp", ascending=False)

    longest_row = df.iloc[0]
    shortest_row = df.iloc[-1]

    # save uploaded file
    save_path = DATA_DIR / filename
    try:
        save_path.write_bytes(decoded)
    except Exception as e:
        return dbc.Alert(f"Parsed, but failed to save file: {e}", color="warning")

    table = dash_table.DataTable(
        columns=[{"name": c, "id": c} for c in df.columns],
        data=df.to_dict("records"),
        page_size=10,
        style_table={"overflowX": "auto"},
    )

    return html.Div([
        html.P(f"File saved to: {save_path}"),
        html.P(f"Total sequences: {len(df)}"),
        html.P(f"Longest: {longest_row['ID']} ({int(longest_row['Length_bp'])} bp, GC {longest_row['GC_percent']}%)"),
        html.P(f"Shortest: {shortest_row['ID']} ({int(shortest_row['Length_bp'])} bp, GC {shortest_row['GC_percent']}%)"),
        html.Hr(),
        html.H6("Per-sequence summary"),
        table
    ])

# ---------- Heatmap callback (unchanged) ----------
@app.callback(
    Output("heatmap-status", "children"),
    Output("heatmap-figure", "figure"),
    Output("heatmap-figure", "style"),
    Input("btn-gen-heatmap", "n_clicks"),
    State("upload-heatmap-fasta", "contents"),
    State("upload-heatmap-fasta", "filename"),
    prevent_initial_call=True
)
def make_heatmap(n_clicks, contents, filename):
    if not contents:
        return dbc.Alert("Please upload a FASTA file.", color="warning"), {}, {"display": "none"}

    try:
        _, payload = contents.split(",", 1)
        decoded = base64.b64decode(payload)
        fasta_text = decoded.decode("utf-8")
        handle = io.StringIO(fasta_text)
        records = list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        return dbc.Alert(f"Failed to read FASTA: {e}", color="danger"), {}, {"display": "none"}

    if not records:
        return dbc.Alert("No sequences found in FASTA.", color="warning"), {}, {"display": "none"}

    # create df: GC% and length per sequence id
    rows = {}
    for r in records:
        s = str(r.seq)
        if len(s) == 0:
            continue
        gc = (s.upper().count("G") + s.upper().count("C")) / len(s) * 100
        rows[r.id] = {"GC_pct": round(gc, 2), "Length": len(s)}
    df = pd.DataFrame.from_dict(rows, orient="index")

    fig = generate_heatmap(df, x_label="Features", y_label="Sequences", title=f"GC% and Length ({filename})")
    style = {"display": "block"}
    status = html.P(f"Heatmap generated for {len(df)} sequences from {filename}.")
    return status, fig, style

# ---------- BLAST callback (local + online preview) ----------
@app.callback(
    Output("blast-status", "children"),
    Output("blast-table-area", "children"),
    Output("blast-raw-area", "children"),
    Input("btn-run-blast", "n_clicks"),
    State("upload-blast", "contents"),
    State("upload-blast", "filename"),
    State("blast-mode", "value"),
    State("blast-db", "value"),
    State("blast-evalue", "value"),
    State("blast-show-raw", "value"),
    prevent_initial_call=True
)
def run_blast_cb(n_clicks, contents, filename, mode, db, evalue, show_raw):
    import tempfile, base64, io, os, subprocess
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIXML
    import pandas as pd
    from datetime import datetime
    import dash_table
    from dash import html
    import dash_bootstrap_components as dbc

    MAX_SEQ_DISPLAY = 5
    RESULTS_DIR = Path("results")
    RESULTS_DIR.mkdir(exist_ok=True)

    # Step 1 — Check upload
    if not contents:
        return dbc.Alert("Please upload a FASTA file to BLAST.", color="warning"), None, None

    try:
        _, payload = contents.split(",", 1)
        decoded = base64.b64decode(payload)
        fasta_text = decoded.decode("utf-8")
    except Exception as e:
        return dbc.Alert(f"Invalid upload: {e}", color="danger"), None, None

    records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
    if not records:
        return dbc.Alert("No sequences found in the uploaded FASTA.", color="warning"), None, None

    # Step 2 — Limit preview for online mode
    preview_records = records[:1] if mode == "online" else records[:MAX_SEQ_DISPLAY]

    status_msgs = []
    raw_outputs = []
    frames = []

    # Step 3 — Local BLAST
    if mode == "local":
        # Check if BLAST+ is installed
        try:
            subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        except Exception:
            return dbc.Alert("Local BLAST selected but `blastn` not found. Install BLAST+ or switch to Online mode.", color="danger"), None, None

        for rec in preview_records:
            # Write sequence to temp file
            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmpf:
                tmpf.write(f">{rec.id}\n{str(rec.seq)}\n")
                tmpq = tmpf.name

            try:
                cmd = [
                    "blastn",
                    "-query", tmpq,
                    "-db", db,
                    "-evalue", str(evalue),
                    "-outfmt", "6 qseqid sseqid stitle pident length evalue bitscore",
                    "-max_target_seqs", "10"
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode != 0:
                    raw_outputs.append(f"[blastn error for {rec.id}]: {proc.stderr}\n")
                else:
                    raw_outputs.append(f"--- {rec.id} (local tabular) ---\n{proc.stdout}\n")
                    if proc.stdout.strip():
                        rows = []
                        for line in proc.stdout.strip().splitlines():
                            parts = line.split("\t")
                            if len(parts) >= 7:
                                rows.append({
                                    "query_id": parts[0],
                                    "subject_id": parts[1],
                                    "subject_title": parts[2],
                                    "identity": float(parts[3]),
                                    "align_len": int(parts[4]),
                                    "evalue": float(parts[5]),
                                    "bitscore": float(parts[6])
                                })
                        if rows:
                            frames.append(pd.DataFrame(rows).assign(query_header=rec.id))
            except Exception as e:
                raw_outputs.append(f"[ERROR for {rec.id}] {e}\n")
            finally:
                try:
                    os.remove(tmpq)
                except Exception:
                    pass

    # Step 4 — Online BLAST
    elif mode == "online":
        rec = preview_records[0]
        try:
            online_db = db if db and db.lower() in {"nt", "refseq_rna", "refseq_genomic"} else "nt"
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database=online_db,
                sequence=str(rec.seq),
                expect=float(evalue),
                hitlist_size=10
            )
            blast_record = NCBIXML.read(result_handle)
            result_handle.close()

            rows = []
            for alignment in blast_record.alignments:
                if not alignment.hsps:
                    continue
                hsp = alignment.hsps[0]
                rows.append({
                    "query_id": rec.id,
                    "subject_id": alignment.hit_id,
                    "subject_title": alignment.hit_def,
                    "identity": round((hsp.identities / hsp.align_length) * 100, 2) if hsp.align_length else 0.0,
                    "align_len": int(hsp.align_length),
                    "evalue": float(hsp.expect),
                    "bitscore": float(hsp.bits)
                })
            if rows:
                frames.append(pd.DataFrame(rows))
            raw_outputs.append(f"--- {rec.id} (online summary) ---\nParsed {len(rows)} hits from NCBI qblast XML.\n")
        except Exception as e:
            raw_outputs.append(f"[ONLINE BLAST ERROR for {rec.id}] {e}\n")

    # Step 5 — Build results table
    if frames:
        df_all = pd.concat(frames, ignore_index=True)
        table = dash_table.DataTable(
            columns=[{"name": c, "id": c} for c in df_all.columns],
            data=df_all.to_dict("records"),
            page_size=10,
            style_table={"overflowX": "auto"}
        )
    else:
        table = html.Div("No parsed results available. Check raw output.")

    # Step 6 — Raw output section
    raw_text = "\n".join(raw_outputs) if raw_outputs else "No raw output."
    raw_area = html.Pre(raw_text) if ("raw" in (show_raw or [])) else html.Pre("(Hidden) Enable 'Show raw BLAST output' to view.")

    # Step 7 — Save raw output
    try:
        out_name = f"blast_preview_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        (RESULTS_DIR / out_name).write_text(raw_text)
        download_link = html.A("Download raw preview results", href=f"/download/{out_name}", target="_blank")
    except Exception:
        download_link = None

    # Step 8 — Status message
    mode_note = "Online mode limited to 1 sequence for preview." if mode == "online" else f"Local mode preview up to {MAX_SEQ_DISPLAY} seqs."
    status_block = html.Div([
        html.P(f"Ran BLAST preview on {len(preview_records)} sequence(s). ({mode_note})"),
        download_link if download_link else html.Span("")
    ])

    return status_block, table, raw_area


# ---------- KEGG callback (enhanced) ----------
@app.callback(
    Output("kegg-results-area", "children"),
    Input("btn-search-kegg", "n_clicks"),
    State("kegg-query", "value"),
    prevent_initial_call=True
)
def kegg_search_cb(n_clicks, query):
    if not query:
        return dbc.Alert(
            "Please type a pathway name/ID (e.g., glycolysis or path:hsa00010) or an EC number (e.g., EC:1.1.1.1).",
            color="warning"
        )

    q = query.strip()

    try:
        # Detect EC number
        if is_ec_number(q):
            norm = q.upper().replace("EC:", "")
            r = requests.get(f"https://rest.kegg.jp/get/ec:{norm}", timeout=15)
            if r.status_code != 200 or not r.text.strip():
                return html.Div(f"No KEGG enzyme entry for EC:{norm}.")

            lines = r.text.strip().splitlines()
            name, definition, sys_name = "", "", ""
            pathways, genes, compounds = [], [], []

            for line in lines:
                if line.startswith("NAME"):
                    name = line.split("NAME", 1)[1].strip()
                elif line.startswith("DEFINITION"):
                    definition = line.split("DEFINITION", 1)[1].strip()
                elif line.startswith("SYSNAME"):
                    sys_name = line.split("SYSNAME", 1)[1].strip()
                elif line.startswith("PATHWAY"):
                    pathways.append(line.split("PATHWAY", 1)[1].strip())
                elif line.startswith("GENES"):
                    genes.append(line.split("GENES", 1)[1].strip())
                elif line.startswith("COMPOUND"):
                    compounds.append(line.split("COMPOUND", 1)[1].strip())

            # Build display
            content = [
                html.H5(f"EC:{norm} — {name}"),
                html.P(f"Systematic name: {sys_name}" if sys_name else ""),
                html.P(f"Definition: {definition}"),
                html.H6("Linked pathways:"),
                html.Ul([html.Li(p) for p in pathways]) if pathways else html.P("None"),
                html.H6("Associated genes:"),
                html.Ul([html.Li(g) for g in genes]) if genes else html.P("None"),
                html.H6("Associated compounds:"),
                html.Ul([html.Li(c) for c in compounds]) if compounds else html.P("None"),
                html.A("View full KEGG entry", href=f"https://www.kegg.jp/entry/ec:{norm}", target="_blank")
            ]
            return html.Div(content)

        else:
            # Pathway search
            # Try to detect if input is already a pathway ID
            pid = q if q.lower().startswith("path:") else None
            if not pid:
                search_res = requests.get(f"https://rest.kegg.jp/find/pathway/{q}", timeout=15)
                if search_res.status_code != 200 or not search_res.text.strip():
                    return html.Div("No KEGG pathway results.")
                # Pick first match for detail
                pid = search_res.text.strip().split("\t")[0]

            # Get detailed pathway info
            r = requests.get(f"https://rest.kegg.jp/get/{pid}", timeout=15)
            if r.status_code != 200 or not r.text.strip():
                return html.Div("No KEGG pathway data found.")

            lines = r.text.strip().splitlines()
            title, description, class_info = "", "", ""
            genes, compounds, enzymes = [], [], []

            section = None
            for line in lines:
                if line.startswith("NAME"):
                    title = line.split("NAME", 1)[1].strip()
                elif line.startswith("DESCRIPTION"):
                    description = line.split("DESCRIPTION", 1)[1].strip()
                elif line.startswith("CLASS"):
                    class_info = line.split("CLASS", 1)[1].strip()
                elif line.startswith("GENE"):
                    section = "GENE"
                    genes.append(line.split("GENE", 1)[1].strip())
                elif line.startswith("COMPOUND"):
                    section = "COMPOUND"
                    compounds.append(line.split("COMPOUND", 1)[1].strip())
                elif line.startswith("ENZYME"):
                    section = "ENZYME"
                    enzymes.append(line.split("ENZYME", 1)[1].strip())
                elif section == "GENE" and line.startswith(" "):
                    genes.append(line.strip())
                elif section == "COMPOUND" and line.startswith(" "):
                    compounds.append(line.strip())
                elif section == "ENZYME" and line.startswith(" "):
                    enzymes.append(line.strip())
                else:
                    section = None

            content = [
                html.H5(f"{pid} — {title}"),
                html.P(f"Description: {description}"),
                html.P(f"Class: {class_info}"),
                html.H6("Genes:"),
                html.Ul([html.Li(g) for g in genes[:50]]) if genes else html.P("None"),
                html.H6("Compounds:"),
                html.Ul([html.Li(c) for c in compounds[:50]]) if compounds else html.P("None"),
                html.H6("Enzymes:"),
                html.Ul([html.Li(e) for e in enzymes[:50]]) if enzymes else html.P("None"),
                html.A("View KEGG pathway", href=f"https://www.kegg.jp/pathway/{pid.replace('path:', '')}", target="_blank"),
                html.Br(),
                html.A("Download KGML", href=f"https://rest.kegg.jp/get/{pid}/kgml", target="_blank"),
                html.Br(),
                html.A("Download PNG", href=f"https://www.kegg.jp/kegg-bin/show_pathway?{pid.replace('path:', '')}", target="_blank"),
            ]
            return html.Div(content)

    except Exception as e:
        return dbc.Alert(f"KEGG REST error: {e}", color="danger")


# ---------- PubMed callback (unchanged behavior) ----------
@app.callback(
    Output("pubmed-results-area", "children"),
    Input("btn-search-pubmed", "n_clicks"),
    State("pubmed-query", "value"),
    prevent_initial_call=True
)
def pubmed_search_cb(n_clicks, query):
    if not query:
        return dbc.Alert("Please enter a search term for PubMed.", color="warning")

    if helix_fetch_pubmed:
        try:
            r = helix_fetch_pubmed(query, max_results=10)
            if r.get("status") == "ok":
                items = []
                for art in r.get("results", []):
                    items.append(html.Div([
                        html.H6(art.get("title", "No Title")),
                        html.Small(f"{art.get('source','')} | {art.get('pubdate','')}"),
                        html.Br(),
                        html.A("View on PubMed", href=f"https://pubmed.ncbi.nlm.nih.gov/{art.get('id','')}/", target="_blank"),
                        html.Hr()
                    ]))
                return html.Div(items)
            else:
                return dbc.Alert(f"PubMed module error: {r.get('error_message','unknown')}", color="danger")
        except Exception as e:
            return dbc.Alert(f"Error calling helixmind.pubmed_module: {e}", color="danger")

    # fallback: NCBI esearch (basic)
    try:
        params = {"db": "pubmed", "term": query, "retmax": 10, "retmode": "json"}
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params=params, timeout=15)
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return html.Div("No PubMed articles found.")
        items = []
        for pmid in ids:
            items.append(html.Div([
                html.H6(f"PubMed Article {pmid}"),
                html.Small("PubMed | Unknown date"),
                html.Br(),
                html.A("View on PubMed", href=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", target="_blank"),
                html.Hr()
            ]))
        return html.Div(items)
    except Exception as e:
        return dbc.Alert(f"PubMed search failed: {e}", color="danger")

# ---------- Structure Viewer (clientside rendering using 3Dmol) ----------
# Stronger validation & clearer messages.
app.clientside_callback(
    """
    function(n_clicks, q) {
        if (!n_clicks) return window.dash_clientside.no_update;
        var el = document.getElementById('viewer-div');
        var out = document.getElementById('structure-output');
        if (!el || !out) return window.dash_clientside.no_update;
        el.innerHTML = '';
        out.innerText = '';
        if (!q) { out.innerText = 'Please enter a query like "pdb:1CRN" or "cid:2244".'; return null; }
        try {
            var viewer = $3Dmol.createViewer(el, {backgroundColor: 'white'});
            var query = q.trim();
            var lower = query.toLowerCase();
            if (lower.startsWith('pdb:')) {
                var id = query.split(':')[1];
                if (!id) { out.innerText = 'Format error: use "pdb:<id>", e.g., pdb:1CRN'; return null; }
                $3Dmol.download('pdb:' + id, viewer, {}, function() {
                    viewer.setStyle({}, {cartoon:{color:'spectrum'}});
                    viewer.zoomTo();
                    viewer.render();
                });
                out.innerText = 'Loaded PDB: ' + id;
            } else if (lower.startsWith('cid:')) {
                var cid = query.split(':')[1];
                if (!cid) { out.innerText = 'Format error: use "cid:<pubchem_cid>", e.g., cid:2244'; return null; }
                var url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + cid + '/SDF?record_type=3d';
                fetch(url).then(resp => {
                    if (!resp.ok) throw new Error('HTTP ' + resp.status);
                    return resp.text();
                }).then(txt => {
                    viewer.addModel(txt, 'sdf');
                    viewer.setStyle({}, {stick:{}});
                    viewer.zoomTo();
                    viewer.render();
                    out.innerText = 'Loaded PubChem CID: ' + cid;
                }).catch(e => {
                    out.innerText = 'Failed to fetch structure for CID ' + cid + ': ' + e;
                });
            } else {
                out.innerText = 'Enter query like "pdb:1CRN" or "cid:2244" (must start with "pdb:" or "cid:").';
            }
        } catch(err) {
            out.innerText = '3Dmol error: ' + err;
        }
        return null;
    }
    """,
    Output("structure-output", "children"),
    Input("btn-load-structure", "n_clicks"),
    State("structure-query", "value"),
    prevent_initial_call=True
)

# ---------- ChatGPT assistant (unchanged) ----------
@app.callback(
    Output("chat-response", "children"),
    Input("btn-ask-chat", "n_clicks"),
    State("chat-prompt", "value"),
    prevent_initial_call=True
)
def chat_assistant(n_clicks, prompt):
    if not prompt:
        return dbc.Alert("Please enter a question.", color="warning")
    # keep user module if present; else stub
    try:
        from modules.chatgpt_module import ask_chatgpt as helix_ask_chatgpt
    except Exception:
        helix_ask_chatgpt = None

    if helix_ask_chatgpt:
        try:
            resp = helix_ask_chatgpt(prompt)
            return html.Pre(resp)
        except Exception as e:
            return dbc.Alert(f"ChatGPT module error: {e}", color="danger")
    else:
        return html.Pre(f"(Stub) AI Response to: {prompt}")

# Run server
if __name__ == "__main__":
    print(f"Starting Helix Mind app from {ROOT}")
    app.run(debug=True, port=8050)
