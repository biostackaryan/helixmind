"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

import requests
from dash import html

def fetch_structure_data(query):
    """Fetch structure data from PDB or PubChem."""
    try:
        if query.startswith("pdb:"):
            pdb_id = query.split(":")[1].upper()
            urls = [
                f"https://files.rcsb.org/download/{pdb_id}.pdb",
                f"https://models.rcsb.org/{pdb_id}.pdb"
            ]
            for url in urls:
                resp = requests.get(url)
                if resp.status_code == 200:
                    return resp.text, "pdb"

        elif query.startswith("cid:"):
            cid = query.split(":")[1]
            urls = [
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d",
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
            ]
            for url in urls:
                resp = requests.get(url)
                if resp.status_code == 200:
                    return resp.text, "sdf"

        else:
            raise ValueError("Query must start with 'pdb:' or 'cid:'")

    except Exception as e:
        return f"[Error] {e}", None

    return None, None


def render_structure(query):
    """Return a Dash HTML component with embedded 3Dmol.js viewer."""

    data, fmt = fetch_structure_data(query)
    if not data:
        return html.Div("Could not load structure.")

    # Escape backticks and backslashes in data for JS template literal safety
    escaped_data = data.replace("\\", "\\\\").replace("`", "\\`")

    js_code = f"""
    <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <script>
    (function() {{
        let viewerDiv = document.getElementById("viewer-div");
        if (!viewerDiv) return;
        let viewer = $3Dmol.createViewer(viewerDiv, {{backgroundColor: "black"}});
        viewer.addModel(`{escaped_data}`, "{fmt}");
        if ("{fmt}" === "pdb") {{
            viewer.setStyle({{"cartoon": {{color: "spectrum"}}}});
            viewer.addStyle({{"hetflag": true}}, {{"stick": {{}}}});
        }} else {{
            viewer.setStyle({{"stick": {{}}}});
        }}
        viewer.zoomTo();
        viewer.render();
    }})();
    </script>
    """

    return html.Div([
        html.Div(id="viewer-div", style={"width": "800px", "height": "600px"}),
        html.Div(dangerously_set_inner_html=js_code)
    ])
