"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

from dash import html, dcc, Input, Output, State
from modules.structure_viewer import render_structure
from modules_dash import app  # import app instance

layout = html.Div([
    html.H3("3D Structure Viewer"),
    html.P("Enter a structure query:"),
    dcc.Input(
        id="structure-query",
        type="text",
        placeholder="Example: pdb:1HHP or cid:2244",
        style={"width": "300px", "marginRight": "10px"}
    ),
    html.Button("Load Structure", id="load-structure-btn", n_clicks=0),
    html.Div(id="structure-output", style={"marginTop": "20px"})
])

@app.callback(
    Output("structure-output", "children"),
    Input("load-structure-btn", "n_clicks"),
    State("structure-query", "value"),
    prevent_initial_call=True
)
def update_structure_viewer(n_clicks, query):
    if not query:
        return html.Div("Please enter a valid query (pdb:XXXX or cid:NNNN).")
    return render_structure(query)
