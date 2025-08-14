# helixmind/heatmap.py
"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""
from typing import Optional, Sequence, Union
import numpy as np
import pandas as pd
import plotly.graph_objects as go

ArrayLike2D = Union[np.ndarray, pd.DataFrame, Sequence[Sequence[float]]]

def _to_matrix_and_labels(
    data: ArrayLike2D,
    xlabels: Optional[Sequence[str]] = None,
    ylabels: Optional[Sequence[str]] = None
):
    """
    Normalize input to (z, x, y) for Heatmap.
    - If DataFrame: rows -> y, columns -> x, values -> z; use index/columns as labels.
    - Else: convert to ndarray; generate default labels if not provided.
    """
    if isinstance(data, pd.DataFrame):
        z = data.values
        x = list(data.columns) if xlabels is None else list(xlabels)
        y = list(data.index) if ylabels is None else list(ylabels)
        return z, x, y

    arr = np.asarray(data)
    if arr.ndim != 2:
        raise ValueError("Heatmap data must be 2D.")
    n_rows, n_cols = arr.shape
    x = list(xlabels) if xlabels is not None else [f"C{j}" for j in range(n_cols)]
    y = list(ylabels) if ylabels is not None else [f"R{i}" for i in range(n_rows)]
    return arr, x, y

def make_heatmap_figure(
    data: ArrayLike2D,
    xlabels: Optional[Sequence[str]] = None,
    ylabels: Optional[Sequence[str]] = None,
    colorscale: str = "Viridis",
    zmin: Optional[float] = None,
    zmax: Optional[float] = None,
    show_annotations: bool = False,
    annotation_format: str = ".2f",  # d3-format for texttemplate when using z values
    custom_text: Optional[Sequence[Sequence[str]]] = None,  # same shape as z if provided
    xgap: int = 0,  # px gap between cells horizontally
    ygap: int = 0,  # px gap between cells vertically
    showscale: bool = True,
    title: Optional[str] = None,
    hovertemplate: Optional[str] = None,  # e.g., "x=%{x}<br>y=%{y}<br>z=%{z:.3f}<extra></extra>"
):
    """
    Build a Plotly Heatmap figure for Dash.

    - data: 2D array-like or DataFrame.
    - xlabels/ylabels: axis labels if data not a DataFrame (or to override DF labels).
    - colorscale: Plotly colorscale name or custom list.
    - zmin/zmax: lock color scale range across multiple heatmaps.
    - show_annotations: overlay text on each cell.
    - annotation_format: d3 format for z text (e.g., '.2f', '.0%', '.3f').
    - custom_text: 2D strings to display instead of z (must match z shape).
    - xgap/ygap: gaps in px to visually separate cells.
    - showscale: show colorbar.
    - title: optional title.
    - hovertemplate: custom hover text template.
    """
    z, x, y = _to_matrix_and_labels(data, xlabels, ylabels)

    # Prepare text + template for annotations if requested
    text = None
    texttemplate = None
    textfont = None

    if show_annotations:
        if custom_text is not None:
            # Use custom text
            text = custom_text
            texttemplate = "%{text}"
        else:
            # Use z values with formatting
            # Note: Heatmap texttemplate supports d3-format via %{z:.2f}, %{z:.0%}, etc.
            texttemplate = f"%{{z:{annotation_format}}}"
        textfont = {"color": "black"}  # ensures visibility on most colorscales

    hm = go.Heatmap(
        z=z,
        x=x,
        y=y,
        colorscale=colorscale,
        zmin=zmin,
        zmax=zmax,
        xgap=xgap,
        ygap=ygap,
        showscale=showscale,
        text=text,
        texttemplate=texttemplate,
        textfont=textfont,
        hovertemplate=hovertemplate,  # if None, default hover is used
        hoverongaps=False,
        colorbar=dict(title="Value")
    )

    fig = go.Figure(data=hm)

    fig.update_layout(
        title=title or "",
        xaxis=dict(title="", tickangle=0, automargin=True),
        yaxis=dict(title="", automargin=True, autorange="reversed"),  # heatmaps often show first row at top
        margin=dict(l=60, r=20, t=60, b=60),
        template="plotly_white"
    )

    return fig
