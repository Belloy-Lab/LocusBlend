"""Gene annotation parsing helpers extracted from app.py.

Pure helpers for parsing GTF attribute fields and for laying out the gene
track rows consumed by the plot layer. They do not touch Streamlit, session
state, the filesystem, or Plotly.
"""

import re

import pandas as pd


def get_attr(attr, key):
    m = re.search(fr'{key} "([^"]+)"', str(attr))
    return m.group(1) if m else None


def assign_gene_rows(df, min_gap=30000):
    if df.empty:
        out = df.copy()
        out["track_row"] = pd.Series(dtype=int)
        return out

    df = df.sort_values("start").copy()
    row_ends = []
    rows = []

    for _, r in df.iterrows():
        placed = False
        for i in range(len(row_ends)):
            if r["start"] > row_ends[i] + min_gap:
                rows.append(i)
                row_ends[i] = r["end"]
                placed = True
                break
        if not placed:
            rows.append(len(row_ends))
            row_ends.append(r["end"])

    df["track_row"] = rows
    return df
