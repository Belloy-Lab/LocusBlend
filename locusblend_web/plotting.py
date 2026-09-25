"""Plotly figure construction for the LocusBlend web app.

Figure builders and figure-modifying helpers extracted from app.py. They never
touch Streamlit, session state, widgets, CSS, or files other than the optional
recombination bigWig. Inputs are already-prepared data frames and static config.
"""

import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pyBigWig
from plotly.subplots import make_subplots

from locusblend_web.config import COLOR_MAPPING, DATA_DIR
from locusblend_web.genes import assign_gene_rows
from locusblend_web.references import normalize_chrom
from locusblend_web.variants import chrom_mask, dedup_columns


def _safe_row_nanmax(arr):
    out = np.full(arr.shape[0], np.nan, dtype=float)
    ok = ~np.all(np.isnan(arr), axis=1)
    if ok.any():
        out[ok] = np.nanmax(arr[ok], axis=1)
    return out


def _pick_bw_chrom(chrom, chroms_dict):
    chrom = str(chrom)
    base = chrom.replace("chr", "")
    candidates = [chrom, base, f"chr{base}"]
    for c in candidates:
        if c in chroms_dict:
            return c
    raise ValueError(f"Cannot find chromosome {chrom} in bigWig.")


def _read_recomb_track_bw(chrom, start, end, bw_path, n_bins=800):
    bw = pyBigWig.open(bw_path)
    try:
        chroms = bw.chroms()
        chrom_name = _pick_bw_chrom(chrom, chroms)
        vals = bw.stats(chrom_name, int(start), int(end), nBins=int(n_bins), type="mean")
    finally:
        bw.close()

    vals = np.array(vals, dtype=float)
    x = np.linspace(start, end, n_bins, endpoint=False) + (end - start) / n_bins / 2

    rec = pd.DataFrame({
        "BP": x,
        "value": vals
    }).dropna()

    return rec


def format_chrom_axis_title(chrom):
    """Return a locus x-axis title with the selected chromosome."""
    chrom = normalize_chrom(chrom)
    return f"Chromosome {chrom} (Mb)"


def add_gene_track_to_subplot(fig, genes_df, row, col=1, min_gap=30000, highlight_names=None):
    if genes_df.empty:
        return fig, 1

    if highlight_names is None:
        highlight_names = set()

    genes_plot = assign_gene_rows(genes_df, min_gap=min_gap)

    for _, r in genes_plot.iterrows():
        y = -int(r["track_row"])
        label = r["gene_name"] if pd.notna(r["gene_name"]) else r["gene_id"]
        is_highlighted = label.lower() in highlight_names

        line_color = "#c45c00" if is_highlighted else "#2f6f7e"
        line_width = 9 if is_highlighted else 6
        font_color = "#c45c00" if is_highlighted else "#2f2f2f"
        font_family = "Arial Black" if is_highlighted else None
        font_size = 11 if is_highlighted else 10

        fig.add_trace(
            go.Scatter(
                x=[r["start"] / 1e6, r["end"] / 1e6],
                y=[y, y],
                mode="lines",
                line=dict(width=line_width, color=line_color),
                hovertemplate=(
                    f"gene: {label}<br>"
                    f"start: {r['start']}<br>"
                    f"end: {r['end']}<br>"
                    f"strand: {r['strand']}<br>"
                    f"type: {r['gene_type']}<extra></extra>"
                ),
                showlegend=False
            ),
            row=row,
            col=col
        )

        _font = dict(size=font_size, color=font_color)
        if font_family:
            _font["family"] = font_family
        fig.add_annotation(
            x=((r["start"] + r["end"]) / 2) / 1e6,
            y=y + 0.18,
            text=label,
            showarrow=False,
            font=_font,
            xanchor="center",
            yanchor="bottom",
            row=row,
            col=col
        )

    n_rows = int(genes_plot["track_row"].max()) + 1
    return fig, n_rows


def get_compare_signal_info(signal, idx1_label, idx2_label, idx3_label):
    mapping = {
        "variant 1": {"r2_col": "r2_1", "color": "#00ffdb", "index_label": idx1_label},
        "variant 2": {"r2_col": "r2_2", "color": "#ff00fa", "index_label": idx2_label},
        "variant 3": {"r2_col": "r2_3", "color": "#ffc900", "index_label": idx3_label},
    }
    return mapping[signal]


def build_compare_data(df_top, df_bottom, signal, idx1_label, idx2_label, idx3_label):
    df_top = dedup_columns(df_top)
    df_bottom = dedup_columns(df_bottom)

    info = get_compare_signal_info(signal, idx1_label, idx2_label, idx3_label)
    r2_col = info["r2_col"]

    top = df_top[["CHR", "BP", "DISPLAY_ID", "P", r2_col, "in_ref", "REF_SNP"]].copy().rename(
        columns={"DISPLAY_ID": "DISPLAY_ID_top", "P": "P_top", r2_col: "r2_use", "in_ref": "in_ref_top", "REF_SNP": "REF_SNP_top"}
    )
    bottom = df_bottom[["CHR", "BP", "DISPLAY_ID", "P", "in_ref", "REF_SNP"]].copy().rename(
        columns={"DISPLAY_ID": "DISPLAY_ID_bottom", "P": "P_bottom", "in_ref": "in_ref_bottom", "REF_SNP": "REF_SNP_bottom"}
    )

    # Prefer REF_SNP for internal alignment; fall back to the position key.
    top["MERGE_KEY"] = np.where(
        top["REF_SNP_top"].notna(),
        top["REF_SNP_top"].astype(str),
        "POS:" + top["CHR"].astype(str) + ":" + top["BP"].astype("Int64").astype(str)
    )
    bottom["MERGE_KEY"] = np.where(
        bottom["REF_SNP_bottom"].notna(),
        bottom["REF_SNP_bottom"].astype(str),
        "POS:" + bottom["CHR"].astype(str) + ":" + bottom["BP"].astype("Int64").astype(str)
    )

    merged = pd.merge(top, bottom, on="MERGE_KEY", how="inner")
    merged["P_top"] = pd.to_numeric(merged["P_top"], errors="coerce")
    merged["P_bottom"] = pd.to_numeric(merged["P_bottom"], errors="coerce")
    merged["r2_use"] = pd.to_numeric(merged["r2_use"], errors="coerce")
    merged["in_ref"] = merged["in_ref_top"].fillna(False).astype(bool)

    merged = merged[
        merged["P_top"].notna() &
        merged["P_bottom"].notna() &
        (merged["P_top"] > 0) &
        (merged["P_bottom"] > 0)
    ].copy()

    merged["logP_top"] = -np.log10(merged["P_top"])
    merged["logP_bottom"] = -np.log10(merged["P_bottom"])
    merged["size"] = np.maximum(8, 4 + 12 * merged["r2_use"].fillna(0))
    return merged, info


def _add_compare_panel(
    fig,
    merged,
    info,
    title_top,
    title_bottom,
    row,
    col,
    x_range,
    y_range,
    ld_labels,
    show_y_title=False
):
    if merged.empty:
        fig.add_annotation(
            x=(x_range[0] + x_range[1]) / 2,
            y=(y_range[0] + y_range[1]) / 2,
            text="No overlapping SNPs",
            showarrow=False,
            font=dict(size=12, color="#666666"),
            row=row,
            col=col
        )
        fig.update_xaxes(
            title_text=f"-log10(P): {title_top}",
            range=x_range,
            zeroline=False,
            showgrid=False,
            row=row,
            col=col
        )
        fig.update_yaxes(
            title_text=f"-log10(P): {title_bottom}" if show_y_title else None,
            range=y_range,
            zeroline=False,
            showgrid=False,
            row=row,
            col=col
        )
        return 0

    missing_ref = merged[~merged["in_ref"]].copy()
    grey = merged[merged["in_ref"] & ((merged["r2_use"] < 0.2) | (merged["r2_use"].isna()))].copy()
    colored = merged[merged["in_ref"] & (merged["r2_use"] >= 0.2)].copy()

    def tooltip_text(r):
        r2_text = "NA" if pd.isna(r["r2_use"]) else f"{r['r2_use']:.3f}"
        ref_text = ld_labels["in_ref"] if bool(r["in_ref"]) else ld_labels["not_in_ref"]
        label = r["DISPLAY_ID_top"] if pd.notna(r["DISPLAY_ID_top"]) else r["MERGE_KEY"]
        return (
            f"SNP: {label}"
            f"<br>-log10(P top): {r['logP_top']:.3f}"
            f"<br>-log10(P bottom): {r['logP_bottom']:.3f}"
            f"<br>r2: {r2_text}"
            f"<br>{ref_text}"
        )

    if not missing_ref.empty:
        missing_ref["tooltip"] = [tooltip_text(r) for _, r in missing_ref.iterrows()]
        fig.add_trace(
            go.Scattergl(
                x=missing_ref["logP_top"],
                y=missing_ref["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="x",
                    color="#d9d9d9",
                    size=4,
                    line=dict(width=0),
                    opacity=0.6,
                ),
                text=missing_ref["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col
        )

    if not grey.empty:
        grey["tooltip"] = [tooltip_text(r) for _, r in grey.iterrows()]
        fig.add_trace(
            go.Scattergl(
                x=grey["logP_top"],
                y=grey["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color="#FFFFFF",
                    line=dict(color="#3c3c3c", width=1),
                    opacity=0.25,
                    size=8,
                ),
                text=grey["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col
        )

    if not colored.empty:
        colored["tooltip"] = [tooltip_text(r) for _, r in colored.iterrows()]
        fig.add_trace(
            go.Scattergl(
                x=colored["logP_top"],
                y=colored["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color=info["color"],
                    line=dict(color="white", width=1),
                    opacity=1.0,
                    size=colored["size"],
                ),
                text=colored["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col
        )

    idx_label = str(info["index_label"]).strip()
    idx_df = merged[merged["DISPLAY_ID_top"].astype(str) == idx_label].copy()
    if idx_df.empty:
        idx_df = merged[merged["DISPLAY_ID_bottom"].astype(str) == idx_label].copy()

    if not idx_df.empty:
        idx_in_ref = bool(idx_df["in_ref"].fillna(False).iloc[0])
        if idx_in_ref:
            fig.add_trace(
                go.Scattergl(
                    x=idx_df["logP_top"],
                    y=idx_df["logP_bottom"],
                    mode="markers",
                    marker=dict(
                        symbol="diamond",
                        color=info["color"],
                        line=dict(color="black", width=2),
                        size=22,
                        opacity=1.0
                    ),
                    text=[f"Index SNP: {idx_label}"] * len(idx_df),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False
                ),
                row=row,
                col=col
            )
        else:
            fig.add_trace(
                go.Scattergl(
                    x=idx_df["logP_top"],
                    y=idx_df["logP_bottom"],
                    mode="markers",
                    marker=dict(
                        symbol="x",
                        color="#d9d9d9",
                        size=4,
                        line=dict(width=0),
                        opacity=0.6,
                    ),
                    text=[f"{ld_labels['index_not_found']}: {idx_label}"] * len(idx_df),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False
                ),
                row=row,
                col=col
            )

    fig.update_xaxes(
        title_text=f"-log10(P): {title_top}",
        range=x_range,
        zeroline=False,
        showgrid=False,
        row=row,
        col=col
    )
    fig.update_yaxes(
        title_text=f"-log10(P): {title_bottom}" if show_y_title else None,
        range=y_range,
        zeroline=False,
        showgrid=False,
        row=row,
        col=col
    )

    return len(merged)


def build_compare_figure_triptych(
    df_top,
    df_bottom,
    idx1_label,
    idx2_label,
    idx3_label,
    title_top,
    title_bottom,
    compare_size,
    ld_labels,
    active_signals=None,
):
    if active_signals is None:
        active_signals = ["variant 1", "variant 2", "variant 3"]
    signals = list(active_signals)
    n_cols = len(signals)
    panel_data = []

    global_xmin = np.inf
    global_xmax = -np.inf
    global_ymin = np.inf
    global_ymax = -np.inf

    for signal in signals:
        merged, info = build_compare_data(
            df_top=df_top,
            df_bottom=df_bottom,
            signal=signal,
            idx1_label=idx1_label,
            idx2_label=idx2_label,
            idx3_label=idx3_label
        )
        panel_data.append((signal, merged, info))

        if not merged.empty:
            global_xmin = min(global_xmin, float(np.floor(np.nanmin(merged["logP_top"]))))
            global_xmax = max(global_xmax, float(np.ceil(np.nanmax(merged["logP_top"]))))
            global_ymin = min(global_ymin, float(np.floor(np.nanmin(merged["logP_bottom"]))))
            global_ymax = max(global_ymax, float(np.ceil(np.nanmax(merged["logP_bottom"]))))

    if not np.isfinite(global_xmin):
        global_xmin, global_xmax = 0.0, 1.0
    if not np.isfinite(global_ymin):
        global_ymin, global_ymax = 0.0, 1.0

    if global_xmax <= global_xmin:
        global_xmax = global_xmin + 1
    if global_ymax <= global_ymin:
        global_ymax = global_ymin + 1

    x_range = [global_xmin, global_xmax]
    y_range = [global_ymin, global_ymax]

    fig = make_subplots(
        rows=1,
        cols=n_cols,
        shared_yaxes=True,
        horizontal_spacing=0.04,
        subplot_titles=[f"Locus compare ({s})" for s in signals]
    )

    counts = {}

    for i, (signal, merged, info) in enumerate(panel_data, start=1):
        counts[signal] = _add_compare_panel(
            fig=fig,
            merged=merged,
            info=info,
            title_top=title_top,
            title_bottom=title_bottom,
            row=1,
            col=i,
            x_range=x_range,
            y_range=y_range,
            ld_labels=ld_labels,
            show_y_title=(i == 1)
        )

    fig.update_layout(
        width=int(compare_size * (1.1 * n_cols)),
        height=int(compare_size),
        dragmode="zoom",
        margin=dict(l=60, r=30, b=60, t=60),
        showlegend=False
    )

    return fig, counts


def build_single_blended_locuscompare(
    df_top,
    df_bottom,
    idx1_ref,
    idx2_ref,
    idx3_ref,
    idx1_label,
    idx2_label,
    idx3_label,
    title_top,
    title_bottom,
    compare_size,
    ld_labels,
    locusblend_mode="Three-index LocusBlend",
):
    df_top = dedup_columns(df_top)
    df_bottom = dedup_columns(df_bottom)

    top_cols = ["CHR", "BP", "DISPLAY_ID", "P", "r2_1", "r2_2", "r2_3", "in_ref", "REF_SNP"]
    top = df_top[[c for c in top_cols if c in df_top.columns]].copy().rename(
        columns={"DISPLAY_ID": "DISPLAY_ID_top", "P": "P_top",
                 "in_ref": "in_ref_top", "REF_SNP": "REF_SNP_top"}
    )
    bottom = df_bottom[["CHR", "BP", "DISPLAY_ID", "P", "in_ref", "REF_SNP"]].copy().rename(
        columns={"DISPLAY_ID": "DISPLAY_ID_bottom", "P": "P_bottom",
                 "in_ref": "in_ref_bottom", "REF_SNP": "REF_SNP_bottom"}
    )

    top["MERGE_KEY"] = np.where(
        top["REF_SNP_top"].notna(),
        top["REF_SNP_top"].astype(str),
        "POS:" + top["CHR"].astype(str) + ":" + top["BP"].astype("Int64").astype(str)
    )
    bottom["MERGE_KEY"] = np.where(
        bottom["REF_SNP_bottom"].notna(),
        bottom["REF_SNP_bottom"].astype(str),
        "POS:" + bottom["CHR"].astype(str) + ":" + bottom["BP"].astype("Int64").astype(str)
    )

    merged = pd.merge(top, bottom, on="MERGE_KEY", how="inner", suffixes=("", "_b"))
    merged = merged.drop_duplicates("MERGE_KEY").copy()

    merged["P_top"] = pd.to_numeric(merged["P_top"], errors="coerce")
    merged["P_bottom"] = pd.to_numeric(merged["P_bottom"], errors="coerce")
    for c in ["r2_1", "r2_2", "r2_3"]:
        if c in merged.columns:
            merged[c] = pd.to_numeric(merged[c], errors="coerce")
        else:
            merged[c] = np.nan
    merged["in_ref"] = merged["in_ref_top"].fillna(False).astype(bool)

    merged = merged[
        merged["P_top"].notna() &
        merged["P_bottom"].notna() &
        (merged["P_top"] > 0) &
        (merged["P_bottom"] > 0)
    ].copy()

    if merged.empty:
        fig = go.Figure()
        fig.update_layout(
            width=compare_size,
            height=compare_size,
            title=dict(text="Locus compare (blended)", x=0.5, xanchor="center")
        )
        return fig, 0

    merged["logP_top"] = -np.log10(merged["P_top"])
    merged["logP_bottom"] = -np.log10(merged["P_bottom"])

    bins = [-np.inf, 0.2, 0.4, 0.6, 0.8, np.inf]
    labels = ["0", "2", "4", "6", "8"]
    merged["ld_bin"] = pd.cut(merged["r2_1"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")

    if locusblend_mode == "Standard locus zoom":
        merged["r2_2_bin"] = "0"
        merged["r2_3_bin"] = "0"
        merged["group_code"] = merged["ld_bin"] + "0" + "0"
    elif locusblend_mode == "Two-index LocusBlend":
        merged["r2_2_bin"] = pd.cut(merged["r2_2"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        merged["r2_3_bin"] = "0"
        merged["group_code"] = merged["ld_bin"] + merged["r2_2_bin"] + "0"
    else:
        merged["r2_2_bin"] = pd.cut(merged["r2_2"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        merged["r2_3_bin"] = pd.cut(merged["r2_3"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        merged["group_code"] = merged["ld_bin"] + merged["r2_2_bin"] + merged["r2_3_bin"]

    if locusblend_mode == "Standard locus zoom":
        arr_any = merged[["r2_1"]].to_numpy(dtype=float)
    elif locusblend_mode == "Two-index LocusBlend":
        arr_any = merged[["r2_1", "r2_2"]].to_numpy(dtype=float)
    else:
        arr_any = merged[["r2_1", "r2_2", "r2_3"]].to_numpy(dtype=float)
    merged["max_r_any"] = _safe_row_nanmax(arr_any)

    index_keys = {
        str(x).strip()
        for x in [idx1_ref, idx2_ref, idx3_ref, idx1_label, idx2_label, idx3_label]
        if x is not None and str(x).strip() != ""
    }

    is_index = pd.Series(False, index=merged.index)
    for col in ["REF_SNP_top", "REF_SNP_bottom", "DISPLAY_ID_top", "DISPLAY_ID_bottom"]:
        if col in merged.columns:
            is_index = is_index | merged[col].astype(str).str.strip().isin(index_keys)
    merged["is_index"] = is_index

    arr_size = arr_any.copy()
    idx_mask = merged["is_index"].to_numpy()
    if idx_mask.any():
        idx_vals = arr_size[idx_mask]
        idx_vals[np.isclose(idx_vals, 1.0, equal_nan=False)] = np.nan
        arr_size[idx_mask] = idx_vals
    merged["max_r_size"] = _safe_row_nanmax(arr_size)
    merged["max_r_size"] = merged["max_r_size"].fillna(0)

    colored = merged[merged["in_ref"] & (merged["max_r_any"] >= 0.2)].copy()
    grey = merged[merged["in_ref"] & ((merged["max_r_any"] < 0.2) | (merged["max_r_any"].isna()))].copy()
    missing_ref = merged[~merged["in_ref"]].copy()

    if not colored.empty:
        colored["fill_hex"] = colored["group_code"].map(COLOR_MAPPING).fillna("#bfbfbf")
        colored["size"] = np.clip(2 + 10 * colored["max_r_size"].fillna(0), 4, 12)
    if not grey.empty:
        grey["size"] = 4
    if not missing_ref.empty:
        missing_ref["size"] = 6

    def tooltip_text(r):
        label = r["DISPLAY_ID_top"] if pd.notna(r.get("DISPLAY_ID_top")) else r["MERGE_KEY"]
        ref_text = ld_labels["in_ref"] if bool(r["in_ref"]) else ld_labels["not_in_ref"]
        parts = [
            f"SNP: {label}",
            f"-log10(P top): {r['logP_top']:.3f}",
            f"-log10(P bottom): {r['logP_bottom']:.3f}",
        ]
        for k in ["r2_1", "r2_2", "r2_3"]:
            v = r.get(k)
            parts.append(f"{k}: {'NA' if pd.isna(v) else f'{v:.3f}'}")
        parts.append(ref_text)
        return "<br>".join(parts)

    if not missing_ref.empty:
        missing_ref["tooltip"] = [tooltip_text(r) for _, r in missing_ref.iterrows()]
    if not grey.empty:
        grey["tooltip"] = [tooltip_text(r) for _, r in grey.iterrows()]
    if not colored.empty:
        colored["tooltip"] = [tooltip_text(r) for _, r in colored.iterrows()]

    fig = go.Figure()

    if not missing_ref.empty:
        fig.add_trace(
            go.Scattergl(
                x=missing_ref["logP_top"],
                y=missing_ref["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="x",
                    color="#d9d9d9",
                    size=4,
                    line=dict(width=0),
                    opacity=0.6,
                ),
                text=missing_ref["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

    if not grey.empty:
        fig.add_trace(
            go.Scattergl(
                x=grey["logP_top"],
                y=grey["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color="#bfbfbf",
                    size=grey["size"],
                    line=dict(width=0),
                    opacity=0.45,
                ),
                text=grey["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

    if not colored.empty:
        fig.add_trace(
            go.Scattergl(
                x=colored["logP_top"],
                y=colored["logP_bottom"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color=colored["fill_hex"],
                    line=dict(color="white", width=1),
                    opacity=1.0,
                    size=colored["size"],
                ),
                text=colored["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            )
        )

    def add_blended_index_marker(fig, idx_ref, idx_label, fill_color):
        if idx_ref is None and (idx_label is None or str(idx_label).strip() == ""):
            return fig

        idx_df = pd.DataFrame()
        if idx_ref is not None:
            for col in ["REF_SNP_top", "REF_SNP_bottom"]:
                if col in merged.columns:
                    cand = merged[merged[col].astype(str) == str(idx_ref)]
                    if not cand.empty:
                        idx_df = cand.copy()
                        break

        if idx_df.empty and idx_label is not None and str(idx_label).strip() != "":
            for col in ["DISPLAY_ID_top", "DISPLAY_ID_bottom"]:
                if col in merged.columns:
                    cand = merged[merged[col].astype(str) == str(idx_label).strip()]
                    if not cand.empty:
                        idx_df = cand.copy()
                        break

        if idx_df.empty:
            return fig

        idx_in_ref = bool(idx_df["in_ref"].fillna(False).iloc[0])
        if idx_in_ref:
            fig.add_trace(
                go.Scattergl(
                    x=idx_df["logP_top"],
                    y=idx_df["logP_bottom"],
                    mode="markers",
                    marker=dict(
                        symbol="diamond",
                        color=fill_color,
                        line=dict(color="black", width=2),
                        size=18,
                        opacity=1.0,
                    ),
                    text=[f"Index SNP: {idx_label}"] * len(idx_df),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                )
            )
        else:
            fig.add_trace(
                go.Scattergl(
                    x=idx_df["logP_top"],
                    y=idx_df["logP_bottom"],
                    mode="markers",
                    marker=dict(
                        symbol="x",
                        color="#d9d9d9",
                        size=4,
                        line=dict(width=0),
                        opacity=0.6,
                    ),
                    text=[f"{ld_labels['index_not_found']}: {idx_label}"] * len(idx_df),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                )
            )
        return fig

    fig = add_blended_index_marker(fig, idx1_ref, idx1_label, "#00ffdb")
    if locusblend_mode != "Standard locus zoom":
        fig = add_blended_index_marker(fig, idx2_ref, idx2_label, "#ff00fa")
    if locusblend_mode == "Three-index LocusBlend":
        fig = add_blended_index_marker(fig, idx3_ref, idx3_label, "#ffc900")

    x_min = float(np.floor(np.nanmin(merged["logP_top"])))
    x_max = float(np.ceil(np.nanmax(merged["logP_top"])))
    y_min = float(np.floor(np.nanmin(merged["logP_bottom"])))
    y_max = float(np.ceil(np.nanmax(merged["logP_bottom"])))
    if x_max <= x_min:
        x_max = x_min + 1
    if y_max <= y_min:
        y_max = y_min + 1

    fig.update_xaxes(
        title_text=f"-log10(P): {title_top}",
        range=[x_min, x_max],
        zeroline=False,
        showgrid=False
    )
    fig.update_yaxes(
        title_text=f"-log10(P): {title_bottom}",
        range=[y_min, y_max],
        zeroline=False,
        showgrid=False
    )
    fig.update_layout(
        title=dict(text="Locus compare (blended)", x=0.5, xanchor="center"),
        width=int(compare_size),
        height=int(compare_size),
        autosize=False,
        dragmode="zoom",
        margin=dict(l=70, r=40, b=70, t=50),
        showlegend=False
    )

    return fig, len(merged)


def make_locus_tooltip(row, ref_status):
    label = row["DISPLAY_ID"] if "DISPLAY_ID" in row.index and pd.notna(row["DISPLAY_ID"]) else row.get("REF_SNP", row.get("SNP", "NA"))
    parts = [
        f"SNP: {label}",
        f"BP: {int(row['BP'])}",
        f"-log10(P): {row['logP']:.3f}",
    ]

    if "A1" in row.index and "A2" in row.index:
        parts.append(f"coding: {row['A1']}/{row['A2']}")
    if "beta" in row.index:
        parts.append(f"beta: {row['beta']:.5g}" if pd.notna(row["beta"]) else "beta: NA")
    if "coding_flipped" in row.index:
        parts.append(f"coding_flipped: {bool(row['coding_flipped'])}")
    if "r2_1" in row.index:
        parts.append(f"r2_1: {row['r2_1']:.3f}" if pd.notna(row["r2_1"]) else "r2_1: NA")
    if "r2_2" in row.index:
        parts.append(f"r2_2: {row['r2_2']:.3f}" if pd.notna(row["r2_2"]) else "r2_2: NA")
    if "r2_3" in row.index:
        parts.append(f"r2_3: {row['r2_3']:.3f}" if pd.notna(row["r2_3"]) else "r2_3: NA")
    if "group_code" in row.index:
        parts.append(f"code: {row['group_code']}")

    parts.append(ref_status)
    return "<br>".join(parts)


def get_plotly_locus_py(
    max_ylim,
    bp,
    window_bp,
    merged_female_withld,
    merged_df,
    idx1_ref,
    idx2_ref,
    idx3_ref,
    idx1_label,
    idx2_label,
    idx3_label,
    y_label,
    ld_labels,
    chrom,
    show_recomb=True,
    bw_path=str(DATA_DIR / "recomb1000GAvg.bw"),
    n_recomb_bins=800,
    locusblend_mode="Three-index LocusBlend",
):
    chrom = normalize_chrom(chrom)

    merged_female_withld = dedup_columns(merged_female_withld)
    merged_df = dedup_columns(merged_df)

    # Subset to selected chromosome first
    mf = merged_female_withld.loc[chrom_mask(merged_female_withld, chrom)].copy()
    md = merged_df.loc[chrom_mask(merged_df, chrom)].copy()

    if mf.empty or md.empty:
        raise ValueError(f"No SNPs found on chromosome {chrom} in the selected dataset.")

    data_bp_min = int(np.nanmin(mf["BP"]))
    data_bp_max = int(np.nanmax(mf["BP"]))

    bp_start = max(data_bp_min, int(bp - window_bp))
    bp_end = min(data_bp_max, int(bp + window_bp))

    d = md.loc[
        (md["BP"] >= bp_start) & (md["BP"] <= bp_end)
    ].copy()
    d = dedup_columns(d)

    if d.empty:
        raise ValueError(f"No SNPs found in chr{chrom}:{bp_start}-{bp_end}.")

    for col in ["BP", "P", "r2_1", "r2_2", "r2_3"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    d["CHR"] = d["CHR"].astype(str)
    d["logP"] = -np.log10(d["P"].clip(lower=np.finfo(float).tiny))
    d["in_ref"] = d["in_ref"].fillna(False).astype(bool)

    bins = [-np.inf, 0.2, 0.4, 0.6, 0.8, np.inf]
    labels = ["0", "2", "4", "6", "8"]

    d["ld_bin"] = pd.cut(d["r2_1"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")

    if locusblend_mode == "Standard locus zoom":
        d["r2_2_bin"] = "0"
        d["r2_3_bin"] = "0"
        d["group_code"] = d["ld_bin"] + "0" + "0"
    elif locusblend_mode == "Two-index LocusBlend":
        d["r2_2_bin"] = pd.cut(d["r2_2"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        d["r2_3_bin"] = "0"
        d["group_code"] = d["ld_bin"] + d["r2_2_bin"] + "0"
    else:
        d["r2_2_bin"] = pd.cut(d["r2_2"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        d["r2_3_bin"] = pd.cut(d["r2_3"], bins=bins, labels=labels, include_lowest=True).astype("string").fillna("0")
        d["group_code"] = d["ld_bin"] + d["r2_2_bin"] + d["r2_3_bin"]

    # ---------- size logic ----------
    # max_r_any: still used for color / grey split
    # max_r_size: used only for marker size
    # only remove self-LD=1 for true index SNP rows
    if locusblend_mode == "Standard locus zoom":
        arr_any = d[["r2_1"]].to_numpy(dtype=float)
    elif locusblend_mode == "Two-index LocusBlend":
        arr_any = d[["r2_1", "r2_2"]].to_numpy(dtype=float)
    else:
        arr_any = d[["r2_1", "r2_2", "r2_3"]].to_numpy(dtype=float)
    d["max_r_any"] = _safe_row_nanmax(arr_any)

    index_keys = {
        str(x).strip()
        for x in [idx1_ref, idx2_ref, idx3_ref, idx1_label, idx2_label, idx3_label]
        if x is not None and str(x).strip() != ""
    }

    d["is_index"] = False
    for col in ["REF_SNP", "DISPLAY_ID", "rsid", "SNP", "REF_MATCH", "uniqueid", "QUERY_UID"]:
        if col in d.columns:
            d["is_index"] = d["is_index"] | d[col].astype(str).str.strip().isin(index_keys)

    arr_size = arr_any.copy()
    idx_mask = d["is_index"].to_numpy()

    if idx_mask.any():
        idx_vals = arr_size[idx_mask]
        idx_vals[np.isclose(idx_vals, 1.0, equal_nan=False)] = np.nan
        arr_size[idx_mask] = idx_vals

    d["max_r_size"] = _safe_row_nanmax(arr_size)
    d["max_r_size"] = d["max_r_size"].fillna(0)

    colored = d.loc[d["in_ref"] & (d["max_r_any"] >= 0.2)].copy()
    grey = d.loc[d["in_ref"] & ((d["max_r_any"] < 0.2) | (d["max_r_any"].isna()))].copy()
    missing_ref = d.loc[~d["in_ref"]].copy()

    for x in [colored, grey, missing_ref]:
        x["x_mb"] = x["BP"] / 1e6

    colored["fill_hex"] = colored["group_code"].map(COLOR_MAPPING).fillna("#bfbfbf")

    # Marker size for regular colored points.
    # Non-index SNPs in full LD with an index (r2=1) still grow normally.
    # Index SNPs do not grow oversized from their own self-LD=1.
    colored["size"] = np.clip(2 + 10 * colored["max_r_size"].fillna(0), 4, 12)

    grey["size"] = 4
    missing_ref["size"] = 6

    if not colored.empty:
        colored["tooltip"] = [make_locus_tooltip(row, ld_labels["in_ref"]) for _, row in colored.iterrows()]
    if not grey.empty:
        grey["tooltip"] = [make_locus_tooltip(row, ld_labels["in_ref"]) for _, row in grey.iterrows()]
    if not missing_ref.empty:
        missing_ref["tooltip"] = [make_locus_tooltip(row, ld_labels["not_in_ref"]) for _, row in missing_ref.iterrows()]

    min_ylim = int(np.floor(np.nanmin(d["logP"])))
    if max_ylim is None:
        max_ylim = int(np.ceil(np.nanmax(d["logP"])))

    if show_recomb and os.path.exists(bw_path):
        try:
            rec = _read_recomb_track_bw(chrom, bp_start, bp_end, bw_path=bw_path, n_bins=n_recomb_bins)
            rec["x_mb"] = rec["BP"] / 1e6
        except Exception:
            rec = pd.DataFrame(columns=["BP", "value", "x_mb"])
    else:
        rec = pd.DataFrame(columns=["BP", "value", "x_mb"])

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # 1) missing_ref (x) first (bottom layer)
    if not missing_ref.empty:
        fig.add_trace(
            go.Scattergl(
                x=missing_ref["x_mb"],
                y=missing_ref["logP"],
                mode="markers",
                marker=dict(
                    symbol="x",
                    color="#d9d9d9",
                    size=4,
                    line=dict(width=0),
                    opacity=0.6,
                ),
                text=missing_ref["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            secondary_y=False,
        )

    # 2) grey points
    if not grey.empty:
        fig.add_trace(
            go.Scattergl(
                x=grey["x_mb"],
                y=grey["logP"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color="#bfbfbf",
                    size=grey["size"],
                    line=dict(width=0),
                    opacity=0.45,
                ),
                text=grey["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            secondary_y=False,
        )

    # 3) colored points
    if not colored.empty:
        fig.add_trace(
            go.Scattergl(
                x=colored["x_mb"],
                y=colored["logP"],
                mode="markers",
                marker=dict(
                    symbol="circle",
                    color=colored["fill_hex"],
                    line=dict(color="white", width=1),
                    opacity=1.0,
                    size=colored["size"],
                ),
                text=colored["tooltip"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ),
            secondary_y=False,
        )

    def add_index_marker(fig, idx_ref, idx_label, fill_color):
        idx = pd.DataFrame()
        if idx_ref is not None:
            idx = d.loc[d["REF_SNP"].astype(str) == str(idx_ref)].copy()

        if idx.empty:
            idx = d.loc[d["DISPLAY_ID"].astype(str) == str(idx_label)].copy()

        if idx.empty:
            return fig

        idx_in_ref = bool(idx["in_ref"].fillna(False).iloc[0])

        if idx_in_ref:
            fig.add_trace(
                go.Scattergl(
                    x=idx["BP"] / 1e6,
                    y=idx["logP"],
                    mode="markers",
                    marker=dict(
                        symbol="diamond",
                        color=fill_color,
                        line=dict(color="black", width=1),
                        size=14,
                        opacity=0.9,
                    ),
                    text=[f"Index SNP: {idx_label}"] * len(idx),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                ),
                secondary_y=False,
            )
        else:
            fig.add_trace(
                go.Scattergl(
                    x=idx["BP"] / 1e6,
                    y=idx["logP"],
                    mode="markers",
                    marker=dict(
                        symbol="x",
                        color="#d9d9d9",
                        size=4,
                        line=dict(width=0),
                        opacity=0.6,
                    ),
                    text=[f"{ld_labels['index_not_found']}: {idx_label}"] * len(idx),
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                ),
                secondary_y=False,
            )
        return fig

    fig = add_index_marker(fig, idx1_ref, idx1_label, "#00ffdb")
    if locusblend_mode != "Standard locus zoom":
        fig = add_index_marker(fig, idx2_ref, idx2_label, "#ff00fa")
    if locusblend_mode == "Three-index LocusBlend":
        fig = add_index_marker(fig, idx3_ref, idx3_label, "#ffc900")

    if not rec.empty:
        fig.add_trace(
            go.Scattergl(
                x=rec["x_mb"],
                y=rec["value"],
                mode="lines",
                line=dict(color="#2d2d2d", width=1),
                hoverinfo="skip",
                showlegend=False,
            ),
            secondary_y=True,
        )

    fig.update_layout(
        title=dict(text=y_label, x=0.5, xanchor="center"),
        dragmode="zoom",
        margin=dict(l=60, r=60, b=50, t=40),
    )

    fig.update_xaxes(title_text=format_chrom_axis_title(chrom), zeroline=False)
    fig.update_yaxes(
        title_text="-log<sub>10</sub>(P)",
        range=[min_ylim, max_ylim],
        zeroline=False,
        secondary_y=False
    )
    fig.update_yaxes(
        title_text="Recombination rate",
        range=[0, 100],
        autorange=False,
        zeroline=False,
        showgrid=False,
        secondary_y=True
    )

    return fig, len(d), bp_start, bp_end


def apply_locusblend_plot_theme(fig):
    """Force Plotly figures to remain readable in Streamlit/browser dark mode.

    Layout-only. Must not change trace data, marker colors, LD color mapping,
    marker sizes, group_code, or index-marker colors.
    """
    if fig is None:
        return fig

    text = "#111827"
    muted = "#374151"
    grid = "#e5e7eb"
    axis = "#9ca3af"
    bg = "#ffffff"

    fig.update_layout(
        template="plotly_white",
        paper_bgcolor=bg,
        plot_bgcolor=bg,
        font=dict(color=text),
        hoverlabel=dict(
            bgcolor=bg,
            bordercolor=axis,
            font=dict(color=text),
        ),
    )

    # Do not overwrite fig.layout.title or title.text.
    # Only set title font color if a title already exists.
    try:
        if getattr(fig.layout, "title", None) is not None:
            existing_title = getattr(fig.layout.title, "text", None)
            if existing_title not in (None, "", "undefined"):
                fig.update_layout(title_font_color=text)
            elif existing_title == "undefined":
                fig.update_layout(title_text=None)
    except Exception:
        pass

    fig.update_xaxes(
        title_font=dict(color=text),
        tickfont=dict(color=muted),
        gridcolor=grid,
        zerolinecolor=grid,
        linecolor=axis,
    )

    fig.update_yaxes(
        title_font=dict(color=text),
        tickfont=dict(color=muted),
        gridcolor=grid,
        zerolinecolor=grid,
        linecolor=axis,
    )

    # Preserve explicit annotation colors, especially highlighted gene labels.
    # Only fill in missing annotation font color.
    for ann in fig.layout.annotations or []:
        try:
            if ann.font is None:
                ann.font = dict(color=text)
            elif not getattr(ann.font, "color", None):
                ann.font.color = text
        except Exception:
            pass

    return fig


def _numeric_values(values):
    """Return finite numeric values from a Plotly trace coordinate array."""
    if values is None:
        return []
    try:
        arr = pd.to_numeric(pd.Series(list(values)), errors="coerce")
        arr = arr[np.isfinite(arr)]
        return arr.astype(float).tolist()
    except Exception:
        return []


def _padded_range(values, include_zero=True, pad_frac=0.12):
    """Compute a safe padded numeric range for Plotly axes.

    Similar to ggplot2 expand: add visual breathing room around the data so
    large diamond/index markers are not clipped by axis limits.
    """
    vals = [float(v) for v in values if np.isfinite(v)]
    if not vals:
        return None

    vmin = min(vals)
    vmax = max(vals)

    if include_zero:
        vmin = min(0.0, vmin)

    span = max(vmax - vmin, 1.0)
    pad = max(span * float(pad_frac), 0.35)

    lower = vmin - pad
    upper = vmax + pad

    if include_zero:
        lower = min(0.0, lower)

    return [lower, upper]


def apply_locus_compare_safe_autoscale(fig, pad_frac=0.12):
    """Set a safe initial axis range for locus compare figures.

    This prevents Plotly 'Reset axes' from returning to a too-tight or
    stale range that clips points or diamond markers.

    Layout-only. Must not change trace data, marker colors, LD colors,
    marker sizes, group_code, or index-marker colors.
    """
    if fig is None:
        return fig

    all_x = []
    all_y = []

    for trace in fig.data:
        all_x.extend(_numeric_values(getattr(trace, "x", None)))
        all_y.extend(_numeric_values(getattr(trace, "y", None)))

    x_range = _padded_range(all_x, include_zero=True, pad_frac=pad_frac)
    y_range = _padded_range(all_y, include_zero=True, pad_frac=pad_frac)

    if x_range is not None:
        fig.update_xaxes(range=x_range, autorange=False)

    if y_range is not None:
        fig.update_yaxes(range=y_range, autorange=False)

    return fig


def get_locus_compare_plotly_config(base_config=None):
    """Return Plotly config for locus compare charts.

    Removes Reset axes because safe autoscale is the intended recovery
    behavior for compare plots. Autoscale should remain available.
    """
    cfg = dict(base_config or {})
    remove = list(cfg.get("modeBarButtonsToRemove", []))

    for button_name in ["resetScale2d"]:
        if button_name not in remove:
            remove.append(button_name)

    cfg["modeBarButtonsToRemove"] = remove
    return cfg
