"""Gene annotation parsing helpers extracted from app.py.

Pure helpers for parsing GTF attribute fields and for laying out the gene
track rows consumed by the plot layer. They do not touch Streamlit, session
state, the filesystem, or Plotly.
"""

import re

import pandas as pd

from locusblend_web.config import GTF_COLS
from locusblend_web.logutil import log
from locusblend_web.references import normalize_chrom


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


def read_gene_table(parquet_path):
    df = pd.read_parquet(parquet_path)

    if "chrom" in df.columns:
        df["chrom"] = df["chrom"].map(normalize_chrom)
    elif "seqname" in df.columns:
        df["chrom"] = df["seqname"].map(normalize_chrom)

    return df


def filter_genes_from_table(genes, chrom, start, end, gene_display_mode="protein_coding"):
    chrom = normalize_chrom(chrom)

    sub = genes[
        (genes["chrom"] == chrom) &
        (genes["end"] >= start) &
        (genes["start"] <= end)
    ].copy()

    if gene_display_mode == "protein_coding" and "gene_type" in sub.columns:
        sub = sub[sub["gene_type"] == "protein_coding"].copy()

    if "seqname" not in sub.columns:
        sub["seqname"] = "chr" + sub["chrom"].astype(str)

    if "gene_name" not in sub.columns:
        sub["gene_name"] = None
    if "gene_id" not in sub.columns:
        sub["gene_id"] = None
    if "gene_type" not in sub.columns:
        sub["gene_type"] = None
    if "strand" not in sub.columns:
        sub["strand"] = None

    return sub


def read_genes_from_gtf(gtf_path, chrom, start, end, gene_display_mode="protein_coding", chunksize=200000):
    log(f"load_genes_from_gtf start: {gtf_path}, chrom={chrom}, start={start}, end={end}, mode={gene_display_mode}")
    keep = []

    chrom = normalize_chrom(chrom)
    chrom_candidates = [chrom, f"chr{chrom}"]

    for chunk in pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLS,
        compression="infer",
        chunksize=chunksize
    ):
        chunk["seqname"] = chunk["seqname"].astype(str)

        sub = chunk[
            (chunk["feature"] == "gene") &
            (chunk["seqname"].isin(chrom_candidates)) &
            (chunk["end"] >= start) &
            (chunk["start"] <= end)
        ].copy()

        if not sub.empty:
            keep.append(sub)

    if not keep:
        return pd.DataFrame(columns=["seqname", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])

    g = pd.concat(keep, ignore_index=True)
    g["gene_id"] = g["attribute"].apply(lambda x: get_attr(x, "gene_id"))
    g["gene_name"] = g["attribute"].apply(lambda x: get_attr(x, "gene_name"))
    g["gene_type"] = g["attribute"].apply(lambda x: get_attr(x, "gene_type"))

    g = g[["seqname", "start", "end", "strand", "gene_id", "gene_name", "gene_type"]].sort_values("start").reset_index(drop=True)

    if gene_display_mode == "protein_coding":
        g = g[g["gene_type"] == "protein_coding"].copy()

    log(f"load_genes_from_gtf done: n_genes={len(g)}")
    return g
