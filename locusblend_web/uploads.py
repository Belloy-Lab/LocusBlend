"""Uploaded summary-statistic parsing and dataset-sync recommendations.

Backend helpers extracted from app.py. This module never imports Streamlit:
the cached wrappers and all session-state mutation stay in app.py.
"""

import hashlib
import os
import time
from io import BytesIO

import numpy as np
import pandas as pd

from locusblend_web.references import (
    chrom_sort_key,
    format_chrom_label,
    get_supported_chromosomes,
    normalize_chrom,
)
from locusblend_web.variants import _clean_locus_df, chrom_mask


def log(msg):
    """Print a timestamped message; mirrors app.py's log() to avoid a circular import."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def get_uploaded_summary_file_kind(file_name):
    """Return the supported uploaded summary-statistic file kind."""
    if not file_name:
        return None
    name = os.path.basename(str(file_name)).lower()
    for suffix in [".csv.gz", ".tsv.gz", ".txt.gz", ".csv", ".tsv", ".txt"]:
        if name.endswith(suffix):
            return suffix.lstrip(".")
    return None


def _read_uploaded_table_once(file_obj, sep, compression=None, engine=None):
    file_obj.seek(0)
    kwargs = {"sep": sep}
    if compression:
        kwargs["compression"] = compression
    if engine:
        kwargs["engine"] = engine
    return pd.read_csv(file_obj, **kwargs)


def read_summary_stats_file(file_obj):
    """Read an uploaded summary-statistic file in CSV, TSV, TXT, or gzip form."""
    kind = get_uploaded_summary_file_kind(getattr(file_obj, "name", ""))
    unsupported_message = (
        "Unsupported summary-statistic file format. Please upload .csv, .tsv, "
        ".txt, .csv.gz, .tsv.gz, or .txt.gz."
    )
    if kind is None:
        raise ValueError(unsupported_message)

    compression = "gzip" if kind.endswith(".gz") else None

    if kind in {"csv", "csv.gz"}:
        return _read_uploaded_table_once(file_obj, sep=",", compression=compression)

    if kind in {"tsv", "tsv.gz"}:
        return _read_uploaded_table_once(file_obj, sep="\t", compression=compression)

    if kind in {"txt", "txt.gz"}:
        attempts = [
            {"sep": "\t", "engine": None},
            {"sep": r"\s+", "engine": "python"},
            {"sep": ",", "engine": None},
        ]
        errors = []
        for attempt in attempts:
            try:
                df = _read_uploaded_table_once(
                    file_obj,
                    sep=attempt["sep"],
                    compression=compression,
                    engine=attempt["engine"],
                )
                if df.shape[1] >= 2:
                    return df
            except Exception as e:
                errors.append(str(e))
        detail = f" Parser errors: {'; '.join(errors)}" if errors else ""
        raise ValueError(
            "Could not parse TXT summary-statistic file. Please use tab-delimited, "
            f"whitespace-delimited, or comma-delimited text with a header row.{detail}"
        )

    raise ValueError(unsupported_message)


def read_locus_csv(path):
    log(f"loading csv from path: {path}")
    df = pd.read_csv(path)
    print("RAW COLUMNS:", [repr(c) for c in df.columns], flush=True)
    df = _clean_locus_df(df, source_name=path)
    log(f"{path} loaded, shape={df.shape}")
    return df


def read_locus_csv_uploaded(file_bytes, file_name):
    log(f"loading uploaded summary-statistic file: {file_name}")
    raw = BytesIO(file_bytes)
    raw.name = file_name
    df = read_summary_stats_file(raw)
    print("RAW COLUMNS:", [repr(c) for c in df.columns], flush=True)
    df = _clean_locus_df(df, source_name=file_name)
    log(f"{file_name} loaded, shape={df.shape}")
    return df


def make_uploaded_dataset_signature(uploaded_top, uploaded_bottom):
    """Return a stable upload signature, or None when no summary file is uploaded."""
    if uploaded_top is None and uploaded_bottom is None:
        return None

    def file_sig(uploaded):
        if uploaded is None:
            return None
        data = uploaded.getvalue()
        return {
            "name": getattr(uploaded, "name", ""),
            "size": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }

    return {
        "top": file_sig(uploaded_top),
        "bottom": file_sig(uploaded_bottom),
    }


def make_dataset_title_from_filename(file_name):
    """Create a plot title from an uploaded summary-stat filename."""
    if not file_name:
        return ""
    title = os.path.basename(str(file_name)).strip()
    if not title:
        return ""
    lower_title = title.lower()
    for suffix in [".csv.gz", ".tsv.gz", ".txt.gz", ".csv", ".tsv", ".txt", ".gz"]:
        if lower_title.endswith(suffix):
            title = title[: -len(suffix)]
            break
    return title.strip()


def summarize_locus_dataset(df, label):
    """Summarize usable locus rows for upload-driven chromosome/BP inference."""
    summary = {
        "label": label,
        "usable": False,
        "chrom_counts": {},
        "best_rows": {},
        "median_bp": {},
        "warnings": [],
    }
    if df is None or df.empty:
        summary["warnings"].append(f"{label} dataset is empty.")
        return summary
    missing = [c for c in ["CHR", "BP"] if c not in df.columns]
    if missing:
        summary["warnings"].append(f"{label} dataset missing columns: {missing}")
        return summary

    work = df.copy()
    work["CHR_NORM"] = work["CHR"].map(normalize_chrom)
    work["BP_NUM"] = pd.to_numeric(work["BP"], errors="coerce")
    if "P" in work.columns:
        work["P_NUM"] = pd.to_numeric(work["P"], errors="coerce")
    else:
        work["P_NUM"] = np.nan

    work = work[
        work["CHR_NORM"].isin(get_supported_chromosomes())
        & work["BP_NUM"].notna()
    ].copy()
    if work.empty:
        summary["warnings"].append(f"{label} dataset has no valid chr1-22/X BP rows.")
        return summary

    summary["usable"] = True
    summary["chrom_counts"] = {
        str(chrom): int(count)
        for chrom, count in work.groupby("CHR_NORM").size().to_dict().items()
    }

    for chrom, chrom_df in work.groupby("CHR_NORM"):
        valid_p = chrom_df[
            chrom_df["P_NUM"].notna()
            & (chrom_df["P_NUM"] > 0)
            & (chrom_df["P_NUM"] <= 1)
        ].copy()
        if not valid_p.empty:
            best = valid_p.sort_values(["P_NUM", "BP_NUM"], ascending=[True, True]).iloc[0]
            best_p = float(best["P_NUM"])
        else:
            best = chrom_df.sort_values("BP_NUM").iloc[len(chrom_df) // 2]
            best_p = None
        summary["best_rows"][str(chrom)] = {
            "bp": int(best["BP_NUM"]),
            "p": best_p,
        }
        summary["median_bp"][str(chrom)] = int(round(float(chrom_df["BP_NUM"].median())))

    return summary


def choose_sync_chromosome(top_summary, bottom_summary):
    """Choose the chromosome most compatible with uploaded/active datasets."""
    usable_summaries = [s for s in [top_summary, bottom_summary] if s and s.get("usable")]
    if not usable_summaries:
        return None, "none"

    top_chroms = set((top_summary or {}).get("chrom_counts", {}).keys())
    bottom_chroms = set((bottom_summary or {}).get("chrom_counts", {}).keys())
    shared = top_chroms & bottom_chroms
    candidate_chroms = shared if shared else set().union(*(set(s["chrom_counts"].keys()) for s in usable_summaries))

    def score(chrom):
        count = 0
        best_p_values = []
        for summary in usable_summaries:
            count += int(summary["chrom_counts"].get(chrom, 0))
            p = summary.get("best_rows", {}).get(chrom, {}).get("p")
            if p is not None:
                best_p_values.append(float(p))
        best_p = min(best_p_values) if best_p_values else 1.0
        shared_bonus = 1 if chrom in shared else 0
        return (shared_bonus, count, -best_p, -chrom_sort_key(chrom))

    chrom = max(candidate_chroms, key=score)
    if chrom in top_chroms and chrom in bottom_chroms:
        source = "shared"
    elif chrom in top_chroms:
        source = "top"
    else:
        source = "bottom"
    return chrom, source


def choose_sync_center_bp(df_top, df_bottom, chrom):
    """Choose center BP from the best valid P row on chrom, falling back to median BP."""
    candidates = []
    medians = []
    for label, df in [("top", df_top), ("bottom", df_bottom)]:
        if df is None or df.empty or "CHR" not in df.columns or "BP" not in df.columns:
            continue
        work = df.loc[chrom_mask(df, chrom)].copy()
        if work.empty:
            continue
        work["BP_NUM"] = pd.to_numeric(work["BP"], errors="coerce")
        work = work[work["BP_NUM"].notna()].copy()
        if work.empty:
            continue
        medians.append((label, int(round(float(work["BP_NUM"].median())))))
        if "P" in work.columns:
            work["P_NUM"] = pd.to_numeric(work["P"], errors="coerce")
            valid_p = work[
                work["P_NUM"].notna()
                & (work["P_NUM"] > 0)
                & (work["P_NUM"] <= 1)
            ].copy()
            if not valid_p.empty:
                best = valid_p.sort_values(["P_NUM", "BP_NUM"], ascending=[True, True]).iloc[0]
                source_priority = 0 if label == "top" else 1
                candidates.append((float(best["P_NUM"]), source_priority, label, int(best["BP_NUM"])))

    if candidates:
        _, _, source_label, bp = min(candidates)
        return bp, source_label
    if medians:
        source_label, bp = medians[0]
        return bp, source_label
    return None, None


def recommend_y_axis_max_for_dataset(
    df,
    chrom=None,
    center_bp=None,
    window_kb=None,
    min_default=7.0,
    pad_frac=0.12,
):
    """Recommend a y-axis max from -log10(P) for the active uploaded locus.

    Lightweight only. Does not run LD, PLINK, clumping, or figure building.
    """
    try:
        if df is None or df.empty or "P" not in df.columns:
            return float(min_default)

        x = df.copy()

        if chrom is not None and "CHR" in x.columns:
            chrom_str = normalize_chrom(chrom)
            x["_CHR_NORM"] = x["CHR"].map(normalize_chrom)
            x = x[x["_CHR_NORM"] == chrom_str]

        if center_bp is not None and window_kb is not None and "BP" in x.columns:
            bp = pd.to_numeric(x["BP"], errors="coerce")
            center = int(center_bp)
            half_window = int(float(window_kb) * 1000)
            x = x[(bp >= center - half_window) & (bp <= center + half_window)]

        p = pd.to_numeric(x["P"], errors="coerce")
        p = p[np.isfinite(p) & (p > 0)]
        if p.empty:
            return float(min_default)

        y = -np.log10(p)
        y = y[np.isfinite(y)]
        if y.empty:
            return float(min_default)

        ymax = float(np.nanmax(y))
        recommended = max(float(min_default), ymax * (1.0 + float(pad_frac)) + 0.5)

        # Keep a clean display value.
        if recommended <= 20:
            return float(np.ceil(recommended))
        return float(np.ceil(recommended / 5.0) * 5.0)
    except Exception:
        return float(min_default)


def infer_uploaded_locus_sync(df_top, df_bottom, uploaded_top, uploaded_bottom):
    """Infer upload-compatible locus controls without touching LD or plot state."""
    top_uploaded = uploaded_top is not None
    bottom_uploaded = uploaded_bottom is not None
    top_sync_df = df_top if top_uploaded else None
    bottom_sync_df = df_bottom if bottom_uploaded else None
    top_summary = summarize_locus_dataset(top_sync_df, "top") if top_uploaded else None
    bottom_summary = summarize_locus_dataset(bottom_sync_df, "bottom") if bottom_uploaded else None

    chrom, chrom_source = choose_sync_chromosome(top_summary, bottom_summary)
    if chrom is None:
        return {
            "status": "warning",
            "applied": False,
            "message": "Uploaded data were detected, but no valid CHR/BP/P rows could be used for automatic locus sync. Existing locus controls were left unchanged.",
        }

    bp, bp_source = choose_sync_center_bp(top_sync_df, bottom_sync_df, chrom)
    if bp is None:
        return {
            "status": "warning",
            "applied": False,
            "message": "Uploaded data were detected, but no valid CHR/BP/P rows could be used for automatic locus sync. Existing locus controls were left unchanged.",
        }

    preferred_source = "top" if top_uploaded and chrom in top_summary.get("chrom_counts", {}) else "bottom"
    if preferred_source == "bottom" and not bottom_uploaded:
        preferred_source = bp_source or "top"
    index_method = (
        "Auto-select by LD clumping from dataset 1 (top)"
        if preferred_source == "top"
        else "Auto-select by LD clumping from dataset 2 (bottom)"
    )

    no_shared = (
        top_uploaded
        and bottom_uploaded
        and chrom_source != "shared"
    )
    if no_shared:
        message = (
            f"Uploaded datasets do not share a chromosome; synced to {format_chrom_label(chrom)} "
            f"from the {preferred_source} dataset. Check settings before updating."
        )
        status = "warning"
    else:
        message = (
            f"New uploaded data detected. Locus controls synced to {format_chrom_label(chrom)}:{int(bp):,}; "
            f"index selection switched to auto-select from {preferred_source} dataset. "
            "Select Update plot to apply."
        )
        status = "info"

    return {
        "status": status,
        "applied": True,
        "chrom": normalize_chrom(chrom),
        "bp": int(bp),
        "index_selection_method": index_method,
        "preferred_source": preferred_source,
        "message": message,
    }


def format_upload_sync_message(sync_result):
    if not sync_result:
        return ""
    return str(sync_result.get("message", ""))
