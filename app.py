import base64
from io import BytesIO
import hashlib
import os
import re
import time
import traceback
import html as html_lib
from mimetypes import guess_type

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import streamlit.components.v1 as components

from locusblend_web.config import (
    ASSET_DIR,
    BIN_DIR,
    COLOR_MAPPING,
    DATA_DIR,
    INTERNAL_1000G_ANCESTRIES,
    INTERNAL_1000G_ANCESTRY_OPTIONS,
    INTERNAL_1000G_DEFAULT_ANCESTRY,
)

from locusblend_web.references import (
    get_gtf_path_for_chrom,
    get_internal_1000g_prefix,
    normalize_chrom,
    normalize_internal_1000g_ancestry,
)

from locusblend_web.variants import (
    _clean_locus_df,
    attach_reference_snp_two_pass,
    attach_uploaded_ld_keys,
    chrom_mask,
    dedup_columns,
    resolve_index_variant_from_input,
)

from locusblend_web.genes import (
    assign_gene_rows,
    filter_genes_from_table,
    get_attr,
    read_gene_table,
    read_genes_from_gtf,
)

from locusblend_web.ld import (
    build_ld_annot_for_window,
    build_ld_maps_with_plink,
    compute_ld_maps_from_uploaded_long,
    compute_ld_maps_from_uploaded_matrix,
    greedy_clump_uploaded_ld_for_auto_indices,
    merge_ld_annot,
    parse_uploaded_ld_long,
    parse_uploaded_ld_matrix,
    prepare_clump_candidates,
    read_reference_bim,
    run_plink_clump_for_auto_indices,
)

from locusblend_web.plotting import (
    add_gene_track_to_subplot,
    apply_locus_compare_safe_autoscale,
    apply_locusblend_plot_theme,
    build_compare_figure_triptych,
    build_single_blended_locuscompare,
    format_chrom_axis_title,
    get_plotly_locus_py,
)

from locusblend_web.export import (
    build_locusblend_export_image,
    clone_plotly_figure,
)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


@st.cache_data(show_spinner=False)
def load_locus_csv(path):
    log(f"loading csv from path: {path}")
    df = pd.read_csv(path)
    print("RAW COLUMNS:", [repr(c) for c in df.columns], flush=True)
    df = _clean_locus_df(df, source_name=path)
    log(f"{path} loaded, shape={df.shape}")
    return df


@st.cache_data(show_spinner=False)
def load_locus_csv_uploaded(file_bytes, file_name):
    log(f"loading uploaded summary-statistic file: {file_name}")
    raw = BytesIO(file_bytes)
    raw.name = file_name
    df = read_summary_stats_file(raw)
    print("RAW COLUMNS:", [repr(c) for c in df.columns], flush=True)
    df = _clean_locus_df(df, source_name=file_name)
    log(f"{file_name} loaded, shape={df.shape}")
    return df


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


def apply_uploaded_dataset_sync_if_new(
    signature,
    sync_result,
    df_top=None,
    df_bottom=None,
    uploaded_top=None,
    uploaded_bottom=None,
):
    """Apply upload-inferred state once per new uploaded dataset signature."""
    if signature is None:
        return False
    if st.session_state.get("last_uploaded_dataset_signature") == signature:
        return False

    if sync_result.get("applied"):
        st.session_state["chrom"] = sync_result["chrom"]
        st.session_state["bp"] = int(sync_result["bp"])
        st.session_state["index_selection_method"] = sync_result["index_selection_method"]
        if uploaded_top is not None:
            st.session_state["title_top"] = make_dataset_title_from_filename(uploaded_top.name)
        if uploaded_bottom is not None:
            st.session_state["title_bottom"] = make_dataset_title_from_filename(uploaded_bottom.name)

        synced_chrom = sync_result["chrom"]
        synced_bp = int(sync_result["bp"])
        synced_window_kb = st.session_state.get("window_kb", 500)

        st.session_state["max_ylim_top"] = recommend_y_axis_max_for_dataset(
            df_top,
            chrom=synced_chrom,
            center_bp=synced_bp,
            window_kb=synced_window_kb,
            min_default=7.0,
            pad_frac=0.12,
        )
        st.session_state["max_ylim_bottom"] = recommend_y_axis_max_for_dataset(
            df_bottom,
            chrom=synced_chrom,
            center_bp=synced_bp,
            window_kb=synced_window_kb,
            min_default=7.0,
            pad_frac=0.12,
        )

        st.session_state["rsid"] = ""
        st.session_state["rsid_2"] = ""
        st.session_state["rsid_3"] = ""

    st.session_state["last_uploaded_dataset_signature"] = signature
    st.session_state["last_upload_sync_result"] = sync_result
    return True


def format_upload_sync_message(sync_result):
    if not sync_result:
        return ""
    return str(sync_result.get("message", ""))


@st.cache_data(show_spinner=False)
def read_uploaded_ld_long(file_bytes, file_name):
    return parse_uploaded_ld_long(file_bytes, file_name)


@st.cache_data(show_spinner=False)
def read_uploaded_ld_matrix(file_bytes, file_name):
    return parse_uploaded_ld_matrix(file_bytes, file_name)


def is_supported_chrom(chrom):
    """Return True for chromosomes supported by the visualizer."""
    return normalize_chrom(chrom) in get_supported_chromosomes()


def get_supported_chromosomes():
    """Return supported autosomes plus chromosome X."""
    return [str(i) for i in range(1, 23)] + ["X"]


def chrom_sort_key(chrom):
    """Sort chromosomes 1-22 numerically, then X."""
    c = normalize_chrom(chrom)
    if c == "X":
        return 23
    try:
        return int(c)
    except (TypeError, ValueError):
        return 10_000


def format_chrom_label(chrom):
    """Return display label such as chr14 or chrX."""
    c = normalize_chrom(chrom)
    return f"chr{c}" if c else "chr"


def format_internal_1000g_ancestry_option(ancestry):
    ancestry = normalize_internal_1000g_ancestry(ancestry)
    label = INTERNAL_1000G_ANCESTRIES[ancestry].split(" / ", 1)[0]
    return f"{ancestry} - {label}"


def get_internal_bfile_prefix_for_chrom(chrom, ancestry=INTERNAL_1000G_DEFAULT_ANCESTRY):
    """Build the chromosome-specific PLINK bfile prefix and verify all three
    files (.bed, .bim, .fam) exist."""
    chrom = normalize_chrom(chrom)
    ancestry = normalize_internal_1000g_ancestry(ancestry)
    prefix = str(get_internal_1000g_prefix(chrom, ancestry))
    missing = []
    for suf in [".bed", ".bim", ".fam"]:
        if not os.path.exists(prefix + suf):
            missing.append(prefix + suf)
    if missing:
        if chrom == "X":
            message = (
                f"Internal 1000G {ancestry} reference files for chromosome X were not found. "
                "Select an ancestry with chrX support or upload your own LD reference."
            )
        else:
            message = f"Internal 1000G {ancestry} reference files for chromosome {chrom} were not found."
        raise FileNotFoundError(
            message + "\n" +
            "\n".join(f"  {m}" for m in missing)
        )
    return prefix


def get_required_n_indices(active_mode):
    """Return the number of index variants required by the active mode."""
    if active_mode == "Standard locus zoom":
        return 1
    if active_mode == "Two-index LocusBlend":
        return 2
    return 3


def auto_select_index_variants_by_clumping(
    df_source_ref,
    selected_chrom,
    center_bp,
    window_bp,
    required_n,
    active_ld_source,
    bfile_prefix=None,
    plink_path=None,
    ld_long_table=None,
    ld_matrix_table=None,
    clump_r2=0.01,
):
    """Dispatch to PLINK or greedy clumping, returning (selected_rows, summary_dict)."""
    candidates = prepare_clump_candidates(
        df_source_ref, selected_chrom, center_bp, window_bp,
    )
    if candidates.empty:
        raise ValueError(
            f"No valid candidate SNPs found in chr{selected_chrom}:"
            f"{center_bp - window_bp}-{center_bp + window_bp}. "
            "Check the chromosome selection, window size, and that the "
            "dataset has been matched to an LD reference."
        )

    start_bp = int(center_bp - window_bp)
    end_bp = int(center_bp + window_bp)

    if active_ld_source == "Use internal 1000G reference":
        selected = run_plink_clump_for_auto_indices(
            candidates=candidates,
            bfile_prefix=bfile_prefix,
            selected_chrom=selected_chrom,
            start_bp=start_bp,
            end_bp=end_bp,
            max_indices=required_n,
            clump_r2=clump_r2,
            plink_path=plink_path,
        )
    else:
        selected = greedy_clump_uploaded_ld_for_auto_indices(
            candidates=candidates,
            max_indices=required_n,
            clump_r2=clump_r2,
            ld_long_table=ld_long_table,
            ld_matrix_table=ld_matrix_table,
        )

    if len(selected) < required_n:
        raise ValueError(
            f"Only {len(selected)} independent index variant(s) were found "
            f"at r² = {clump_r2} within chr{selected_chrom}:"
            f"{start_bp}-{end_bp}. "
            f"The '{st.session_state.get('active_locusblend_mode', 'selected')}' mode "
            f"requires {required_n}. Please enlarge the window, switch to "
            f"a mode requiring fewer index variants, or use manual input."
        )

    n_candidates = len(candidates)
    summary = {
        "n_candidates": n_candidates,
        "n_selected": len(selected),
        "clump_r2": clump_r2,
        "selected_chrom": selected_chrom,
        "start_bp": start_bp,
        "end_bp": end_bp,
    }
    return selected, summary


@st.cache_data(show_spinner=False)
def load_reference_bim(bfile_prefix):
    return read_reference_bim(bfile_prefix)


@st.cache_data(show_spinner=False)
def compute_ld_maps_with_plink(
    bfile_prefix,
    chrom,
    start,
    end,
    window_snps,
    idx1_ref,
    idx2_ref,
    idx3_ref,
    plink_path=str(BIN_DIR / "plink")
):
    return build_ld_maps_with_plink(
        bfile_prefix,
        chrom,
        start,
        end,
        window_snps,
        idx1_ref,
        idx2_ref,
        idx3_ref,
        plink_path,
        load_reference_bim,
    )


@st.cache_data(show_spinner=False)
def load_gene_table(parquet_path):
    return read_gene_table(parquet_path)


def load_genes_from_table(chrom, start, end, parquet_path, gene_display_mode="protein_coding"):
    genes = load_gene_table(parquet_path)
    return filter_genes_from_table(genes, chrom, start, end, gene_display_mode)


@st.cache_data(show_spinner=False)
def load_genes_from_gtf(gtf_path, chrom, start, end, gene_display_mode="protein_coding", chunksize=200000):
    return read_genes_from_gtf(gtf_path, chrom, start, end, gene_display_mode, chunksize)


def get_ld_reference_labels(active_ld_source, ancestry=INTERNAL_1000G_DEFAULT_ANCESTRY):
    """Return tooltip/caption labels keyed to the currently active LD source.

    The downstream builders never hard-code the reference name; they read
    from this dict so the same plot machinery serves the 1000G mode and the
    two uploaded-LD modes.
    """
    if active_ld_source == "Use internal 1000G reference":
        ancestry = normalize_internal_1000g_ancestry(ancestry)
        source_name = f"1000G {ancestry} LD"
        return {
            "source_name": source_name,
            "in_ref": f"reference: in {source_name}",
            "not_in_ref": f"reference: not found in {source_name}",
            "index_not_found": f"Index SNP not found in {source_name}",
            "summary_label": f"1000G {ancestry} window SNPs",
        }
    return {
        "source_name": "user uploaded LD",
        "in_ref": "reference: in user uploaded LD",
        "not_in_ref": "reference: not found in user uploaded LD",
        "index_not_found": "Index SNP not found in user uploaded LD",
        "summary_label": "Uploaded LD SNPs",
    }


def update_progress(progress_bar, status_box, pct, msg):
    if progress_bar is not None:
        progress_bar.progress(pct)
    if status_box is not None:
        status_box.markdown(f"**{msg}**")
    log(msg)


def image_to_data_uri(path):
    mime = guess_type(str(path))[0] or "image/png"
    data = base64.b64encode(path.read_bytes()).decode("utf-8")
    return f"data:{mime};base64,{data}"


def set_streamlit_chrome_minimal():
    """Best-effort app chrome minimization.

    This reduces access to the Streamlit top-right menu in supported
    Streamlit versions. It must never trigger reruns and must never break
    the app if the option is unavailable.
    """
    try:
        st.set_option("client.toolbarMode", "minimal")
    except Exception:
        pass


def inject_locusblend_css():
    st.markdown(
        """
        <style>
        :root {
            color-scheme: light;
            --lb-bg: #ffffff;
            --lb-surface: #ffffff;
            --lb-sidebar-bg: #f8fafc;
            --lb-text: #111827;
            --lb-muted: #4b5563;
            --lb-border: #e5e7eb;
            --lb-input-bg: #ffffff;
            --lb-input-border: #d1d5db;
        }

        html, body, .stApp {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
            color-scheme: light !important;
        }

        [data-testid="stAppViewContainer"],
        [data-testid="stHeader"],
        [data-testid="stToolbar"] {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-testid="stSidebar"],
        [data-testid="stSidebarContent"] {
            background: var(--lb-sidebar-bg) !important;
            color: var(--lb-text) !important;
        }

        .block-container {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
        }

        h1, h2, h3, h4, h5, h6,
        p, li, label,
        [data-testid="stMarkdownContainer"],
        [data-testid="stMarkdownContainer"] p,
        [data-testid="stMarkdownContainer"] li,
        [data-testid="stCaptionContainer"],
        [data-testid="stWidgetLabel"],
        [data-testid="stExpander"],
        [data-testid="stExpander"] details,
        [data-testid="stExpander"] summary {
            color: var(--lb-text) !important;
        }

        [data-testid="stCaptionContainer"],
        .stCaptionContainer {
            color: var(--lb-muted) !important;
        }

        input,
        textarea,
        select,
        [data-baseweb="input"] input,
        [data-baseweb="textarea"] textarea,
        [data-baseweb="select"] div,
        [data-baseweb="select"] span {
            background-color: var(--lb-input-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-baseweb="input"],
        [data-baseweb="textarea"],
        [data-baseweb="select"] {
            background-color: var(--lb-input-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-testid="stDataFrame"],
        [data-testid="stTable"],
        [data-testid="stMetric"],
        [data-testid="stMetricLabel"],
        [data-testid="stMetricValue"] {
            color: var(--lb-text) !important;
        }

        [data-testid="stFileUploader"] {
            color: var(--lb-text) !important;
        }

        [data-testid="stFileUploader"] section {
            background: var(--lb-surface) !important;
            color: var(--lb-text) !important;
            border-color: var(--lb-border) !important;
        }

        #MainMenu {
            visibility: hidden !important;
        }

        [data-testid="stToolbar"] {
            display: none !important;
            visibility: hidden !important;
            height: 0 !important;
        }

        [data-testid="stDecoration"] {
            display: none !important;
            visibility: hidden !important;
        }

        [data-testid="stDeployButton"] {
            display: none !important;
            visibility: hidden !important;
        }

        /* Force Streamlit dialogs and modal content to light mode. */
        [data-testid="stDialog"],
        [data-testid="stDialog"] *,
        div[role="dialog"],
        div[role="dialog"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        [data-testid="stDialog"] button,
        div[role="dialog"] button {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
        }

        /* Force expanders to remain readable in browser/Streamlit dark mode. */
        [data-testid="stExpander"],
        [data-testid="stExpander"] details,
        [data-testid="stExpander"] summary,
        [data-testid="stExpander"] summary *,
        [data-testid="stExpander"] div,
        [data-testid="stExpander"] p,
        [data-testid="stExpander"] label,
        [data-testid="stExpander"] span {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Streamlit/BaseWeb controls, menus and dropdown popovers. */
        [data-baseweb="popover"],
        [data-baseweb="popover"] *,
        [data-baseweb="menu"],
        [data-baseweb="menu"] *,
        [data-baseweb="select"],
        [data-baseweb="select"] *,
        [data-baseweb="input"],
        [data-baseweb="input"] *,
        [data-baseweb="textarea"],
        [data-baseweb="textarea"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Inputs and buttons should remain visible, not hidden. */
        input,
        textarea,
        select,
        button {
            color-scheme: light !important;
        }

        button {
            color: #111827 !important;
        }

        /* Keep Streamlit alerts readable. */
        [data-testid="stAlert"],
        [data-testid="stAlert"] *,
        .stAlert,
        .stAlert * {
            color: #111827 !important;
        }

        /* Keep export/download controls readable. */
        [data-testid="stDownloadButton"],
        [data-testid="stDownloadButton"] *,
        [data-testid="stButton"],
        [data-testid="stButton"] * {
            color: #111827 !important;
        }

        /* --- Force buttons and number steppers to light mode --- */

        /* Normal Streamlit buttons */
        [data-testid="stButton"] button,
        [data-testid="stDownloadButton"] button,
        [data-testid="stFormSubmitButton"] button {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Streamlit primary buttons, including Update plot / Prepare export when primary */
        button[data-testid="baseButton-primary"],
        [data-testid="baseButton-primary"],
        [data-testid="stButton"] button[data-testid="baseButton-primary"],
        [data-testid="stFormSubmitButton"] button[data-testid="baseButton-primary"] {
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Hover states */
        [data-testid="stButton"] button:hover,
        [data-testid="stDownloadButton"] button:hover,
        [data-testid="stFormSubmitButton"] button:hover {
            background-color: #f9fafb !important;
            color: #111827 !important;
            border-color: #9ca3af !important;
        }

        button[data-testid="baseButton-primary"]:hover,
        [data-testid="baseButton-primary"]:hover,
        [data-testid="stButton"] button[data-testid="baseButton-primary"]:hover,
        [data-testid="stFormSubmitButton"] button[data-testid="baseButton-primary"]:hover {
            background-color: #ff3333 !important;
            color: #ffffff !important;
            border-color: #ff3333 !important;
        }

        /* Disabled buttons should be readable, not black */
        [data-testid="stButton"] button:disabled,
        [data-testid="stDownloadButton"] button:disabled,
        [data-testid="stFormSubmitButton"] button:disabled,
        button:disabled,
        button[disabled] {
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* BaseWeb number input stepper buttons, including +/- controls */
        [data-baseweb="input"] button,
        [data-baseweb="input"] [role="button"],
        [data-testid="stNumberInput"] button,
        [data-testid="stNumberInput"] [role="button"] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Icons inside number steppers */
        [data-baseweb="input"] button svg,
        [data-baseweb="input"] [role="button"] svg,
        [data-testid="stNumberInput"] button svg,
        [data-testid="stNumberInput"] [role="button"] svg {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Disabled number stepper buttons */
        [data-baseweb="input"] button:disabled,
        [data-baseweb="input"] [role="button"][aria-disabled="true"],
        [data-testid="stNumberInput"] button:disabled,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] {
            background-color: #f3f4f6 !important;
            color: #9ca3af !important;
            border-color: #d1d5db !important;
            opacity: 1 !important;
        }

        /* Icons inside disabled number steppers */
        [data-baseweb="input"] button:disabled svg,
        [data-baseweb="input"] [role="button"][aria-disabled="true"] svg,
        [data-testid="stNumberInput"] button:disabled svg,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] svg {
            color: #9ca3af !important;
            fill: #9ca3af !important;
            stroke: #9ca3af !important;
        }

        /* Keep selectbox controls light and readable */
        [data-baseweb="select"] > div,
        [data-baseweb="select"] div[role="button"],
        [data-baseweb="select"] svg {
            background-color: #ffffff !important;
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
            color-scheme: light !important;
        }

        /* Checkbox should remain visible on light background */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] span,
        [data-testid="stCheckbox"] div {
            color: #111827 !important;
            color-scheme: light !important;
        }

        .lb-sidebar-visual-card {
            container-type: inline-size;
            background: #ffffff !important;
            border: 1px solid var(--lb-border);
            border-radius: 12px;
            padding: 8px;
            margin: 0.25rem 0 0.75rem 0;
            box-shadow: 0 1px 2px rgba(15, 23, 42, 0.06);
            overflow: hidden;
            width: 100%;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-row {
            display: flex;
            flex-direction: row;
            flex-wrap: nowrap;
            align-items: center;
            justify-content: center;
            gap: 8px;
            width: 100%;
            min-width: 0;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-pane {
            min-width: 0;
            max-width: 100%;
            overflow: hidden;
            display: flex;
            align-items: center;
            justify-content: center;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-pane--legend {
            flex: 1.08 1 0;
        }

        .lb-sidebar-visual-pane--diagram {
            flex: 0.92 1 0;
        }

        .lb-sidebar-visual-pane--single {
            flex: 1 1 auto;
        }

        .lb-sidebar-visual-img {
            display: block;
            width: 100%;
            max-width: 100%;
            min-width: 0;
            height: auto;
            max-height: 125px;
            object-fit: contain;
            box-sizing: border-box;
        }

        @container (max-width: 280px) {
            .lb-sidebar-visual-card {
                padding: 6px;
            }

            .lb-sidebar-visual-row {
                gap: 4px;
            }

            .lb-sidebar-visual-img {
                max-height: 115px;
            }
        }

        @container (max-width: 220px) {
            .lb-sidebar-visual-card {
                padding: 4px;
            }

            .lb-sidebar-visual-row {
                gap: 3px;
            }

            .lb-sidebar-visual-img {
                max-height: 100px;
            }
        }

        /* --- Exact Streamlit 1.50+ stBaseButton override --- */

        /* Exact Streamlit button selectors observed in Chrome DevTools.
           Do not target st-emotion-cache-* classes because they are generated. */

        button[data-testid="stBaseButton-secondary"],
        button[data-testid="stBaseButton-tertiary"],
        button[kind="secondary"],
        button[kind="tertiary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Text, spans, and icons inside secondary/tertiary buttons */
        button[data-testid="stBaseButton-secondary"] *,
        button[data-testid="stBaseButton-tertiary"] *,
        button[kind="secondary"] *,
        button[kind="tertiary"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Primary Streamlit buttons */
        button[data-testid="stBaseButton-primary"],
        button[kind="primary"] {
            background: #ff4b4b !important;
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Text, spans, and icons inside primary buttons */
        button[data-testid="stBaseButton-primary"] *,
        button[kind="primary"] * {
            color: #ffffff !important;
            fill: #ffffff !important;
            stroke: #ffffff !important;
        }

        /* Generic stBaseButton fallback */
        button[data-testid^="stBaseButton"] {
            color-scheme: light !important;
            box-shadow: none !important;
        }

        /* Disabled buttons: make them light gray, not black */
        button[data-testid^="stBaseButton"]:disabled,
        button[data-testid^="stBaseButton"][disabled],
        button[data-testid^="stBaseButton"][aria-disabled="true"],
        button[kind]:disabled,
        button[kind][disabled],
        button[kind][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            cursor: not-allowed !important;
            color-scheme: light !important;
        }

        /* Disabled button inner text/icons */
        button[data-testid^="stBaseButton"]:disabled *,
        button[data-testid^="stBaseButton"][disabled] *,
        button[data-testid^="stBaseButton"][aria-disabled="true"] *,
        button[kind]:disabled *,
        button[kind][disabled] *,
        button[kind][aria-disabled="true"] * {
            color: #6b7280 !important;
            fill: #6b7280 !important;
            stroke: #6b7280 !important;
        }

        /* Hover states for enabled secondary/tertiary buttons */
        button[data-testid="stBaseButton-secondary"]:not(:disabled):hover,
        button[data-testid="stBaseButton-tertiary"]:not(:disabled):hover,
        button[kind="secondary"]:not(:disabled):hover,
        button[kind="tertiary"]:not(:disabled):hover {
            background: #f9fafb !important;
            background-color: #f9fafb !important;
            color: #111827 !important;
            border-color: #9ca3af !important;
        }

        /* Hover states for enabled primary buttons */
        button[data-testid="stBaseButton-primary"]:not(:disabled):hover,
        button[kind="primary"]:not(:disabled):hover {
            background: #ff3333 !important;
            background-color: #ff3333 !important;
            color: #ffffff !important;
            border-color: #ff3333 !important;
        }

        /* File uploader buttons and inner labels */
        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] button[kind],
        [data-testid="stFileUploader"] button {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"] *,
        [data-testid="stFileUploader"] button[kind] *,
        [data-testid="stFileUploader"] button * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* File uploader dropzone and text */
        [data-testid="stFileUploader"] section,
        [data-testid="stFileUploaderDropzone"],
        [data-testid="stFileUploadDropzone"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] section *,
        [data-testid="stFileUploaderDropzone"] *,
        [data-testid="stFileUploadDropzone"] * {
            color: #111827 !important;
        }

        /* Number input text field */
        [data-testid="stNumberInput"] input,
        [data-baseweb="input"] input {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        /* Number input stepper controls.
           Streamlit/BaseWeb can use buttons, role=button, or aria-label wrappers. */
        [data-testid="stNumberInput"] button,
        [data-testid="stNumberInput"] button[data-testid^="stBaseButton"],
        [data-testid="stNumberInput"] button[kind],
        [data-testid="stNumberInput"] [role="button"],
        [data-testid="stNumberInput"] div[aria-label],
        [data-testid="stNumberInput"] span[aria-label],
        [data-baseweb="input"] button,
        [data-baseweb="input"] button[data-testid^="stBaseButton"],
        [data-baseweb="input"] button[kind],
        [data-baseweb="input"] [role="button"],
        [data-baseweb="input"] div[aria-label],
        [data-baseweb="input"] span[aria-label] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Number input stepper icons */
        [data-testid="stNumberInput"] button *,
        [data-testid="stNumberInput"] button[data-testid^="stBaseButton"] *,
        [data-testid="stNumberInput"] button[kind] *,
        [data-testid="stNumberInput"] [role="button"] *,
        [data-testid="stNumberInput"] div[aria-label] *,
        [data-testid="stNumberInput"] span[aria-label] *,
        [data-baseweb="input"] button *,
        [data-baseweb="input"] button[data-testid^="stBaseButton"] *,
        [data-baseweb="input"] button[kind] *,
        [data-baseweb="input"] [role="button"] *,
        [data-baseweb="input"] div[aria-label] *,
        [data-baseweb="input"] span[aria-label] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Disabled number input steppers */
        [data-testid="stNumberInput"] button:disabled,
        [data-testid="stNumberInput"] button[disabled],
        [data-testid="stNumberInput"] button[aria-disabled="true"],
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"],
        [data-baseweb="input"] button:disabled,
        [data-baseweb="input"] button[disabled],
        [data-baseweb="input"] button[aria-disabled="true"],
        [data-baseweb="input"] [role="button"][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #9ca3af !important;
            border-color: #d1d5db !important;
            opacity: 1 !important;
        }

        /* Disabled number input stepper icons */
        [data-testid="stNumberInput"] button:disabled *,
        [data-testid="stNumberInput"] button[disabled] *,
        [data-testid="stNumberInput"] button[aria-disabled="true"] *,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] *,
        [data-baseweb="input"] button:disabled *,
        [data-baseweb="input"] button[disabled] *,
        [data-baseweb="input"] button[aria-disabled="true"] *,
        [data-baseweb="input"] [role="button"][aria-disabled="true"] * {
            color: #9ca3af !important;
            fill: #9ca3af !important;
            stroke: #9ca3af !important;
        }

        /* Download button */
        [data-testid="stDownloadButton"] button[data-testid^="stBaseButton"],
        [data-testid="stDownloadButton"] button[kind],
        [data-testid="stDownloadButton"] button {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stDownloadButton"] button[data-testid^="stBaseButton"] *,
        [data-testid="stDownloadButton"] button[kind] *,
        [data-testid="stDownloadButton"] button * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Radio buttons and checkbox text should remain readable */
        [data-testid="stRadio"] *,
        [data-testid="stCheckbox"] * {
            color-scheme: light !important;
        }

        [data-testid="stRadio"] label,
        [data-testid="stRadio"] span,
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] span {
            color: #111827 !important;
        }

        /* --- Export button and checkbox repair --- */

        /* Repair normal/secondary Streamlit buttons, including Prepare export file. */
        button[data-testid="stBaseButton-secondary"],
        button[kind="secondary"],
        [data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stButton"] button[kind="secondary"],
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stFormSubmitButton"] button[kind="secondary"],
        [data-testid="stDownloadButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stDownloadButton"] button[kind="secondary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Repair text inside secondary buttons.
           Do not force SVG fill/stroke here; broad SVG rules can break checkbox marks. */
        button[data-testid="stBaseButton-secondary"] span,
        button[kind="secondary"] span,
        [data-testid="stButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stButton"] button[kind="secondary"] span,
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stFormSubmitButton"] button[kind="secondary"] span,
        [data-testid="stDownloadButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stDownloadButton"] button[kind="secondary"] span {
            color: #111827 !important;
        }

        /* Keep primary buttons readable if any remain. */
        button[data-testid="stBaseButton-primary"],
        button[kind="primary"],
        [data-testid="stButton"] button[data-testid="stBaseButton-primary"],
        [data-testid="stButton"] button[kind="primary"],
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-primary"],
        [data-testid="stFormSubmitButton"] button[kind="primary"] {
            background: #ff4b4b !important;
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        button[data-testid="stBaseButton-primary"] span,
        button[kind="primary"] span,
        [data-testid="stButton"] button[data-testid="stBaseButton-primary"] span,
        [data-testid="stButton"] button[kind="primary"] span,
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-primary"] span,
        [data-testid="stFormSubmitButton"] button[kind="primary"] span {
            color: #ffffff !important;
        }

        /* Disabled buttons should be light gray and readable. */
        button[data-testid^="stBaseButton"]:disabled,
        button[data-testid^="stBaseButton"][disabled],
        button[data-testid^="stBaseButton"][aria-disabled="true"],
        button[kind]:disabled,
        button[kind][disabled],
        button[kind][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        button[data-testid^="stBaseButton"]:disabled span,
        button[data-testid^="stBaseButton"][disabled] span,
        button[data-testid^="stBaseButton"][aria-disabled="true"] span,
        button[kind]:disabled span,
        button[kind][disabled] span,
        button[kind][aria-disabled="true"] span {
            color: #6b7280 !important;
        }

        /* Checkbox repair:
           Only control the label text. Do NOT override checkbox internal spans/SVGs/marks,
           because Streamlit/BaseWeb uses those to draw the checkmark. */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] label p,
        [data-testid="stCheckbox"] label span:not([data-baseweb]) {
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Restore checkbox box visibility without forcing the inner checkmark SVG. */
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] {
            color-scheme: light !important;
        }

        /* Checkbox input should remain selectable and visible. */
        [data-testid="stCheckbox"] input[type="checkbox"] {
            accent-color: #ff4b4b !important;
        }

        /* Remove overly broad checkbox icon overrides from winning.
           This intentionally avoids fill/stroke rules on checkbox descendants. */
        [data-testid="stCheckbox"] svg {
            color: revert !important;
            fill: revert !important;
            stroke: revert !important;
        }

        /* File uploader Browse files button remains secondary-style. */
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stFileUploader"] button[kind="secondary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
        }

        [data-testid="stFileUploader"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stFileUploader"] button[kind="secondary"] span {
            color: #111827 !important;
        }

        /* --- Markdown code and documentation readability --- */

        /* Inline markdown code: prevent black background / green text in browser dark mode. */
        [data-testid="stMarkdownContainer"] code,
        [data-testid="stMarkdownContainer"] p code,
        [data-testid="stMarkdownContainer"] li code,
        [data-testid="stMarkdownContainer"] td code,
        [data-testid="stMarkdownContainer"] th code,
        [data-testid="stExpander"] code,
        [data-testid="stExpander"] p code,
        [data-testid="stExpander"] li code,
        code {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #111827 !important;
            border: 1px solid #e5e7eb !important;
            border-radius: 4px !important;
            padding: 0.08rem 0.28rem !important;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace !important;
            font-size: 0.88em !important;
            white-space: break-spaces !important;
            color-scheme: light !important;
        }

        /* Code blocks should also stay light. */
        [data-testid="stMarkdownContainer"] pre,
        [data-testid="stMarkdownContainer"] pre code,
        [data-testid="stExpander"] pre,
        [data-testid="stExpander"] pre code,
        pre,
        pre code {
            background: #f8fafc !important;
            background-color: #f8fafc !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        /* Markdown tables in the User Guide should remain readable. */
        [data-testid="stMarkdownContainer"] table,
        [data-testid="stMarkdownContainer"] thead,
        [data-testid="stMarkdownContainer"] tbody,
        [data-testid="stMarkdownContainer"] tr,
        [data-testid="stMarkdownContainer"] th,
        [data-testid="stMarkdownContainer"] td {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        /* --- Compact sidebar sections --- */

        /* Reduce sidebar inner padding and vertical gaps. */
        [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
            padding-top: 1.0rem !important;
            padding-bottom: 1.0rem !important;
        }

        /* Keep sidebar markdown/caption spacing tighter. */
        [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] {
            margin-bottom: 0.25rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stCaptionContainer"],
        [data-testid="stSidebar"] .stCaptionContainer {
            margin-top: 0.15rem !important;
            margin-bottom: 0.45rem !important;
            line-height: 1.35 !important;
        }

        /* Compact sidebar expanders while keeping them easy to select. */
        [data-testid="stSidebar"] [data-testid="stExpander"] {
            margin-top: 0.25rem !important;
            margin-bottom: 0.45rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stExpander"] details {
            border-radius: 8px !important;
        }

        [data-testid="stSidebar"] [data-testid="stExpander"] summary {
            min-height: 2.35rem !important;
            padding-top: 0.45rem !important;
            padding-bottom: 0.45rem !important;
        }

        /* Reduce oversized blank gaps between Streamlit element containers in sidebar. */
        [data-testid="stSidebar"] [data-testid="stVerticalBlock"] {
            gap: 0.35rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stElementContainer"] {
            margin-bottom: 0.20rem !important;
        }

        /* Make section helper text compact. */
        [data-testid="stSidebar"] p {
            line-height: 1.35 !important;
        }

        /* Keep form submit button close to final section. */
        [data-testid="stSidebar"] [data-testid="stFormSubmitButton"] {
            margin-top: 0.35rem !important;
            margin-bottom: 0.35rem !important;
        }

        /* --- Checkbox and uploader control visibility --- */

        /* Make Streamlit/BaseWeb checkbox boxes and ticks visible in forced light mode. */
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] > div,
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] div[role="checkbox"],
        [data-testid="stCheckbox"] div[role="checkbox"] {
            background-color: #ffffff !important;
            border-color: #ff4b4b !important;
            color-scheme: light !important;
        }

        /* Checked checkbox state. */
        [data-testid="stCheckbox"] input[type="checkbox"]:checked + div,
        [data-testid="stCheckbox"] [aria-checked="true"],
        [data-testid="stCheckbox"] div[role="checkbox"][aria-checked="true"] {
            background-color: #ff4b4b !important;
            border-color: #ff4b4b !important;
        }

        /* Checkbox tick/check icon. Avoid broad rules on all checkbox descendants. */
        [data-testid="stCheckbox"] [aria-checked="true"] svg,
        [data-testid="stCheckbox"] div[role="checkbox"][aria-checked="true"] svg {
            color: #ffffff !important;
            fill: #ffffff !important;
            stroke: #ffffff !important;
        }

        /* Checkbox label remains readable. */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] label p,
        [data-testid="stCheckbox"] label span {
            color: #111827 !important;
        }

        /* Browser-native checkbox fallback. */
        [data-testid="stCheckbox"] input[type="checkbox"] {
            accent-color: #ff4b4b !important;
        }

        /* Uploaded file row and remove controls should not appear black in browser dark mode. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Uploaded-file remove button / small icon buttons. */
        [data-testid="stFileUploader"] button,
        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] button[kind],
        [data-testid="stFileUploader"] [role="button"] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] button svg,
        [data-testid="stFileUploader"] button *,
        [data-testid="stFileUploader"] [role="button"] svg,
        [data-testid="stFileUploader"] [role="button"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* If Streamlit renders a small dark delete/remove pill, force it light. */
        [data-testid="stFileUploader"] [aria-label*="remove" i],
        [data-testid="stFileUploader"] [aria-label*="delete" i],
        [data-testid="stFileUploader"] [title*="remove" i],
        [data-testid="stFileUploader"] [title*="delete" i] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] [aria-label*="remove" i] *,
        [data-testid="stFileUploader"] [aria-label*="delete" i] *,
        [data-testid="stFileUploader"] [title*="remove" i] *,
        [data-testid="stFileUploader"] [title*="delete" i] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* --- Uploaded-file remove-control fix --- */

        /* Uploaded file row containers. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] div {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Uploaded file row text and file-size text. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] span,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] p,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] small,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] span,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] p,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] small {
            color: #111827 !important;
        }

        /* Uploaded-file row buttons, including remove/delete controls. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button[kind],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button[kind],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            border-radius: 8px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Uploaded-file row button icons. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button *,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"] svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button svg,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] svg,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Directly target common uploader remove/delete button labels. */
        [data-testid="stFileUploader"] button[aria-label*="remove" i],
        [data-testid="stFileUploader"] button[aria-label*="delete" i],
        [data-testid="stFileUploader"] button[aria-label*="clear" i],
        [data-testid="stFileUploader"] button[title*="remove" i],
        [data-testid="stFileUploader"] button[title*="delete" i],
        [data-testid="stFileUploader"] button[title*="clear" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="remove" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="delete" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="clear" i],
        [data-testid="stFileUploader"] [role="button"][title*="remove" i],
        [data-testid="stFileUploader"] [role="button"][title*="delete" i],
        [data-testid="stFileUploader"] [role="button"][title*="clear" i] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            border-radius: 8px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Icons inside direct remove/delete targets. */
        [data-testid="stFileUploader"] button[aria-label*="remove" i] *,
        [data-testid="stFileUploader"] button[aria-label*="delete" i] *,
        [data-testid="stFileUploader"] button[aria-label*="clear" i] *,
        [data-testid="stFileUploader"] button[title*="remove" i] *,
        [data-testid="stFileUploader"] button[title*="delete" i] *,
        [data-testid="stFileUploader"] button[title*="clear" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="remove" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="delete" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="clear" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="remove" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="delete" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="clear" i] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Some Streamlit versions render the file-row action as the last child.
        Keep this scoped to the uploaded file row only. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div:last-child,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div:last-child *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div:last-child,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div:last-child * {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        /* Do not let the file row action inherit dark theme surfaces. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [data-baseweb],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [data-baseweb] {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* If this still fails, inspect the dark control in DevTools and add its stable data-testid or aria-label selector. Avoid st-emotion-cache-* classes. */

        /* --- Uploaded-file internal scrollbar fix --- */

        /* The remaining black vertical pill in uploaded-file rows is likely an
        internal scrollbar thumb, not a button. Keep this scoped to file uploader. */
        [data-testid="stFileUploader"],
        [data-testid="stFileUploader"] *,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] {
            scrollbar-color: #cbd5e1 #ffffff !important;
            scrollbar-width: thin !important;
            color-scheme: light !important;
        }

        /* WebKit / Chrome scrollbar track inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] *::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar {
            width: 8px !important;
            height: 8px !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
        }

        /* WebKit / Chrome scrollbar thumb inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] *::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar-thumb {
            background: #cbd5e1 !important;
            background-color: #cbd5e1 !important;
            border: 2px solid #ffffff !important;
            border-radius: 999px !important;
        }

        /* WebKit / Chrome scrollbar corner inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] *::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar-corner {
            background: #ffffff !important;
            background-color: #ffffff !important;
        }

        /* Keep the uploaded file row surface light even when the scrollbar is present. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Do not make the file row action black when the browser creates overlay scrollbars. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] * {
            color-scheme: light !important;
        }

        /* Optional: make uploader file rows less likely to create tiny internal scrollbars. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] {
            overflow: visible !important;
        }

        /* Fallback for Streamlit file-row wrappers that create a tiny scrollable box. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div {
            scrollbar-color: #cbd5e1 #ffffff !important;
            scrollbar-width: thin !important;
            color-scheme: light !important;
        }

        /* --- Uploaded-file delete button fix --- */

        /* DevTools-confirmed target:
        div[data-testid="stFileUploaderDeleteBtn"]
        > button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] {
            background: transparent !important;
            background-color: transparent !important;
            border: none !important;
            box-shadow: none !important;
            color-scheme: light !important;
            display: flex !important;
            align-items: center !important;
            justify-content: center !important;
        }

        /* Exact remove button. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button[data-testid="stBaseButton-minimal"],
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "],
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] {
            width: 24px !important;
            height: 24px !important;
            min-width: 24px !important;
            min-height: 24px !important;
            max-width: 24px !important;
            max-height: 24px !important;
            padding: 0 !important;
            margin: 0 !important;
            display: inline-flex !important;
            align-items: center !important;
            justify-content: center !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #64748b !important;
            border: 1px solid #cbd5e1 !important;
            border-radius: 999px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Hover/focus state: keep light, not black. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:active,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:focus,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:active {
            background: #f8fafc !important;
            background-color: #f8fafc !important;
            color: #334155 !important;
            border: 1px solid #94a3b8 !important;
            box-shadow: none !important;
            outline: none !important;
        }

        /* Remove icon SVG sizing and color. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] svg,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] svg,
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] svg {
            width: 14px !important;
            height: 14px !important;
            color: #64748b !important;
            fill: none !important;
            stroke: #64748b !important;
            background: transparent !important;
            background-color: transparent !important;
        }

        /* Remove icon path. Critical: do not give the path a white/black background box. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] svg path,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] svg path,
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] svg path {
            color: #64748b !important;
            fill: none !important;
            stroke: #64748b !important;
            background: transparent !important;
            background-color: transparent !important;
        }

        /* Hover/focus icon color. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover svg path,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus svg path,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover svg,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover svg path {
            color: #334155 !important;
            fill: none !important;
            stroke: #334155 !important;
        }

        /* Keep uploaded file row layout stable. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] {
            align-items: center !important;
        }

        /* Avoid older broad rules making the delete icon into a dark block. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] *,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] * {
            box-shadow: none !important;
            text-shadow: none !important;
        }

        /* --- WashU Medicine header --- */

        .lb-washu-header {
            width: 100%;
            background: #A51417 !important;
            background-color: #A51417 !important;
            color: #ffffff !important;
            border-radius: 0;
            margin: -0.5rem 0 1.25rem 0;
            padding: 0;
            box-shadow: none;
            color-scheme: light !important;
        }

        .lb-washu-header-inner {
            min-height: 44px;
            display: flex;
            align-items: center;
            justify-content: flex-start;
            padding: 0.45rem 1.25rem;
        }

        .lb-washu-logo {
            display: block;
            height: 30px;
            max-width: 260px;
            object-fit: contain;
        }

        .lb-washu-logo-fallback {
            font-size: 1.25rem;
            font-weight: 700;
            letter-spacing: 0;
            color: #ffffff !important;
        }

        @media (max-width: 700px) {
            .lb-washu-header {
                margin-top: -0.25rem;
                margin-bottom: 1rem;
            }

            .lb-washu-header-inner {
                min-height: 40px;
                padding: 0.4rem 0.85rem;
            }

            .lb-washu-logo {
                height: 25px;
                max-width: 220px;
            }
        }

        /* --- Full-width fixed WashU Medicine header --- */

        :root {
            --lb-washu-header-height: 54px;
            --lb-washu-red: #A51417;
        }

        /* Fixed global header across sidebar + main content. */
        .lb-washu-global-header {
            position: fixed !important;
            top: 0 !important;
            left: 0 !important;
            right: 0 !important;
            width: 100vw !important;
            height: var(--lb-washu-header-height) !important;
            z-index: 999990 !important;
            background: var(--lb-washu-red) !important;
            background-color: var(--lb-washu-red) !important;
            color: #ffffff !important;
            display: flex !important;
            align-items: center !important;
            justify-content: flex-start !important;
            margin: 0 !important;
            padding: 0 !important;
            border: none !important;
            border-radius: 0 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Header content width. Keep logo aligned with main app content, but header background spans all. */
        .lb-washu-global-inner {
            width: 100% !important;
            height: var(--lb-washu-header-height) !important;
            display: flex !important;
            align-items: center !important;
            justify-content: flex-start !important;
            padding: 0 1.5rem !important;
            box-sizing: border-box !important;
        }

        /* Logo in fixed header. */
        .lb-washu-header-link {
            display: inline-flex !important;
            align-items: center !important;
            text-decoration: none !important;
            color: inherit !important;
        }

        .lb-washu-header-link:visited,
        .lb-washu-header-link:hover,
        .lb-washu-header-link:active {
            text-decoration: none !important;
            color: inherit !important;
        }

        .lb-washu-header-link img {
            display: block !important;
        }

        .lb-washu-global-header .lb-washu-logo {
            display: block !important;
            height: 32px !important;
            max-width: 280px !important;
            object-fit: contain !important;
        }

        /* Fallback text if local logo is unavailable. */
        .lb-washu-global-header .lb-washu-logo-fallback {
            font-size: 1.25rem !important;
            font-weight: 700 !important;
            color: #ffffff !important;
            letter-spacing: 0.01em !important;
        }

        /* Hide/neutralize the legacy non-fixed header container if still rendered. */
        .lb-washu-header {
            display: none !important;
        }

        /* Push the Streamlit app content below the fixed header. */
        [data-testid="stAppViewContainer"] {
            padding-top: var(--lb-washu-header-height) !important;
        }

        /* Push sidebar content below the fixed header. */
        [data-testid="stSidebar"] {
            padding-top: var(--lb-washu-header-height) !important;
        }

        /* Ensure sidebar background begins below header, while header still spans above it. */
        [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
            padding-top: 1rem !important;
        }

        /* Main block should not add another huge top gap. */
        [data-testid="stMain"] .block-container {
            padding-top: 2rem !important;
        }

        /* Streamlit's own top header can otherwise create a blank strip.
           Keep it transparent and visually minimized without removing app controls. */
        [data-testid="stHeader"] {
            background: transparent !important;
            height: 0 !important;
            min-height: 0 !important;
        }

        /* Keep Streamlit sidebar controls selectable above the banner if present. */
        [data-testid="stSidebarCollapseButton"] {
            z-index: 1000000 !important;
        }

        /* Responsive header. */
        @media (max-width: 700px) {
            :root {
                --lb-washu-header-height: 48px;
            }

            .lb-washu-global-inner {
                padding: 0 1rem !important;
            }

            .lb-washu-global-header .lb-washu-logo {
                height: 27px !important;
                max-width: 230px !important;
            }

            [data-testid="stMain"] .block-container {
                padding-top: 1.5rem !important;
            }
        }

        /* --- Keep sidebar expanded for review/demo --- */

        /* Do not attempt to restyle Streamlit's collapsed/reopen button.
        For this review build, prevent users from collapsing the sidebar.
        This avoids the known issue where the sidebar can be hard to reopen. */

        /* Hide only the native "close/collapse sidebar" control when the sidebar is expanded.
        Keep the selector narrow and scoped to the sidebar header. */
        [data-testid="stSidebar"] [data-testid="stSidebarCollapseButton"],
        [data-testid="stSidebar"] [data-testid="stSidebarHeader"] [data-testid="stSidebarCollapseButton"] {
            display: none !important;
            visibility: hidden !important;
            pointer-events: none !important;
        }

        /* Keep the sidebar itself visible and normal. */
        [data-testid="stSidebar"] {
            visibility: visible !important;
            opacity: 1 !important;
        }

        /* Preserve the app-wide WashU header.
        Do not touch native sidebar reopen controls here. */

        /* --- Remove empty sidebar header gap --- */

        /* In this review/demo build, the sidebar collapse button is intentionally hidden.
        Streamlit's sidebar header container then becomes empty but still occupies
        space above the Legend. Collapse only that empty header area. */
        [data-testid="stSidebar"] [data-testid="stSidebarHeader"] {
            height: 0 !important;
            min-height: 0 !important;
            max-height: 0 !important;
            padding: 0 !important;
            margin: 0 !important;
            overflow: hidden !important;
        }

        /* Keep the actual collapse button hidden. */

        /* Do not change the global WashU header or Streamlit main/header layers here. */

        /* --- Sidebar app title --- */

        /* Keep the native Streamlit sidebar collapse mechanics untouched here.
           This build adds an app-level sidebar title instead of trying to restyle
           native sidebar controls. */

        /* Use the sidebar top area for an intentional app title. */
        .lb-sidebar-app-header {
            margin: 0 0 0.85rem 0 !important;
            padding: 0.85rem 0.9rem !important;
            border-radius: 10px !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
            border: 1px solid #e5e7eb !important;
            color: #111827 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        .lb-sidebar-app-title {
            font-size: 1.05rem !important;
            font-weight: 800 !important;
            line-height: 1.2 !important;
            color: #111827 !important;
            margin: 0 !important;
            padding: 0 !important;
        }

        .lb-sidebar-app-subtitle {
            font-size: 0.78rem !important;
            font-weight: 500 !important;
            line-height: 1.25 !important;
            color: #6b7280 !important;
            margin-top: 0.25rem !important;
            padding: 0 !important;
        }

        /* Reduce only the visible top spacing inside sidebar user content.
           Do not touch native sidebar buttons. */
        [data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
            padding-top: 0.75rem !important;
        }

        /* Keep the empty Streamlit sidebar header compact. */

        /* Keep the sidebar collapse button hidden for this review/demo build. */

        /* --- Cleaner sidebar top label --- */

        /* Use a compact label instead of a heavy sidebar title card. */
        .lb-sidebar-app-header {
            display: none !important;
        }

        /* Compact sidebar label with a WashU-red accent. */
        .lb-sidebar-app-label {
            margin: 0.25rem 0 1.1rem 0 !important;
            padding: 0.15rem 0 0.15rem 0.7rem !important;
            border-left: 4px solid #A51417 !important;
            background: transparent !important;
            color: #111827 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        .lb-sidebar-app-label span {
            font-size: 0.95rem !important;
            font-weight: 800 !important;
            line-height: 1.2 !important;
            color: #111827 !important;
            letter-spacing: 0.01em !important;
        }

        /* Reduce only the visible top spacing inside sidebar user content.
           Do not touch native sidebar buttons. */
        [data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
            padding-top: 0.65rem !important;
        }

        /* Keep the empty Streamlit sidebar header compact. */

        /* Keep the sidebar collapse button hidden for this review/demo build. */

        /* --- Top update button hint --- */

        .lb-sidebar-update-hint {
            margin: 0.25rem 0 0.9rem 0 !important;
            padding: 0 !important;
            font-size: 0.76rem !important;
            line-height: 1.25 !important;
            color: #6b7280 !important;
        }

        /* --- Top-level Visualizer / Documentation switcher --- */

        [data-testid="stRadio"] {
            color-scheme: light !important;
        }

        /* Keep the page switcher visually compact. */
        .lb-doc-note {
            color: #6b7280 !important;
            font-size: 0.9rem !important;
            line-height: 1.45 !important;
        }

        /* Documentation page readability. */
        .lb-doc-section {
            margin-top: 1.25rem !important;
            margin-bottom: 1.25rem !important;
        }

        .lb-doc-section h2,
        .lb-doc-section h3 {
            color: #111827 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_sidebar_app_header():
    """Render a compact app-level label at the top of the sidebar."""
    st.sidebar.markdown(
        """
<div class="lb-sidebar-app-label">
  <span>LocusBlend Controls</span>
</div>
        """,
        unsafe_allow_html=True,
    )


def render_sidebar_visual_strip():
    legend_overlay_path = ASSET_DIR / "legend_overlay.png"
    drawing3_path = ASSET_DIR / "Drawing3.png"

    items = []
    if legend_overlay_path.exists():
        items.append(("legend", "LD legend", legend_overlay_path))
    if drawing3_path.exists():
        items.append(("diagram", "LocusBlend diagram", drawing3_path))

    if not items:
        return

    single = len(items) == 1
    html_chunks = []
    for kind, alt, path in items:
        src = image_to_data_uri(path)
        safe_alt = html_lib.escape(alt, quote=True)
        pane_class = "lb-sidebar-visual-pane--single" if single else f"lb-sidebar-visual-pane--{kind}"
        html_chunks.append(
            f'''
            <div class="lb-sidebar-visual-pane {pane_class}">
                <img src="{src}" alt="{safe_alt}" class="lb-sidebar-visual-img" />
            </div>
            '''
        )

    visual_html = f"""
    <!doctype html>
    <html>
    <head>
    <style>
    html, body {{
        margin: 0;
        padding: 0;
        background: transparent;
        overflow: hidden;
        color-scheme: light;
        font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}

    .lb-sidebar-visual-card {{
        width: 100%;
        box-sizing: border-box;
        background: #ffffff;
        border: 1px solid #e5e7eb;
        border-radius: 12px;
        padding: 8px;
        box-shadow: 0 1px 2px rgba(15, 23, 42, 0.06);
        overflow: hidden;
    }}

    .lb-sidebar-visual-row {{
        display: flex;
        flex-direction: row;
        flex-wrap: nowrap;
        align-items: center;
        justify-content: center;
        gap: 8px;
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
    }}

    .lb-sidebar-visual-pane {{
        min-width: 0;
        max-width: 100%;
        overflow: hidden;
        display: flex;
        align-items: center;
        justify-content: center;
        box-sizing: border-box;
    }}

    .lb-sidebar-visual-pane--legend {{
        flex: 1.08 1 0;
    }}

    .lb-sidebar-visual-pane--diagram {{
        flex: 0.92 1 0;
    }}

    .lb-sidebar-visual-pane--single {{
        flex: 1 1 auto;
    }}

    .lb-sidebar-visual-img {{
        display: block;
        width: 100%;
        max-width: 100%;
        min-width: 0;
        height: auto;
        max-height: 125px;
        object-fit: contain;
        box-sizing: border-box;
    }}

    @media (max-width: 280px) {{
        .lb-sidebar-visual-card {{
            padding: 6px;
        }}

        .lb-sidebar-visual-row {{
            gap: 4px;
        }}

        .lb-sidebar-visual-img {{
            max-height: 115px;
        }}
    }}

    @media (max-width: 220px) {{
        .lb-sidebar-visual-card {{
            padding: 4px;
        }}

        .lb-sidebar-visual-row {{
            gap: 3px;
        }}

        .lb-sidebar-visual-img {{
            max-height: 100px;
        }}
    }}
    </style>
    </head>
    <body>
        <div class="lb-sidebar-visual-card">
            <div class="lb-sidebar-visual-row">
                {''.join(html_chunks)}
            </div>
        </div>
    </body>
    </html>
    """

    st.sidebar.markdown("### LD legend")
    with st.sidebar:
        components.html(visual_html, height=165, scrolling=False)


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


def make_locusblend_export_bytes(
    output_format,
    width_in,
    height_in,
    dpi,
    include_legends,
):
    locus_fig = st.session_state.get("last_locus_fig")
    compare_fig = st.session_state.get("last_compare_fig")
    export_context = st.session_state.get("last_export_context", {})
    compare_mode = st.session_state.get("last_compare_mode", export_context.get("compare_mode", ""))

    if locus_fig is None or compare_fig is None:
        raise ValueError("No cached plot is available yet. Select Update plot first.")

    export_context = dict(export_context or {})
    export_context["compare_mode"] = compare_mode

    canvas = build_locusblend_export_image(
        locus_fig=locus_fig,
        compare_fig=compare_fig,
        export_context=export_context,
        output_width_in=float(width_in),
        output_height_in=float(height_in),
        dpi=int(dpi),
        include_legends=bool(include_legends),
    )

    buf = BytesIO()
    fmt = str(output_format).upper()

    if fmt == "PNG":
        canvas.save(buf, format="PNG")
        return buf.getvalue(), "image/png", "locusblend_export.png"

    if fmt == "PDF":
        canvas.convert("RGB").save(buf, format="PDF", resolution=int(dpi))
        return buf.getvalue(), "application/pdf", "locusblend_export.pdf"

    raise ValueError(f"Unsupported export format: {output_format}")


def render_export_controls():
    if "last_locus_fig" not in st.session_state or "last_compare_fig" not in st.session_state:
        return

    st.markdown("### Export / Print current figure")
    with st.expander("Export settings", expanded=False):
        st.caption(
            "Exports the currently rendered cached plot. "
            "Sidebar edits are not applied until you select Update plot."
        )

        output_format = st.selectbox(
            "Output format",
            ["PNG", "PDF"],
            index=0,
            key="export_output_format",
        )

        paper_preset = st.selectbox(
            "Page size",
            [
                "US Letter portrait (8.5 x 11 in)",
                "US Letter landscape (11 x 8.5 in)",
                "Custom",
            ],
            index=0,
            key="export_paper_preset",
        )

        if paper_preset == "US Letter portrait (8.5 x 11 in)":
            default_w, default_h = 8.5, 11.0
        elif paper_preset == "US Letter landscape (11 x 8.5 in)":
            default_w, default_h = 11.0, 8.5
        else:
            default_w = float(st.session_state.get("export_custom_width_in", 8.5))
            default_h = float(st.session_state.get("export_custom_height_in", 11.0))

        if paper_preset == "Custom":
            col_w, col_h = st.columns(2)
            with col_w:
                width_in = st.number_input(
                    "Width (in)",
                    min_value=4.0,
                    max_value=30.0,
                    value=float(default_w),
                    step=0.25,
                    key="export_custom_width_in",
                )
            with col_h:
                height_in = st.number_input(
                    "Height (in)",
                    min_value=4.0,
                    max_value=30.0,
                    value=float(default_h),
                    step=0.25,
                    key="export_custom_height_in",
                )
        else:
            width_in, height_in = default_w, default_h
            st.caption(f"Output size: {width_in:g} x {height_in:g} inches")

        col_dpi, col_leg = st.columns(2)
        with col_dpi:
            dpi = st.number_input(
                "DPI",
                min_value=72,
                max_value=600,
                value=300,
                step=25,
                key="export_dpi",
            )
        with col_leg:
            _export_legend_default = bool(st.session_state.get("export_include_legends", True))
            _export_legend_choice = st.selectbox(
                "Export legend images",
                options=["Include", "Do not include"],
                index=0 if _export_legend_default else 1,
                key="export_include_legends_choice",
                help="Include or omit the LD legend images in the exported PNG/PDF.",
            )
            include_legends = _export_legend_choice == "Include"
            st.session_state["export_include_legends"] = include_legends

        current_signature = (
            output_format,
            paper_preset,
            float(width_in),
            float(height_in),
            int(dpi),
            bool(include_legends),
            id(st.session_state.get("last_locus_fig")),
            id(st.session_state.get("last_compare_fig")),
        )

        prepare_export = st.button(
            "Prepare export file",
            type="secondary",
            key="prepare_locusblend_export",
        )

        if prepare_export:
            try:
                data, mime, filename = make_locusblend_export_bytes(
                    output_format=output_format,
                    width_in=float(width_in),
                    height_in=float(height_in),
                    dpi=int(dpi),
                    include_legends=bool(include_legends),
                )
                st.session_state["last_export_bytes"] = data
                st.session_state["last_export_mime"] = mime
                st.session_state["last_export_filename"] = filename
                st.session_state["last_export_signature"] = current_signature
                st.session_state.pop("last_export_error", None)
                st.success(f"Export prepared: {filename}")
            except Exception as e:
                st.session_state["last_export_error"] = str(e)
                st.session_state.pop("last_export_bytes", None)
                st.session_state.pop("last_export_mime", None)
                st.session_state.pop("last_export_filename", None)
                st.session_state.pop("last_export_signature", None)

        if st.session_state.get("last_export_error"):
            st.error(st.session_state["last_export_error"])

        if (
            st.session_state.get("last_export_bytes") is not None
            and st.session_state.get("last_export_signature") == current_signature
        ):
            st.download_button(
                "Download export",
                data=st.session_state["last_export_bytes"],
                file_name=st.session_state.get("last_export_filename", "locusblend_export.png"),
                mime=st.session_state.get("last_export_mime", "application/octet-stream"),
                key="download_locusblend_export",
                use_container_width=True,
            )
        elif st.session_state.get("last_export_bytes") is not None:
            st.info("Export options changed. Select Prepare export file again.")


def render_locusblend_welcome_card():
    """Render a simple non-modal welcome panel.

    This replaces the old automatic welcome dialog. It must not change
    session_state, must not trigger recomputation, and must not open a modal.
    """
    with st.container(border=True):
        st.markdown("### Welcome to LocusBlend")
        st.markdown(
            """
LocusBlend visualizes GWAS summary statistics with multi-variant linkage
disequilibrium (LD) coloring, gene-track annotation, and locus compare plots.

**Basic workflow:**

1. Choose a visualization mode at the top of the page.
2. Configure data, LD source, locus, index variants, and display settings in the sidebar.
3. Select **Update plot** to apply changes.

Sidebar edits are intentionally not applied until **Update plot** is selected.
"""
        )

        c1, c2, c3 = st.columns(3)
        c1.markdown("**Modes**  \nStandard single-index, two-index, or three-index LocusBlend")
        c2.markdown("**Data**  \nUse default datasets or upload CSV, TSV, or TXT summary-statistic files")
        c3.markdown("**Output**  \nInteract with plots or export the current cached figure")


def render_washu_header():
    """Render a fixed full-width WashU Medicine-style app header."""
    logo_path = ASSET_DIR / "washulogo.png"
    try:
        logo_uri = image_to_data_uri(logo_path) if logo_path.exists() else None
    except Exception:
        logo_uri = None

    if logo_uri:
        logo_html = f'<img src="{logo_uri}" alt="WashU Medicine" class="lb-washu-logo" />'
    else:
        logo_html = '<span class="lb-washu-logo-fallback">WashU Medicine</span>'

    st.markdown(
        f"""
<div class="lb-washu-global-header" role="banner">
  <a class="lb-washu-header-link" href="https://medicine.washu.edu/" target="_blank" rel="noopener noreferrer" aria-label="WashU Medicine home">
    <div class="lb-washu-global-inner">
      {logo_html}
    </div>
  </a>
</div>
        """,
        unsafe_allow_html=True,
    )


def render_documentation_page():
    """Render the full LocusBlend user guide in the main content area."""
    st.markdown("## Documentation")
    st.markdown(
        '<div class="lb-doc-note">Use this page as the persistent LocusBlend reference while keeping the visualizer focused on plots and controls.</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        """
## Overview

**LocusBlend** visualizes GWAS summary statistics with multi-variant linkage
disequilibrium (LD) coloring, gene-track annotation, and locus compare plots.

The app supports:

* **Standard locus zoom**: one index variant
* **Two-index LocusBlend**: two index variants with blended LD coloring
* **Three-index LocusBlend**: three index variants with blended LD coloring

## Quick start

1. Select a mode at the top of the page.
2. Open the sidebar sections and configure:

   * datasets
   * LD source
   * locus window
   * index variants
   * plot display
   * gene track
   * locus compare
3. Select **Update plot**.
4. Use hover / zoom / pan in Plotly, or use **Export / Print current figure**.

## Important interaction model

Changing sidebar widgets reruns the Streamlit script, but it does **not**
recompute LD or rebuild the figures immediately.

The plot updates only when:

* **Update plot** is selected, or
* the app is loaded for the first time and no cached figure exists yet.

This design prevents expensive LD / PLINK computation from running after every
small parameter edit.

## Input file requirements

Supported uploaded summary-statistic file formats: `.csv`, `.tsv`, `.txt`,
`.csv.gz`, `.tsv.gz`, and `.txt.gz`.

Required genome build: GRCh38 / hg38. LocusBlend does not perform liftover. If
your summary-statistic file uses GRCh37 / hg19 or hg18 coordinates, lift it over
to GRCh38 / hg38 before upload. Otherwise internal 1000G LD matching, gene
annotation, locus windows, and index-variant matching may be incorrect or fail.

Required columns:

* `CHR`
* `BP`
* `P`
* `A1`
* `A2`

Recommended columns:

* `rsid`
* `BETA`
* `SE`
* `A1FREQ`
* `N`

Accepted aliases:

| Concept            | Accepted column names                               |
| ------------------ | --------------------------------------------------- |
| Chromosome         | `CHR`, `chr`, `#chr`, `chrom`, `chromosome`         |
| Base-pair position | `BP`, `bp`, `pos`, `position`, `base_pair_location` |
| P value            | `P`, `p`, `pval`, `pvalue`, `P-value`               |
| Variant ID         | `rsid`, `RSID`, `SNP`, `MarkerName`, `ID`           |
| Effect allele      | `A1`, `EA`, `effect_allele`, `ALLELE1`              |
| Other allele       | `A2`, `NEA`, `other_allele`, `ALLELE0`              |
| Effect size        | `BETA`, `beta`, `Effect`, `estimate`                |
| Standard error     | `SE`, `StdErr`, `stderr`                            |
| Allele frequency   | `A1FREQ`, `EAF`, `MAF`, `freq`                      |
| Sample size        | `N`, `n`, `samplesize`                              |

Example:

| CHR |       BP | rsid       |    P | A1 | A2 | BETA |   SE |
| --- | -------: | ---------- | ---: | -- | -- | ---: | ---: |
| 14  | 73238768 | rs11159021 | 1e-6 | A  | G  | 0.12 | 0.03 |

Supported chromosomes are **1-22 and X**. Chromosome X may be provided in
uploaded data as `X`, `chrX`, `23`, or `chr23`; the app normalizes these values
to `X` internally. X-based allele IDs use the same `CHR:BP:A1:A2` pattern as
autosomes, for example `X:13189:A:G`.

## LD source options

### Internal 1000G reference

No upload is required. Select one of the five 1000 Genomes super-populations:
AFR, AMR, EAS, EUR, or SAS. The selected ancestry determines which
chromosome-specific PLINK binary reference panel is used for LD calculation.
The default is EUR.

Uploaded summary-statistic coordinates must use GRCh38 / hg38; no liftover is
performed. Each selected ancestry requires local PLINK files in the app `data/`
directory using this prefix pattern:

`1000g_{AFR|AMR|EAS|EUR|SAS}_hg38_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_include_MHC_ch{chrom}`

### Uploaded LD matrix

Use a square R2 matrix with SNP IDs as both row and column names. Cells should
contain pairwise R2 values from 0 to 1.

Example:

| SNP | rs1 | rs2 | rs3 |
| --- | --: | --: | --: |
| rs1 | 1.0 | 0.7 | 0.1 |
| rs2 | 0.7 | 1.0 | 0.3 |
| rs3 | 0.1 | 0.3 | 1.0 |

### Uploaded LD long table

Required columns:

* `SNP_A`
* `SNP_B`
* `R2`

Accepted aliases:

* `SNP_A`: `SNP_A`, `SNP1`, `ID1`
* `SNP_B`: `SNP_B`, `SNP2`, `ID2`
* `R2`: `R2`, `r2`, `rsq`

The long table must include self-pair rows where `SNP_A == SNP_B` and `R2 == 1`.
These rows allow the app to distinguish variants present in the LD reference
from variants that are missing.

Example:

| SNP_A | SNP_B |  R2 |
| ----- | ----- | --: |
| rs1   | rs1   | 1.0 |
| rs2   | rs2   | 1.0 |
| rs1   | rs2   | 0.7 |

## Locus settings

Use the sidebar to choose:

* chromosome (1-22 or X)
* center base-pair position
* window size in kb

The window is interpreted as `center +/- window_kb`.

## Index variants

Manual input:

* Enter one, two, or three variant IDs depending on the active mode.

Auto-select:

* The app can select independent index variants from either dataset 1 (top) or
  dataset 2 (bottom) using LD clumping.
* Auto-selection uses LD clumping at r2 = 0.01.
* The number of selected variants depends on the active mode.

## Plot display

Use plot display settings to control:

* dataset 1 (top) and dataset 2 (bottom) plot titles
* y-axis limits
* combined plot height
* vertical spacing
* recombination-rate display
* recombination y-axis maximum

## Gene track

Gene display modes:

* `protein_coding`
* `all`

Gene highlighting:

* Enter comma- or semicolon-separated gene names.
* Matching is case-insensitive.
* Highlighted genes are drawn more prominently in the gene track.

Example:

`PSEN1, PAPLN, HEATR4`

## Locus compare

The locus compare panel supports:

* **Three separate compare plots**: one panel per active index variant
* **Single blended compare plot**: all variants overlaid in one blended plot

The number of active panels depends on the active mode:

* Standard: variant 1
* Two-index: variants 1 and 2
* Three-index: variants 1, 2, and 3

## Export / Print current figure

The export panel exports the currently rendered cached plot. It does not apply
pending sidebar edits.

Supported formats:

* PNG
* PDF

Default page size:

* US Letter portrait, 8.5 x 11 inches

You can also choose:

* US Letter landscape
* Custom page size
* DPI
* whether to include the legend images

If you change export options after preparing a file, select **Prepare export
file** again before downloading.

## Troubleshooting

**I changed a parameter but the plot did not change.**
Select **Update plot**. Sidebar edits are intentionally not applied immediately.

**I changed LD source but the plot did not change.**
First select **Apply LD source**, then select **Update plot**.

**A variant is shown as a grey X.**
The variant was not found in the active LD reference or uploaded LD universe.

**Uploaded LD long table fails.**
Check that self-pair rows are included: `SNP_A == SNP_B`, `R2 == 1`.

**Export fails.**
Server-side export requires Plotly static image export dependencies. The same
Python environment that runs Streamlit needs `plotly`, `kaleido`, and `pillow`.
Kaleido v1 may also require server-side Chrome or Chromium.
"""
    )


def render_page_switcher():
    """Render a lightweight top-level page switcher."""
    return st.radio(
        "Page view",
        options=["Visualizer", "Documentation"],
        index=0,
        horizontal=True,
        label_visibility="collapsed",
        key="page_view",
    )


def load_locusblend_page_icon():
    icon_path = ASSET_DIR / "WashU-SHIELD-Red_RGB.png"
    if icon_path.exists():
        try:
            from PIL import Image

            return Image.open(icon_path)
        except Exception:
            pass
    return "🧬"


st.set_page_config(
    page_title="LocusBlend",
    page_icon=load_locusblend_page_icon(),
    layout="wide",
    initial_sidebar_state="expanded",
)
set_streamlit_chrome_minimal()
inject_locusblend_css()
render_washu_header()
st.title("LocusBlend")
st.caption("Flexible multi-index regional visualization of genomic association signals")
page_view = render_page_switcher()

if page_view == "Documentation":
    render_documentation_page()
    st.stop()

render_locusblend_welcome_card()

st.session_state.setdefault("active_locusblend_mode", "Three-index LocusBlend")

_locusblend_modes = ["Standard locus zoom", "Two-index LocusBlend", "Three-index LocusBlend"]

if (
    "pending_locusblend_mode" not in st.session_state
    or st.session_state["pending_locusblend_mode"] not in _locusblend_modes
):
    st.session_state["pending_locusblend_mode"] = st.session_state["active_locusblend_mode"]

st.radio(
    "Mode",
    _locusblend_modes,
    horizontal=True,
    key="pending_locusblend_mode",
    label_visibility="collapsed",
)

if st.session_state["pending_locusblend_mode"] != st.session_state["active_locusblend_mode"]:
    st.info("Select **Update plot** to apply the selected mode.")

active_mode = st.session_state["active_locusblend_mode"]

# The old automatic welcome dialog was replaced by the non-modal
# welcome card and the top-level Documentation view.
progress_bar = None
status_box = None

try:
    render_sidebar_app_header()
    render_sidebar_visual_strip()

    with st.sidebar.expander("1. Data upload", expanded=False):
        st.markdown("Upload dataset 1 (top) and dataset 2 (bottom) summary-statistic files, or use the default example datasets.")
        st.caption("Full input file format details are in the top Documentation view.")
        uploaded_top = st.file_uploader(
            "Upload dataset 1 (top) summary-statistic file",
            type=["csv", "tsv", "txt", "gz"],
            key="uploaded_top",
        )
        st.caption(
            "Required: CHR, BP, P, A1, A2  -  Recommended: rsid, BETA, SE  -  "
            "Genome build: GRCh38 / hg38 only; no liftover is performed.  -  "
            "Supported: .csv, .tsv, .txt, .csv.gz, .tsv.gz, .txt.gz"
        )
        uploaded_bottom = st.file_uploader(
            "Upload dataset 2 (bottom) summary-statistic file",
            type=["csv", "tsv", "txt", "gz"],
            key="uploaded_bottom",
        )
        st.caption(
            "Required: CHR, BP, P, A1, A2  -  Recommended: rsid, BETA, SE  -  "
            "Genome build: GRCh38 / hg38 only; no liftover is performed.  -  "
            "Supported: .csv, .tsv, .txt, .csv.gz, .tsv.gz, .txt.gz"
        )

    with st.sidebar.expander("2. LD source", expanded=False):
        st.markdown("Choose the active LD source.")
        st.caption("Full LD format details are in the top Documentation view.")
        _ld_source_options = [
            "Use internal 1000G reference",
            "Upload LD matrix",
            "Upload LD long table",
        ]
        # Sticky/applied state — used by the downstream pipeline
        st.session_state.setdefault("active_ld_source", _ld_source_options[0])
        st.session_state.setdefault("active_ld_matrix_bytes", None)
        st.session_state.setdefault("active_ld_matrix_name", "")
        st.session_state.setdefault("active_ld_long_bytes", None)
        st.session_state.setdefault("active_ld_long_name", "")
        if "internal_1000g_ancestry" not in st.session_state:
            st.session_state["internal_1000g_ancestry"] = INTERNAL_1000G_DEFAULT_ANCESTRY
        else:
            normalized_ancestry = normalize_internal_1000g_ancestry(
                st.session_state.get("internal_1000g_ancestry")
            )
            if normalized_ancestry != st.session_state.get("internal_1000g_ancestry"):
                st.session_state["internal_1000g_ancestry"] = normalized_ancestry
        # Pending UI state — bound to the radio widget; only copies to active on Apply
        st.session_state.setdefault("pending_ld_source", st.session_state["active_ld_source"])

        pending_ld_source = st.radio(
            "LD source",
            _ld_source_options,
            index=_ld_source_options.index(st.session_state["pending_ld_source"])
                if st.session_state["pending_ld_source"] in _ld_source_options else 0,
            key="pending_ld_source",
            label_visibility="collapsed",
            captions=[
                "Embedded 1000G PLINK reference (default)",
                "Named R² matrix, SNP IDs as row and column names",
                "Long-format pairwise table with SNP_A, SNP_B, R2",
            ],
        )

        _matrix_widget = None
        _long_widget = None
        if pending_ld_source == "Use internal 1000G reference":
            st.selectbox(
                "1000G ancestry",
                options=INTERNAL_1000G_ANCESTRY_OPTIONS,
                index=INTERNAL_1000G_ANCESTRY_OPTIONS.index(
                    st.session_state["internal_1000g_ancestry"]
                ),
                key="internal_1000g_ancestry",
                format_func=format_internal_1000g_ancestry_option,
                help=(
                    "Select the 1000 Genomes super-population used for LD calculation. "
                    "Uploaded coordinates must be GRCh38 / hg38."
                ),
            )
        elif pending_ld_source == "Upload LD matrix":
            _matrix_widget = st.file_uploader(
                "Upload LD matrix",
                type=["tsv", "csv", "txt", "ld"],
                key="uploaded_ld_matrix",
            )
            st.caption("Square R² matrix with SNP IDs as row and column names.")
        elif pending_ld_source == "Upload LD long table":
            _long_widget = st.file_uploader(
                "Upload LD long table",
                type=["tsv", "csv", "txt", "ld"],
                key="uploaded_ld_long",
            )
            st.caption("Columns: SNP_A, SNP_B, R2. Include self-pairs with R2=1.")

        _apply_clicked = st.button("Apply LD source", key="apply_ld_source")
        if _apply_clicked:
            if pending_ld_source == "Upload LD matrix" and _matrix_widget is None:
                st.error(
                    "LD matrix file is required for 'Upload LD matrix'. "
                    "Please upload a file first, then select 'Apply LD source'."
                )
            elif pending_ld_source == "Upload LD long table" and _long_widget is None:
                st.error(
                    "LD long table file is required for 'Upload LD long table'. "
                    "Please upload a file first, then select 'Apply LD source'."
                )
            else:
                st.session_state["active_ld_source"] = pending_ld_source
                if pending_ld_source == "Upload LD matrix" and _matrix_widget is not None:
                    st.session_state["active_ld_matrix_bytes"] = _matrix_widget.getvalue()
                    st.session_state["active_ld_matrix_name"] = _matrix_widget.name
                if pending_ld_source == "Upload LD long table" and _long_widget is not None:
                    st.session_state["active_ld_long_bytes"] = _long_widget.getvalue()
                    st.session_state["active_ld_long_name"] = _long_widget.name
                if pending_ld_source == "Use internal 1000G reference":
                    _ancestry = normalize_internal_1000g_ancestry(st.session_state["internal_1000g_ancestry"])
                    st.success(f"LD source applied: Internal 1000G reference ({_ancestry})")
                else:
                    _applied_name = (
                        st.session_state["active_ld_matrix_name"]
                        if pending_ld_source == "Upload LD matrix"
                        else st.session_state["active_ld_long_name"]
                    )
                    st.success(f"LD source applied: {pending_ld_source} ({_applied_name})")

    ld_source = st.session_state["active_ld_source"]
    active_internal_1000g_ancestry = normalize_internal_1000g_ancestry(
        st.session_state.get("internal_1000g_ancestry", INTERNAL_1000G_DEFAULT_ANCESTRY)
    )
    if ld_source == "Use internal 1000G reference":
        st.sidebar.caption(f"Active LD source: Internal 1000G reference ({active_internal_1000g_ancestry})")
    else:
        st.sidebar.caption(f"Active LD source: {ld_source}")

    update_progress(progress_bar, status_box, 2, "Preparing app...")

    if uploaded_top is not None:
        update_progress(progress_bar, status_box, 5, f"Reading dataset 1 (top): {uploaded_top.name}")
        df_top = load_locus_csv_uploaded(uploaded_top.getvalue(), uploaded_top.name)
        top_source_name = uploaded_top.name
    else:
        update_progress(progress_bar, status_box, 5, "Reading default dataset 1 (top)...")
        df_top = load_locus_csv(str(DATA_DIR / "merged_female_withld_PSEN1.csv"))
        top_source_name = "merged_female_withld_PSEN1.csv"

    if uploaded_bottom is not None:
        update_progress(progress_bar, status_box, 10, f"Reading dataset 2 (bottom): {uploaded_bottom.name}")
        df_bottom = load_locus_csv_uploaded(uploaded_bottom.getvalue(), uploaded_bottom.name)
        bottom_source_name = uploaded_bottom.name
    else:
        update_progress(progress_bar, status_box, 10, "Reading default dataset 2 (bottom)...")
        df_bottom = load_locus_csv(str(DATA_DIR / "Nedelec_monocytes_PSEN1_example.csv"))
        bottom_source_name = "Nedelec_monocytes_PSEN1_example.csv"

    st.session_state.setdefault("chrom", "14")
    st.session_state.setdefault("index_selection_method", "Manual input")
    st.session_state.setdefault("bp", 73238768)
    st.session_state.setdefault("window_kb", 500)
    st.session_state.setdefault("rsid", "rs11159021")
    st.session_state.setdefault("rsid_2", "rs3742825")
    st.session_state.setdefault("rsid_3", "rs74986264")
    st.session_state.setdefault("title_top", "Alzheimer's disease Female GWAS")
    st.session_state.setdefault("title_bottom", "Monocyte PSEN1 eQTL")
    st.session_state.setdefault("show_recomb", True)
    st.session_state.setdefault("combined_height", 980)
    st.session_state.setdefault("vertical_spacing", 0.04)
    st.session_state.setdefault("recomb_max", 100)
    st.session_state.setdefault("gtf_path", str(DATA_DIR / "gencode.v49.annotation.chr14.gtf.gz"))
    st.session_state.setdefault("gene_display_mode", "protein_coding")
    st.session_state.setdefault("gene_track_gap", 30000)
    st.session_state.setdefault("highlight_genes", "")
    st.session_state.setdefault("compare_size", 560)
    st.session_state.setdefault("locuscompare_mode", "Three separate compare plots")
    st.session_state.setdefault(
        "ld_bfile_prefix",
        str(get_internal_1000g_prefix("14", INTERNAL_1000G_DEFAULT_ANCESTRY))
    )
    st.session_state.setdefault("plink_path", str(BIN_DIR / "plink"))

    upload_signature = make_uploaded_dataset_signature(uploaded_top, uploaded_bottom)
    if upload_signature is not None:
        upload_sync_result = infer_uploaded_locus_sync(
            df_top,
            df_bottom,
            uploaded_top,
            uploaded_bottom,
        )
        apply_uploaded_dataset_sync_if_new(
            upload_signature,
            upload_sync_result,
            df_top=df_top,
            df_bottom=df_bottom,
            uploaded_top=uploaded_top,
            uploaded_bottom=uploaded_bottom,
        )

    active_upload_sync_result = st.session_state.get("last_upload_sync_result")
    if (
        upload_signature is not None
        and st.session_state.get("last_uploaded_dataset_signature") == upload_signature
        and active_upload_sync_result
    ):
        _upload_sync_message = format_upload_sync_message(active_upload_sync_result)
        if _upload_sync_message:
            if active_upload_sync_result.get("status") == "warning":
                st.sidebar.warning(_upload_sync_message)
            else:
                st.sidebar.info(_upload_sync_message)

    current_bp = int(st.session_state["bp"])
    current_window_bp = int(st.session_state["window_kb"] * 1000)

    update_progress(progress_bar, status_box, 15, "Estimating default y-axis ranges...")

    _est_chrom = normalize_chrom(st.session_state.get("chrom", "14"))
    window_top = df_top[
        chrom_mask(df_top, _est_chrom)
        & (df_top["BP"] >= current_bp - current_window_bp)
        & (df_top["BP"] <= current_bp + current_window_bp)
    ].copy()
    if len(window_top) == 0:
        window_top = df_top.loc[chrom_mask(df_top, _est_chrom)].copy()
        if len(window_top) == 0:
            window_top = df_top.copy()
    default_max_ylim_top = int(np.ceil((-np.log10(window_top["P"])).max() + 1))

    window_bottom = df_bottom[
        chrom_mask(df_bottom, _est_chrom)
        & (df_bottom["BP"] >= current_bp - current_window_bp)
        & (df_bottom["BP"] <= current_bp + current_window_bp)
    ].copy()
    if len(window_bottom) == 0:
        window_bottom = df_bottom.loc[chrom_mask(df_bottom, _est_chrom)].copy()
        if len(window_bottom) == 0:
            window_bottom = df_bottom.copy()
    default_max_ylim_bottom = int(np.ceil((-np.log10(window_bottom["P"])).max() + 1))

    if "max_ylim_top" not in st.session_state:
        st.session_state["max_ylim_top"] = default_max_ylim_top
    if "max_ylim_bottom" not in st.session_state:
        st.session_state["max_ylim_bottom"] = default_max_ylim_bottom

    with st.sidebar.form("plot_controls"):
        submitted_top = False

        with st.expander("3. Locus window", expanded=False):
            _chrom_options = get_supported_chromosomes()
            _current_chrom = normalize_chrom(st.session_state.get("chrom", "14"))
            if not is_supported_chrom(_current_chrom):
                _current_chrom = "14"
            st.session_state["chrom"] = _current_chrom
            chrom = st.selectbox(
                "Chromosome",
                _chrom_options,
                index=_chrom_options.index(_current_chrom),
                key="chrom",
            )
            bp = st.number_input("Center position (bp)", value=st.session_state["bp"], step=1, key="bp")
            window_kb = st.number_input("Window size (+/- kb)", value=st.session_state["window_kb"], step=50, min_value=1, key="window_kb")
            window_bp = int(window_kb * 1000)

        with st.expander("4. Index variants", expanded=False):
            _index_method_aliases = {
                "Auto-select by LD clumping from " + "top " + "dataset": "Auto-select by LD clumping from dataset 1 (top)",
                "Auto-select by LD clumping from " + "bottom " + "dataset": "Auto-select by LD clumping from dataset 2 (bottom)",
            }
            if st.session_state.get("index_selection_method") in _index_method_aliases:
                st.session_state["index_selection_method"] = _index_method_aliases[
                    st.session_state["index_selection_method"]
                ]
            _index_methods = [
                "Manual input",
                "Auto-select by LD clumping from dataset 1 (top)",
                "Auto-select by LD clumping from dataset 2 (bottom)",
            ]
            _index_method_idx = _index_methods.index(st.session_state["index_selection_method"]) \
                if st.session_state["index_selection_method"] in _index_methods else 0
            index_selection_method = st.selectbox(
                "Index selection method",
                _index_methods,
                index=_index_method_idx,
                key="index_selection_method",
            )
            st.caption(
                "Manual: use the variant IDs below. "
                "Auto: the app selects the most significant independent variants "
                "within the locus window using LD clumping at r² = 0.01."
            )

            st.markdown("##### Manual index inputs")
            if index_selection_method != "Manual input":
                st.caption("Manual index fields are ignored when auto-selection is enabled.")
            rsid = st.text_input("Index variant 1", value=st.session_state["rsid"], key="rsid")
            _pending_mode = st.session_state.get("pending_locusblend_mode", "Three-index LocusBlend")
            if _pending_mode != "Standard locus zoom":
                rsid_2 = st.text_input("Index variant 2", value=st.session_state["rsid_2"], key="rsid_2")
            else:
                rsid_2 = ""
            if _pending_mode == "Three-index LocusBlend":
                rsid_3 = st.text_input("Index variant 3", value=st.session_state["rsid_3"], key="rsid_3")
            else:
                rsid_3 = ""

        with st.expander("5. Plot display", expanded=False):
            title_top = st.text_input("Dataset 1 (top) plot title", value=st.session_state["title_top"], key="title_top")
            title_bottom = st.text_input("Dataset 2 (bottom) plot title", value=st.session_state["title_bottom"], key="title_bottom")
            max_ylim_top = st.number_input("Y-axis max (top)", value=int(st.session_state["max_ylim_top"]), step=1, key="max_ylim_top")
            max_ylim_bottom = st.number_input("Y-axis max (bottom)", value=int(st.session_state["max_ylim_bottom"]), step=1, key="max_ylim_bottom")
            combined_height = st.number_input("Combined plot height", value=int(st.session_state["combined_height"]), step=40, key="combined_height")
            vertical_spacing = st.number_input("Vertical spacing", value=float(st.session_state["vertical_spacing"]), step=0.01, format="%.3f", key="vertical_spacing")
            _recomb_default = bool(st.session_state.get("show_recomb", True))
            _recomb_choice = st.selectbox(
                "Recombination rate",
                options=["Show", "Hide"],
                index=0 if _recomb_default else 1,
                key="show_recomb_choice",
                help="Show or hide the recombination-rate track on locus zoom plots.",
            )
            show_recomb = _recomb_choice == "Show"
            st.session_state["show_recomb"] = show_recomb
            recomb_max = st.number_input("Recombination Y max", value=int(st.session_state["recomb_max"]), step=10, key="recomb_max")

        with st.expander("6. Gene track", expanded=False):
            gene_display_mode = st.selectbox(
                "Gene display mode",
                ["protein_coding", "all"],
                index=["protein_coding", "all"].index(st.session_state["gene_display_mode"]),
                key="gene_display_mode"
            )
            gene_track_gap = st.number_input("Gene track min gap (bp)", value=int(st.session_state["gene_track_gap"]), step=5000, key="gene_track_gap")
            st.text_input("Highlight genes", value=st.session_state.get("highlight_genes", ""), placeholder="PSEN1, PAPLN, HEATR4", key="highlight_genes")
            st.caption("Comma or semicolon separated gene names (case-insensitive)")

        with st.expander("7. Locus compare", expanded=False):
            compare_size = st.number_input("Compare plot size", value=int(st.session_state["compare_size"]), step=20, key="compare_size")
            _locuscompare_options = ["Three separate compare plots", "Single blended compare plot"]
            locuscompare_mode = st.selectbox(
                "Locus compare mode",
                _locuscompare_options,
                index=_locuscompare_options.index(st.session_state["locuscompare_mode"])
                    if st.session_state["locuscompare_mode"] in _locuscompare_options else 0,
                key="locuscompare_mode"
            )
            st.caption("Three separate: one panel per index variant | Single blended: all variants overlaid in one plot")

        submitted_bottom = st.form_submit_button(
            "Update plot",
            use_container_width=True,
            type="primary",
        )

        submitted = submitted_top or submitted_bottom

    plink_path = st.session_state["plink_path"]
    ld_bfile_prefix = st.session_state["ld_bfile_prefix"]
    gtf_path = st.session_state["gtf_path"]

    st.caption(f"Dataset 1 (top) source: {top_source_name} | Dataset 2 (bottom) source: {bottom_source_name}")

    # Widget changes (radio, file uploader, number inputs) trigger Streamlit
    # reruns but should not trigger the LD pipeline or figure rebuild — those
    # only happen when the user selects the "Update plot" button, or on the
    # very first load (no cached figure yet).
    should_compute = submitted or "last_locus_fig" not in st.session_state

    if should_compute:
        st.session_state["active_locusblend_mode"] = st.session_state["pending_locusblend_mode"]
        progress_bar = st.progress(0)
        status_box = st.empty()
        update_progress(progress_bar, status_box, 5, "Preparing LD pipeline...")

        active_mode = st.session_state["active_locusblend_mode"]
        active_internal_1000g_ancestry = normalize_internal_1000g_ancestry(
            st.session_state.get("internal_1000g_ancestry", INTERNAL_1000G_DEFAULT_ANCESTRY)
        )
        ld_labels = get_ld_reference_labels(ld_source, active_internal_1000g_ancestry)

        selected_chrom = normalize_chrom(chrom)
        if not is_supported_chrom(selected_chrom):
            raise ValueError(
                f"Unsupported chromosome {chrom}. Supported chromosomes are 1-22 and X."
            )

        # Validate selected chromosome data presence
        _top_mask = chrom_mask(df_top, selected_chrom)
        _bottom_mask = chrom_mask(df_bottom, selected_chrom)
        if not _top_mask.any() and not _bottom_mask.any():
            raise ValueError(
                f"No rows found for chromosome {selected_chrom} in either dataset. "
                "Check the chromosome selector or upload data for this chromosome."
            )

        # Filter to selected chromosome for the LD pipeline
        df_top_chr = df_top.loc[_top_mask].copy()
        df_bottom_chr = df_bottom.loc[_bottom_mask].copy()

        # Resolve references dynamically
        if ld_source == "Use internal 1000G reference":
            ld_bfile_prefix = get_internal_bfile_prefix_for_chrom(selected_chrom, active_internal_1000g_ancestry)
        gtf_path = get_gtf_path_for_chrom(selected_chrom)

        chrom_for_ld = selected_chrom
        ld_bp_start = int(bp - window_bp)
        ld_bp_end = int(bp + window_bp)

        ld_long_table = None
        ld_matrix_table = None
        ld_uploaded_universe = None

        if ld_source == "Upload LD long table":
            _long_bytes = st.session_state.get("active_ld_long_bytes")
            _long_name = st.session_state.get("active_ld_long_name") or "ld_long"
            if _long_bytes is None:
                raise ValueError(
                    "Active LD source is 'Upload LD long table' but no LD file has been applied. "
                    "Upload a file and select 'Apply LD source' in the sidebar."
                )
            update_progress(progress_bar, status_box, 16, f"Parsing uploaded LD long table: {_long_name}")
            ld_long_table, ld_uploaded_universe = read_uploaded_ld_long(_long_bytes, _long_name)
        elif ld_source == "Upload LD matrix":
            _matrix_bytes = st.session_state.get("active_ld_matrix_bytes")
            _matrix_name = st.session_state.get("active_ld_matrix_name") or "ld_matrix"
            if _matrix_bytes is None:
                raise ValueError(
                    "Active LD source is 'Upload LD matrix' but no LD file has been applied. "
                    "Upload a file and select 'Apply LD source' in the sidebar."
                )
            update_progress(progress_bar, status_box, 16, f"Parsing uploaded LD matrix: {_matrix_name}")
            ld_matrix_table, ld_uploaded_universe = read_uploaded_ld_matrix(_matrix_bytes, _matrix_name)

        if ld_source == "Use internal 1000G reference":
            update_progress(progress_bar, status_box, 18, f"Computing LD from 1000G {active_internal_1000g_ancestry} with PLINK...")
            bim_ref = load_reference_bim(ld_bfile_prefix)
            df_top_ref = attach_reference_snp_two_pass(df_top_chr, bim_ref)
            df_bottom_ref = attach_reference_snp_two_pass(df_bottom_chr, bim_ref)
        else:
            update_progress(progress_bar, status_box, 18, "Matching summary stats against uploaded LD universe...")
            df_top_ref = attach_uploaded_ld_keys(df_top_chr, ld_uploaded_universe)
            df_bottom_ref = attach_uploaded_ld_keys(df_bottom_chr, ld_uploaded_universe)

        # Index variant resolution — manual or auto
        idx1_label_for_plot = rsid
        idx2_label_for_plot = rsid_2
        idx3_label_for_plot = rsid_3
        auto_index_table = None
        auto_index_summary = ""

        if index_selection_method == "Manual input":
            idx1_row = resolve_index_variant_from_input(df_top_ref, df_bottom_ref, rsid)
            idx1_ref = idx1_row["REF_SNP"]

            idx2_ref = None
            idx3_ref = None

            if active_mode != "Standard locus zoom" and str(rsid_2).strip() != "":
                idx2_row = resolve_index_variant_from_input(df_top_ref, df_bottom_ref, rsid_2)
                idx2_ref = idx2_row["REF_SNP"]

            if active_mode == "Three-index LocusBlend" and str(rsid_3).strip() != "":
                idx3_row = resolve_index_variant_from_input(df_top_ref, df_bottom_ref, rsid_3)
                idx3_ref = idx3_row["REF_SNP"]

        else:
            update_progress(progress_bar, status_box, 19, "Running auto index selection by LD clumping...")
            required_n = get_required_n_indices(active_mode)

            if index_selection_method == "Auto-select by LD clumping from dataset 1 (top)":
                _source_ref = df_top_ref
                _source_label = "dataset 1 (top)"
            else:
                _source_ref = df_bottom_ref
                _source_label = "dataset 2 (bottom)"

            auto_selected, auto_summary = auto_select_index_variants_by_clumping(
                df_source_ref=_source_ref,
                selected_chrom=selected_chrom,
                center_bp=bp,
                window_bp=window_bp,
                required_n=required_n,
                active_ld_source=ld_source,
                bfile_prefix=ld_bfile_prefix if ld_source == "Use internal 1000G reference" else None,
                plink_path=plink_path,
                ld_long_table=ld_long_table,
                ld_matrix_table=ld_matrix_table,
                clump_r2=0.01,
            )

            idx1_ref = auto_selected.iloc[0]["REF_SNP"]
            idx1_label_for_plot = str(auto_selected.iloc[0]["DISPLAY_ID"])
            idx2_ref = None
            idx3_ref = None
            idx2_label_for_plot = ""
            idx3_label_for_plot = ""

            if required_n >= 2 and len(auto_selected) >= 2:
                idx2_ref = auto_selected.iloc[1]["REF_SNP"]
                idx2_label_for_plot = str(auto_selected.iloc[1]["DISPLAY_ID"])
            if required_n >= 3 and len(auto_selected) >= 3:
                idx3_ref = auto_selected.iloc[2]["REF_SNP"]
                idx3_label_for_plot = str(auto_selected.iloc[2]["DISPLAY_ID"])

            auto_index_table = auto_selected.copy()
            auto_index_table["_source"] = _source_label
            auto_index_summary = (
                f"Auto-selected from the {_source_label} within "
                f"chr{selected_chrom}:{auto_summary['start_bp']}-{auto_summary['end_bp']} "
                f"using LD clumping at r² = {auto_summary['clump_r2']}. "
                f"({auto_summary['n_candidates']:,} candidates → "
                f"{auto_summary['n_selected']} independent variants)"
            )

        top_window = df_top_ref[
            chrom_mask(df_top_ref, selected_chrom)
            & (df_top_ref["BP"] >= ld_bp_start)
            & (df_top_ref["BP"] <= ld_bp_end)
        ][["REF_SNP", "CHR", "BP"]].copy()

        bottom_window = df_bottom_ref[
            chrom_mask(df_bottom_ref, selected_chrom)
            & (df_bottom_ref["BP"] >= ld_bp_start)
            & (df_bottom_ref["BP"] <= ld_bp_end)
        ][["REF_SNP", "CHR", "BP"]].copy()

        window_union = pd.concat([top_window, bottom_window], ignore_index=True).drop_duplicates("REF_SNP")
        window_snps = tuple(sorted(window_union["REF_SNP"].dropna().astype(str).unique().tolist()))

        if ld_source == "Use internal 1000G reference":
            ld_maps, index_status, ref_snps = compute_ld_maps_with_plink(
                bfile_prefix=ld_bfile_prefix,
                chrom=chrom_for_ld,
                start=ld_bp_start,
                end=ld_bp_end,
                window_snps=window_snps,
                idx1_ref=idx1_ref,
                idx2_ref=idx2_ref,
                idx3_ref=idx3_ref,
                plink_path=plink_path
            )
        elif ld_source == "Upload LD long table":
            ld_maps, index_status, ref_snps = compute_ld_maps_from_uploaded_long(
                ld_long=ld_long_table,
                window_snps=window_snps,
                idx1_ref=idx1_ref,
                idx2_ref=idx2_ref,
                idx3_ref=idx3_ref,
            )
        else:
            ld_maps, index_status, ref_snps = compute_ld_maps_from_uploaded_matrix(
                ld_matrix=ld_matrix_table,
                window_snps=window_snps,
                idx1_ref=idx1_ref,
                idx2_ref=idx2_ref,
                idx3_ref=idx3_ref,
            )

        ld_annot = build_ld_annot_for_window(
            window_union_df=window_union,
            ld_maps=ld_maps,
            ref_snps=ref_snps,
            idx1_ref=idx1_ref,
            idx2_ref=idx2_ref,
            idx3_ref=idx3_ref
        )

        df_top_plot = merge_ld_annot(df_top_ref, ld_annot)
        df_bottom_plot = merge_ld_annot(df_bottom_ref, ld_annot)

        _src = ld_labels["source_name"]
        _status_parts = [
            f"variant 1: {('found in ' + _src) if index_status['variant 1'] else ('not in ' + _src + ' -> grey X')}",
        ]
        if active_mode != "Standard locus zoom":
            _status_parts.append(
                f"variant 2: {('found in ' + _src) if index_status['variant 2'] else ('not in ' + _src + ' -> grey X')}",
            )
        if active_mode == "Three-index LocusBlend":
            _status_parts.append(
                f"variant 3: {('found in ' + _src) if index_status['variant 3'] else ('not in ' + _src + ' -> grey X')}",
            )
        ld_status_caption = " | ".join(_status_parts)

        _ref_info = ""
        if ld_source == "Use internal 1000G reference":
            _ref_info = f"{active_internal_1000g_ancestry} {format_chrom_label(selected_chrom)}"
        ld_status_caption += f"  |  Chromosome: {format_chrom_label(selected_chrom)}"
        if _ref_info:
            ld_status_caption += f"  |  Ref: 1000G {_ref_info}"

        update_progress(progress_bar, status_box, 25, "Building top locus plot...")
        fig_top, n_top, bp_start_top, bp_end_top = get_plotly_locus_py(
            max_ylim=max_ylim_top,
            bp=bp,
            window_bp=window_bp,
            merged_female_withld=df_top_plot,
            merged_df=df_top_plot,
            idx1_ref=idx1_ref,
            idx2_ref=idx2_ref,
            idx3_ref=idx3_ref,
            idx1_label=idx1_label_for_plot,
            idx2_label=idx2_label_for_plot,
            idx3_label=idx3_label_for_plot,
            y_label=title_top,
            ld_labels=ld_labels,
            chrom=selected_chrom,
            show_recomb=show_recomb,
            bw_path=str(DATA_DIR / "recomb1000GAvg.bw"),
            locusblend_mode=active_mode,
        )

        update_progress(progress_bar, status_box, 40, "Building bottom locus plot...")
        fig_bottom, n_bottom, bp_start_bottom, bp_end_bottom = get_plotly_locus_py(
            max_ylim=max_ylim_bottom,
            bp=bp,
            window_bp=window_bp,
            merged_female_withld=df_bottom_plot,
            merged_df=df_bottom_plot,
            idx1_ref=idx1_ref,
            idx2_ref=idx2_ref,
            idx3_ref=idx3_ref,
            idx1_label=idx1_label_for_plot,
            idx2_label=idx2_label_for_plot,
            idx3_label=idx3_label_for_plot,
            y_label=title_bottom,
            ld_labels=ld_labels,
            chrom=selected_chrom,
            show_recomb=show_recomb,
            bw_path=str(DATA_DIR / "recomb1000GAvg.bw"),
            locusblend_mode=active_mode,
        )

        update_progress(progress_bar, status_box, 55, f"Loading gene track ({gene_display_mode})...")
        gene_bp_start = min(bp_start_top, bp_start_bottom)
        gene_bp_end = max(bp_end_top, bp_end_bottom)

        if gtf_path.endswith(".parquet"):
            genes_df = load_genes_from_table(
                chrom=selected_chrom,
                start=gene_bp_start,
                end=gene_bp_end,
                parquet_path=gtf_path,
                gene_display_mode=gene_display_mode
            )
        else:
            genes_df = load_genes_from_gtf(
                gtf_path=gtf_path,
                chrom=selected_chrom,
                start=gene_bp_start,
                end=gene_bp_end,
                gene_display_mode=gene_display_mode
            )

        update_progress(progress_bar, status_box, 70, "Assembling shared locus figure...")
        locus_fig = make_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=vertical_spacing,
            row_heights=[0.36, 0.36, 0.28],
            specs=[
                [{"secondary_y": True}],
                [{"secondary_y": True}],
                [{"secondary_y": False}]
            ],
            subplot_titles=(title_top, title_bottom, f"GENCODE gene track ({gene_display_mode})")
        )

        for tr in fig_top.data:
            tr_json = tr.to_plotly_json()
            tr_json.pop("xaxis", None)
            tr_json.pop("yaxis", None)
            if tr_json.get("type") == "scattergl":
                tr2 = go.Scattergl(**tr_json)
            else:
                tr2 = go.Scatter(**tr_json)
            is_secondary = getattr(tr, "mode", None) == "lines"
            locus_fig.add_trace(tr2, row=1, col=1, secondary_y=is_secondary)

        for tr in fig_bottom.data:
            tr_json = tr.to_plotly_json()
            tr_json.pop("xaxis", None)
            tr_json.pop("yaxis", None)
            if tr_json.get("type") == "scattergl":
                tr2 = go.Scattergl(**tr_json)
            else:
                tr2 = go.Scatter(**tr_json)
            is_secondary = getattr(tr, "mode", None) == "lines"
            locus_fig.add_trace(tr2, row=2, col=1, secondary_y=is_secondary)

        _raw_highlight = st.session_state.get("highlight_genes", "")
        highlight_names = {g.strip().lower() for g in re.split(r"[,;\s]+", _raw_highlight) if g.strip()}

        locus_fig, gene_track_rows = add_gene_track_to_subplot(
            locus_fig,
            genes_df,
            row=3,
            col=1,
            min_gap=gene_track_gap,
            highlight_names=highlight_names
        )

        locus_fig.update_yaxes(
            title_text="-log<sub>10</sub>(P)",
            range=list(fig_top.layout.yaxis.range),
            zeroline=False,
            row=1,
            col=1,
            secondary_y=False
        )
        locus_fig.update_yaxes(
            title_text="Recombination rate",
            range=[0, recomb_max],
            autorange=False,
            zeroline=False,
            showgrid=False,
            row=1,
            col=1,
            secondary_y=True
        )

        locus_fig.update_yaxes(
            title_text="-log<sub>10</sub>(P)",
            range=list(fig_bottom.layout.yaxis.range),
            zeroline=False,
            row=2,
            col=1,
            secondary_y=False
        )
        locus_fig.update_yaxes(
            title_text="Recombination rate",
            range=[0, recomb_max],
            autorange=False,
            zeroline=False,
            showgrid=False,
            row=2,
            col=1,
            secondary_y=True
        )

        locus_fig.update_yaxes(
            title_text="Genes",
            range=[-gene_track_rows + 0.5, 0.8],
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            row=3,
            col=1
        )

        locus_fig.update_xaxes(
            range=[gene_bp_start / 1e6, gene_bp_end / 1e6],
            showticklabels=False,
            zeroline=False,
            row=1,
            col=1
        )
        locus_fig.update_xaxes(
            range=[gene_bp_start / 1e6, gene_bp_end / 1e6],
            showticklabels=False,
            zeroline=False,
            row=2,
            col=1
        )
        locus_fig.update_xaxes(
            range=[gene_bp_start / 1e6, gene_bp_end / 1e6],
            title_text=format_chrom_axis_title(selected_chrom),
            zeroline=False,
            row=3,
            col=1
        )

        if hasattr(locus_fig.layout, "xaxis2"):
            locus_fig.layout.xaxis2.matches = "x"
        if hasattr(locus_fig.layout, "xaxis3"):
            locus_fig.layout.xaxis3.matches = "x"

        if hasattr(locus_fig.layout, "yaxis2"):
            locus_fig.layout.yaxis2.update(
                range=[0, recomb_max],
                autorange=False,
                tickmode="linear",
                dtick=20
            )
        if hasattr(locus_fig.layout, "yaxis4"):
            locus_fig.layout.yaxis4.update(
                range=[0, recomb_max],
                autorange=False,
                tickmode="linear",
                dtick=20
            )

        locus_fig.update_layout(
            height=combined_height,
            dragmode="zoom",
            showlegend=False,
            margin=dict(l=60, r=60, b=45, t=60),
        )
        locus_fig = apply_locusblend_plot_theme(locus_fig)

        _active_compare_signals = ["variant 1"]
        if active_mode != "Standard locus zoom":
            _active_compare_signals.append("variant 2")
        if active_mode == "Three-index LocusBlend":
            _active_compare_signals.append("variant 3")

        if locuscompare_mode == "Single blended compare plot":
            update_progress(progress_bar, status_box, 85, "Building single blended compare plot...")
            compare_fig, compare_blended_count = build_single_blended_locuscompare(
                df_top=df_top_plot,
                df_bottom=df_bottom_plot,
                idx1_ref=idx1_ref,
                idx2_ref=idx2_ref,
                idx3_ref=idx3_ref,
                idx1_label=idx1_label_for_plot,
                idx2_label=idx2_label_for_plot,
                idx3_label=idx3_label_for_plot,
                title_top=title_top,
                title_bottom=title_bottom,
                compare_size=compare_size,
                ld_labels=ld_labels,
                locusblend_mode=active_mode,
            )
            compare_counts = None
        else:
            update_progress(progress_bar, status_box, 85, "Building compare plots...")
            compare_fig, compare_counts = build_compare_figure_triptych(
                df_top=df_top_plot,
                df_bottom=df_bottom_plot,
                idx1_label=idx1_label_for_plot,
                idx2_label=idx2_label_for_plot,
                idx3_label=idx3_label_for_plot,
                title_top=title_top,
                title_bottom=title_bottom,
                compare_size=compare_size,
                ld_labels=ld_labels,
                active_signals=_active_compare_signals,
            )
            compare_blended_count = None

        compare_fig = apply_locusblend_plot_theme(compare_fig)
        compare_fig = apply_locus_compare_safe_autoscale(compare_fig)

        update_progress(progress_bar, status_box, 95, "Rendering figures...")

        n_flip_top = int(df_top_plot["coding_flipped"].fillna(False).sum()) if "coding_flipped" in df_top_plot.columns else 0
        n_flip_bottom = int(df_bottom_plot["coding_flipped"].fillna(False).sum()) if "coding_flipped" in df_bottom_plot.columns else 0

        if locuscompare_mode == "Single blended compare plot":
            compare_summary = f"Compare points (blended): {compare_blended_count:,}"
        else:
            _compare_parts = [f"v1: {compare_counts.get('variant 1', 0):,}"]
            if active_mode != "Standard locus zoom":
                _compare_parts.append(f"v2: {compare_counts.get('variant 2', 0):,}")
            if active_mode == "Three-index LocusBlend":
                _compare_parts.append(f"v3: {compare_counts.get('variant 3', 0):,}")
            compare_summary = f"Compare points ({'/'.join(_compare_parts)})"

        summary_text = (
            f"Dataset 1 SNP count: {n_top:,} | "
            f"Dataset 2 SNP count: {n_bottom:,} | "
            f"Genes shown ({gene_display_mode}): {len(genes_df):,} | "
            f"{compare_summary} | "
            f"{ld_labels['summary_label']}: {len(ref_snps):,} | "
            f"Coding flipped (dataset 1/dataset 2): {n_flip_top:,}/{n_flip_bottom:,}"
        )

        for _export_key in [
            "last_export_bytes",
            "last_export_filename",
            "last_export_mime",
            "last_export_signature",
            "last_export_error",
        ]:
            st.session_state.pop(_export_key, None)

        st.session_state["last_export_context"] = {
            "mode": active_mode,
            "compare_mode": locuscompare_mode,
            "chromosome": selected_chrom,
            "center_bp": int(bp),
            "window_kb": int(window_kb),
            "title_top": str(title_top),
            "title_bottom": str(title_bottom),
            "highlight_genes": str(st.session_state.get("highlight_genes", "")),
            "ld_status_caption": str(ld_status_caption),
        }

        st.session_state["last_locus_fig"] = locus_fig
        st.session_state["last_compare_fig"] = compare_fig
        st.session_state["last_compare_mode"] = locuscompare_mode
        st.session_state["last_compare_size"] = int(compare_size)
        st.session_state["last_summary_text"] = summary_text
        st.session_state["last_summary_metrics"] = {
            "n_top": n_top,
            "n_bottom": n_bottom,
            "n_genes": len(genes_df),
            "gene_mode": gene_display_mode,
            "compare_summary": compare_summary,
            "n_ref_snps": len(ref_snps),
            "summary_label": ld_labels["summary_label"],
            "n_flip_top": n_flip_top,
            "n_flip_bottom": n_flip_bottom,
        }
        if index_selection_method != "Manual input":
            st.session_state["last_auto_index_table"] = auto_index_table
            st.session_state["last_auto_index_caption"] = auto_index_summary
        else:
            st.session_state["last_auto_index_table"] = None
            st.session_state["last_auto_index_caption"] = ""
        st.session_state["last_index_selection_method"] = index_selection_method
        st.session_state["last_ld_status_caption"] = ld_status_caption
        st.session_state["last_progress_success"] = "Plots updated successfully."

        update_progress(progress_bar, status_box, 100, "Done.")
        status_box.success(st.session_state["last_progress_success"])

    # Display block — always runs, reads from sticky session_state. On
    # widget-only reruns (radio, uploader, params edits) this redraws the
    # previously cached figures without re-running any heavy pipeline work.
    st.caption(st.session_state["last_ld_status_caption"])

    config = {
        "toImageButtonOptions": {
            "format": "png",  # or "svg"
            "filename": "locus_plot",
            "width": 1800,
            "height": 1200,
            "scale": 3,
        }
    }

    st.plotly_chart(st.session_state["last_locus_fig"], use_container_width=True, config=config)

    st.markdown("### Locus compare")
    _last_compare_mode = st.session_state["last_compare_mode"]
    _last_compare_size = int(st.session_state["last_compare_size"])
    _compare_fig = st.session_state["last_compare_fig"]
    compare_config = get_locus_compare_plotly_config(config)

    if _last_compare_mode == "Single blended compare plot":
        _plot_size = _last_compare_size
        _compare_fig_display = clone_plotly_figure(_compare_fig)
        _compare_fig_display.update_layout(
            width=_plot_size,
            height=_plot_size,
            autosize=False,
            margin=dict(l=70, r=40, t=50, b=70),
        )
        blended_config = {
            "responsive": False,
            "toImageButtonOptions": {
                "format": "png",
                "filename": "locus_compare_blended",
                "width": _plot_size,
                "height": _plot_size,
                "scale": 3,
            },
        }
        blended_config = get_locus_compare_plotly_config(blended_config)
        _left, _center, _right = st.columns([1, 2, 1])
        with _center:
            st.plotly_chart(
                _compare_fig_display,
                use_container_width=False,
                config=blended_config,
            )
    else:
        st.plotly_chart(_compare_fig, use_container_width=True, config=compare_config)

    # Auto index table
    if st.session_state.get("last_auto_index_table") is not None:
        _auto_tbl = st.session_state["last_auto_index_table"]
        _auto_cap = st.session_state.get("last_auto_index_caption", "")
        if _auto_cap:
            st.caption(_auto_cap)
        st.dataframe(
            _auto_tbl.rename(columns={
                "rank": "Rank",
                "DISPLAY_ID": "Variant",
                "REF_SNP": "REF_SNP",
                "CHR": "CHR",
                "BP": "BP",
                "P": "P",
                "_source": "Source",
            }),
            use_container_width=True,
            hide_index=True,
        )

    # Summary display
    if "last_summary_metrics" in st.session_state:
        m = st.session_state["last_summary_metrics"]
        _index_method = st.session_state.get("last_index_selection_method", "Manual input")
        with st.container(border=True):
            c1, c2, c3, c4 = st.columns(4)
            c1.markdown(f"**Dataset 1 SNPs**<br><span style='font-size:1.3em'>{m['n_top']:,}</span>", unsafe_allow_html=True)
            c2.markdown(f"**Dataset 2 SNPs**<br><span style='font-size:1.3em'>{m['n_bottom']:,}</span>", unsafe_allow_html=True)
            c3.markdown(f"**Genes** ({m['gene_mode']})<br><span style='font-size:1.3em'>{m['n_genes']:,}</span>", unsafe_allow_html=True)
            c4.markdown(f"**{m['summary_label']}**<br><span style='font-size:1.3em'>{m['n_ref_snps']:,}</span>", unsafe_allow_html=True)
            st.caption(
                f"Index selection: {_index_method}  |  "
                f"{m['compare_summary']}  |  "
                f"Coding flipped (dataset 1/dataset 2): {m['n_flip_top']:,}/{m['n_flip_bottom']:,}"
            )
    elif "last_summary_text" in st.session_state:
        st.write(st.session_state["last_summary_text"])

    render_export_controls()

    if not should_compute:
        st.info("Using previously rendered plots. Select Update plot to apply new settings.")

except Exception as e:
    _progress_bar = locals().get("progress_bar", None)
    _status_box = locals().get("status_box", None)

    if _progress_bar is not None:
        try:
            _progress_bar.empty()
        except Exception:
            pass

    if _status_box is not None:
        try:
            _status_box.empty()
        except Exception:
            pass

    log(traceback.format_exc())
    log(f"APP FAILED: {repr(e)}")
    st.exception(e)
