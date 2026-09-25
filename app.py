from io import BytesIO
import re
import traceback

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

from locusblend_web.config import (
    BIN_DIR,
    COLOR_MAPPING,
    DATA_DIR,
    INTERNAL_1000G_ANCESTRY_OPTIONS,
    INTERNAL_1000G_DEFAULT_ANCESTRY,
)

from locusblend_web.references import (
    chrom_sort_key,
    format_chrom_label,
    format_internal_1000g_ancestry_option,
    get_internal_bfile_prefix_for_chrom,
    get_gtf_path_for_chrom,
    get_internal_1000g_prefix,
    get_supported_chromosomes,
    is_supported_chrom,
    normalize_chrom,
    normalize_internal_1000g_ancestry,
)

from locusblend_web.variants import (
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
    get_ld_reference_labels,
    merge_ld_annot,
    parse_uploaded_ld_long,
    parse_uploaded_ld_matrix,
    prepare_clump_candidates,
    read_reference_bim,
    run_plink_clump_for_auto_indices,
)

from locusblend_web.logutil import log

from locusblend_web.plotting import (
    add_gene_track_to_subplot,
    apply_locus_compare_safe_autoscale,
    apply_locusblend_plot_theme,
    build_compare_figure_triptych,
    build_single_blended_locuscompare,
    format_chrom_axis_title,
    get_locus_compare_plotly_config,
    get_plotly_locus_py,
)

from locusblend_web.export import (
    build_locusblend_export_image,
    clone_plotly_figure,
)

from locusblend_web.styles import inject_locusblend_css

from locusblend_web.ui import (
    load_locusblend_page_icon,
    render_documentation_page,
    render_locusblend_welcome_card,
    render_page_switcher,
    render_sidebar_app_header,
    render_sidebar_visual_strip,
    render_washu_header,
    set_streamlit_chrome_minimal,
    update_progress,
)

from locusblend_web.uploads import (
    format_upload_sync_message,
    infer_uploaded_locus_sync,
    make_dataset_title_from_filename,
    make_uploaded_dataset_signature,
    read_locus_csv,
    read_locus_csv_uploaded,
    recommend_y_axis_max_for_dataset,
)


@st.cache_data(show_spinner=False)
def load_locus_csv(path):
    return read_locus_csv(path)


@st.cache_data(show_spinner=False)
def load_locus_csv_uploaded(file_bytes, file_name):
    return read_locus_csv_uploaded(file_bytes, file_name)


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


@st.cache_data(show_spinner=False)
def read_uploaded_ld_long(file_bytes, file_name):
    return parse_uploaded_ld_long(file_bytes, file_name)


@st.cache_data(show_spinner=False)
def read_uploaded_ld_matrix(file_bytes, file_name):
    return parse_uploaded_ld_matrix(file_bytes, file_name)


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
