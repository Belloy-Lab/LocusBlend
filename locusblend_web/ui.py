"""Static Streamlit UI/rendering helpers for the LocusBlend web app.

Rendering-only helpers extracted from app.py: page chrome, sidebar header,
asset data URIs, the welcome card, the documentation view, the page switcher,
and the page icon. They create Streamlit output but never touch session state.
"""

import base64
import html as html_lib
import time
from mimetypes import guess_type

import streamlit as st
import streamlit.components.v1 as components

from locusblend_web.config import ASSET_DIR


def log(msg):
    """Print a timestamped message; mirrors app.py's log() to avoid a circular import."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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
