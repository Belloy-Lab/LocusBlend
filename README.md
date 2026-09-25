# LocusBlend

LocusBlend is a browser-based visualization tool for comparing and integrating
locus-level genetic association signals, LD information, gene annotations, and
related genomic context.

## Web interface and Python package

This repository contains the LocusBlend Streamlit web application. The
Streamlit layer lives in `app.py`, and the reusable logic lives in the
`locusblend_web` package.

LocusBlend can also be used through a separate Python package. The web
interface and the Python package provide complementary ways to generate
LocusBlend visualizations.

## Features

- One or two summary-statistic datasets shown as stacked locus panels
  (dataset 1 / top and dataset 2 / bottom).
- Three visualization modes: Standard locus zoom (one index variant),
  Two-index LocusBlend, and Three-index LocusBlend (blended LD coloring).
- Locus compare views: three separate compare plots, or a single blended plot.
- Manual index-variant entry or automatic selection by LD clumping from either
  dataset.
- LD from internal 1000 Genomes reference panels, from an uploaded LD long
  table, or from an uploaded LD matrix.
- Internal 1000G reference support for AFR, AMR, EAS, EUR, and SAS, with EUR as
  the default ancestry.
- GRCh38 / hg38 coordinates only; chromosomes 1-22 and X.
- GENCODE gene annotation track, with `protein_coding` or all-gene display and
  gene highlighting by name.
- Recombination-rate track from a local BigWig file.
- Export of the current figure to PNG or PDF, with US Letter or custom page
  size, configurable DPI, and optional legend images.

## Repository structure

```text
app.py                  Streamlit composition root: bootstrap, widgets, session
                        state, compute orchestration, display
locusblend_web/
    config.py           Static constants, paths, 1000G ancestries, color table
    references.py       Chromosome normalization, reference paths, PLINK lookup
    variants.py         Column canonicalization, UID and allele normalization
    genes.py            GENCODE loading, gene filtering, gene track layout
    ld.py               Uploaded-LD parsing, PLINK LD and clumping, BIM loading
    plotting.py         Plotly figure construction and figure theming
    export.py           Kaleido/Pillow rendering, PNG and PDF composition
    uploads.py          Summary-statistic parsing, dataset metadata, sync hints
    ui.py               Static Streamlit UI helpers (headers, docs, page icon)
    styles.py           Application stylesheet
    logutil.py          Shared timestamped stdout logger
assets/                 Logos, legend and diagram images used by the app
.streamlit/config.toml  Streamlit theme and toolbar configuration
requirements.txt        Pinned Python dependencies
```

## Installation

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate
pip install -r requirements.txt
```

## Running locally

Provision the reference data described below, then start the application from
the repository root:

```bash
streamlit run app.py
```

Open the local URL printed by Streamlit in a browser.

## Reference data

Large runtime reference resources are intentionally not stored in Git. The
`data/` and `bin/` directories are excluded from version control and must be
provisioned locally next to `app.py`, relative to the repository root. No
server-specific paths are required.

- Internal 1000 Genomes PLINK binary reference panels, one set per ancestry and
  chromosome, using the prefix
  `1000g_{ANCESTRY}_hg38_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_include_MHC_ch{CHROM}`
  with matching `.bed`, `.bim`, and `.fam` files, where `ANCESTRY` is one of
  `AFR`, `AMR`, `EAS`, `EUR`, `SAS` and `CHROM` is `1`-`22` or `X`.
- GENCODE gene annotation, one gzip-compressed GTF per chromosome, named
  `gencode.v49.annotation.chr{CHROM}.gtf.gz` (chromosome X may fall back to a
  whole-genome `gencode.v49.annotation.gtf.gz`).
- Recombination-rate BigWig file named `recomb1000GAvg.bw`.
- PLINK executable at `bin/plink` (fallbacks: `./plink`, `plink` on `PATH`).

When no summary-statistic files are uploaded, the application loads two example
datasets from `data/`: `merged_female_withld_PSEN1.csv` and
`Nedelec_monocytes_PSEN1_example.csv` (provisioned locally, not stored in Git).

## Input data

Supported summary-statistic upload formats are `.csv`, `.tsv`, `.txt`,
`.csv.gz`, `.tsv.gz`, and `.txt.gz`; gzip-compressed files are read directly and
comma, tab, and whitespace delimiters are supported.

Required fields are chromosome, base-pair position, p-value, and both alleles.
Recognized names include `CHR` (`chr`, `#chr`, `chrom`, `chromosome`), `BP`
(`bp`, `pos`, `position`, `base_pair_location`), `P` (`pval`, `pvalue`,
`p_value`, `P-value`), and `A1` / `A2` (`EA` / `NEA`, `effect_allele` /
`other_allele`, `ALLELE1` / `ALLELE0`), plus optional `rsid` (`SNP`,
`MarkerName`, `ID`), `BETA`, `SE`, `A1FREQ`, and `N`.

Coordinates must use GRCh38 / hg38. Chromosome X may be provided as `X`,
`chrX`, `23`, or `chr23`.

## LD sources

1. Internal 1000 Genomes reference. LD is computed with the local PLINK
   executable against the chromosome-specific panel for the selected ancestry.
2. Uploaded LD long table. Pairwise table with columns `SNP_A`, `SNP_B`, and
   `R2` (aliases such as `SNP1`/`ID1`, `SNP2`/`ID2`, and `r2` are accepted).
   Self-pair rows (`SNP_A` equal to `SNP_B` with `R2` = 1) are required so that
   low LD can be distinguished from missing LD.
3. Uploaded LD matrix. Square R2 matrix with SNP IDs as the row and column
   names.

## Development

- `app.py` contains Streamlit-specific orchestration: page bootstrap, widgets,
  session state, compute flow, and display flow.
- The `locusblend_web` package contains the reusable backend logic for
  configuration, reference handling, variant normalization, gene annotation,
  LD, plotting, export, and uploads, together with the presentation modules
  `ui.py` and `styles.py`.
- Most backend modules do not import Streamlit and can be imported and exercised
  independently; only `app.py`, `ui.py`, and `styles.py` depend on it.
- Each package module starts with a short docstring describing its
  responsibility. This structure keeps the Streamlit layer thin and makes the
  data, LD, plotting, and export logic easier to test and maintain.

## Notes and limitations

- GRCh38 / hg38 coordinates only; no liftover is performed. GRCh37 / hg19 or
  hg18 inputs must be lifted over before upload.
- Reference resources are not distributed with the repository and must be
  provisioned locally. Missing files are reported at runtime with the expected
  paths. Internal-reference LD and automatic clumping require local PLINK.
- Results depend on the reference data used: for the internal reference, on the
  variant coverage of the selected ancestry panel; for uploaded LD, on the
  uploaded table or matrix.
- Static PNG and PDF export requires the Plotly image export dependencies
  (`kaleido`, and `pillow` for composition). Kaleido v1 may additionally
  require server-side Chrome or Chromium.

## Citation

If you use LocusBlend in your work, please cite the associated LocusBlend
publication. Citation details will be updated here upon publication.

## Contact

LocusBlend is developed at Washington University in St. Louis. For questions,
bug reports, or feature requests, please use the repository's issue tracker.
