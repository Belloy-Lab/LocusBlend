"""Chromosome normalization and reference-resource path resolution.

Helpers extracted from app.py. They are import-safe and depend only on the
static values in locusblend_web.config, not on Streamlit or runtime data.
"""

import glob
import os
import shutil
import time

import pandas as pd

from locusblend_web.config import (
    BIN_DIR,
    DATA_DIR,
    INTERNAL_1000G_ANCESTRIES,
    INTERNAL_1000G_DEFAULT_ANCESTRY,
    PROJECT_ROOT,
)


def log(msg):
    """Print a timestamped message; mirrors app.py's log() to avoid a circular import."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def normalize_chrom(chrom):
    """Return a canonical chromosome string for autosomes and chromosome X."""
    if chrom is None:
        return ""
    try:
        if pd.isna(chrom):
            return ""
    except (TypeError, ValueError):
        pass

    s = str(chrom).strip()
    if not s:
        return ""
    if s.lower().startswith("chr"):
        s = s[3:].strip()
    if s.endswith(".0"):
        s = s[:-2]

    if s.upper() == "X" or s == "23":
        return "X"
    if s.isdigit():
        n = int(s)
        if 1 <= n <= 22:
            return str(n)
    return s


def normalize_internal_1000g_ancestry(value):
    """Return a supported 1000G super-population code, defaulting to EUR."""
    ancestry = str(value or INTERNAL_1000G_DEFAULT_ANCESTRY).strip().upper()
    return ancestry if ancestry in INTERNAL_1000G_ANCESTRIES else INTERNAL_1000G_DEFAULT_ANCESTRY


def get_internal_1000g_prefix(chrom, ancestry=INTERNAL_1000G_DEFAULT_ANCESTRY):
    ancestry = normalize_internal_1000g_ancestry(ancestry)
    chrom = normalize_chrom(chrom)
    return DATA_DIR / f"1000g_{ancestry}_hg38_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_include_MHC_ch{chrom}"


def get_gtf_path_for_chrom(chrom):
    """Build the chromosome-specific GENCODE GTF path and verify it exists."""
    chrom = normalize_chrom(chrom)
    path = DATA_DIR / f"gencode.v49.annotation.chr{chrom}.gtf.gz"
    if chrom == "X" and not path.is_file():
        full_path = DATA_DIR / "gencode.v49.annotation.gtf.gz"
        if full_path.is_file():
            return str(full_path)
        raise FileNotFoundError(
            f"Missing GENCODE annotation for chromosome X. Add "
            f"{DATA_DIR / 'gencode.v49.annotation.chrX.gtf.gz'} or "
            f"{DATA_DIR / 'gencode.v49.annotation.gtf.gz'}."
        )
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing GENCODE annotation for chromosome {chrom}: {path}"
        )
    return str(path)


def _normalize_bfile_prefix(p):
    p = str(p).strip()
    for suf in [".bed", ".bim", ".fam"]:
        if p.endswith(suf):
            p = p[:-len(suf)]
    return p


def _resolve_bfile_prefix(bfile_prefix):
    raw = str(bfile_prefix).strip()
    prefix = _normalize_bfile_prefix(raw)

    search_dirs = [DATA_DIR, PROJECT_ROOT]

    candidates = []
    if prefix:
        candidates.append(prefix)
        for d in search_dirs:
            candidates.append(str(d / os.path.basename(prefix)))

    for d in search_dirs:
        for bim in sorted(glob.glob(str(d / "*.bim"))):
            candidates.append(bim[:-4])

    seen = set()
    uniq = []
    for c in candidates:
        if c not in seen:
            uniq.append(c)
            seen.add(c)

    for c in uniq:
        bed = c + ".bed"
        bim = c + ".bim"
        fam = c + ".fam"
        if os.path.exists(bed) and os.path.exists(bim) and os.path.exists(fam):
            log(f"Resolved PLINK prefix: {c}")
            return c

    debug_files = []
    for d in search_dirs:
        debug_files.extend(sorted(glob.glob(str(d / "*"))))
    raise FileNotFoundError(
        "Could not resolve a valid PLINK bfile prefix. "
        f"Input was: {raw!r}. "
        f"Searched directories: {[str(d) for d in search_dirs]}. "
        f"Available files include: {[os.path.basename(x) for x in debug_files[:50]]}"
    )


def _find_plink_exec(plink_path=str(BIN_DIR / "plink")):
    candidates = [str(plink_path).strip(), str(BIN_DIR / "plink"), "./plink", "plink"]
    for c in candidates:
        if c and (os.path.exists(c) or shutil.which(c)):
            return c
    raise FileNotFoundError("PLINK executable not found. Please install plink first.")
