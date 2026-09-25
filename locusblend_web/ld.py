"""Uploaded-LD parsing and pure LD-map computation.

Helpers extracted from app.py. Nothing in this module invokes PLINK or
touches Streamlit, session state, the filesystem, or the UI.
"""

from io import BytesIO

import numpy as np
import pandas as pd

from locusblend_web.variants import dedup_columns


def parse_uploaded_ld_long(file_bytes, file_name):
    """
    Parse an uploaded LD long-format table. Returns:
      ld_long:     DataFrame with canonical columns SNP_A, SNP_B, R2 (R2 in [0,1])
      ld_universe: tuple of SNP IDs known to the LD reference, taken from
                   self-pair rows where SNP_A == SNP_B and R2 == 1.
    Raises ValueError if required columns or self-pair rows are missing.
    """
    raw = BytesIO(file_bytes)
    name = (file_name or "").lower()
    if name.endswith((".tsv", ".tab", ".txt", ".ld")):
        df = pd.read_csv(raw, sep=r"\s+", engine="python")
    else:
        try:
            df = pd.read_csv(raw, sep=None, engine="python")
        except Exception:
            raw.seek(0)
            df = pd.read_csv(raw)
    df = dedup_columns(df)

    alias = {
        "snp_a": "SNP_A", "snp1": "SNP_A", "id1": "SNP_A", "snpa": "SNP_A",
        "snp_b": "SNP_B", "snp2": "SNP_B", "id2": "SNP_B", "snpb": "SNP_B",
        "r2": "R2", "rsq": "R2", "r^2": "R2",
    }
    rename = {c: alias[c.lower()] for c in df.columns if c.lower() in alias}
    df = df.rename(columns=rename)

    missing = [c for c in ["SNP_A", "SNP_B", "R2"] if c not in df.columns]
    if missing:
        raise ValueError(
            f"Uploaded LD long table missing required columns: {missing}. "
            "Expected SNP_A, SNP_B, R2 (or aliases SNP1/ID1/SNP2/ID2/r2)."
        )

    df = df[["SNP_A", "SNP_B", "R2"]].copy()
    df["SNP_A"] = df["SNP_A"].astype(str).str.strip()
    df["SNP_B"] = df["SNP_B"].astype(str).str.strip()
    df["R2"] = pd.to_numeric(df["R2"], errors="coerce")
    df = df.dropna(subset=["SNP_A", "SNP_B", "R2"])
    df = df[(df["SNP_A"] != "") & (df["SNP_B"] != "")]
    df["R2"] = df["R2"].clip(lower=0.0, upper=1.0)

    self_pairs = df[(df["SNP_A"] == df["SNP_B"]) & np.isclose(df["R2"], 1.0)]
    if self_pairs.empty:
        raise ValueError(
            "Uploaded LD long table must include self-pair rows "
            "SNP_A=SNP_B, R2=1 to distinguish low LD from missing LD."
        )

    ld_universe = tuple(sorted(set(self_pairs["SNP_A"].astype(str))))
    return df, ld_universe


def parse_uploaded_ld_matrix(file_bytes, file_name):
    """
    Parse an uploaded LD matrix. First column = row SNP IDs, first row =
    column SNP IDs, cells = R^2. Returns the matrix as a DataFrame with
    SNP IDs on both axes, plus a tuple ld_universe = sorted(rows | cols).
    """
    raw = BytesIO(file_bytes)
    name = (file_name or "").lower()
    if name.endswith((".tsv", ".tab", ".txt", ".ld")):
        df = pd.read_csv(raw, sep=r"\s+", engine="python", index_col=0)
    else:
        try:
            df = pd.read_csv(raw, sep=None, engine="python", index_col=0)
        except Exception:
            raw.seek(0)
            df = pd.read_csv(raw, index_col=0)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.apply(pd.to_numeric, errors="coerce")

    row_ids = [s for s in df.index.tolist() if s and s.lower() != "nan"]
    col_ids = [s for s in df.columns.tolist() if s and s.lower() != "nan"]
    if not row_ids and not col_ids:
        raise ValueError("Uploaded LD matrix has no SNP IDs in row or column names.")

    ld_universe = tuple(sorted(set(row_ids) | set(col_ids)))
    return df, ld_universe


def _safe_idx_str(x):
    if x is None:
        return None
    try:
        if pd.isna(x):
            return None
    except (TypeError, ValueError):
        pass
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return None
    return s


def compute_ld_maps_from_uploaded_long(ld_long, window_snps, idx1_ref, idx2_ref, idx3_ref):
    """
    Build ld_maps from a standardized long-format LD table. Output
    structure matches compute_ld_maps_with_plink so the downstream
    ld_annot / merge_ld_annot pipeline is unchanged.
    """
    df = ld_long[["SNP_A", "SNP_B", "R2"]].copy()
    df["SNP_A"] = df["SNP_A"].astype(str)
    df["SNP_B"] = df["SNP_B"].astype(str)

    self_pairs = df[(df["SNP_A"] == df["SNP_B"]) & np.isclose(df["R2"], 1.0)]
    ld_universe = set(self_pairs["SNP_A"].astype(str))

    window_set = {str(x) for x in (window_snps or ()) if pd.notna(x)}
    ref_snps = tuple(sorted(ld_universe & window_set)) if window_set else tuple(sorted(ld_universe))

    def in_universe(idx):
        s = _safe_idx_str(idx)
        return s is not None and s in ld_universe

    index_status = {
        "variant 1": in_universe(idx1_ref),
        "variant 2": in_universe(idx2_ref),
        "variant 3": in_universe(idx3_ref),
    }

    def build_map(idx_ref):
        s = _safe_idx_str(idx_ref)
        if s is None or s not in ld_universe:
            return {}
        a = df[df["SNP_A"] == s][["SNP_B", "R2"]].rename(columns={"SNP_B": "OTHER"})
        b = df[df["SNP_B"] == s][["SNP_A", "R2"]].rename(columns={"SNP_A": "OTHER"})
        m = pd.concat([a, b], axis=0, ignore_index=True).dropna(subset=["OTHER", "R2"])
        m["OTHER"] = m["OTHER"].astype(str)
        m = m.groupby("OTHER", as_index=False)["R2"].max()
        out = dict(zip(m["OTHER"], m["R2"].astype(float)))
        out[s] = 1.0
        return out

    ld_maps = {
        "r2_1": build_map(idx1_ref),
        "r2_2": build_map(idx2_ref),
        "r2_3": build_map(idx3_ref),
    }
    return ld_maps, index_status, ref_snps


def compute_ld_maps_from_uploaded_matrix(ld_matrix, window_snps, idx1_ref, idx2_ref, idx3_ref):
    """
    Build ld_maps from a square-ish LD matrix (rows/cols = SNP IDs,
    cells = R^2). Output structure matches compute_ld_maps_with_plink.
    """
    rows = {str(x) for x in ld_matrix.index.tolist()}
    cols = {str(x) for x in ld_matrix.columns.tolist()}
    ld_universe = rows | cols

    window_set = {str(x) for x in (window_snps or ()) if pd.notna(x)}
    ref_snps = tuple(sorted(ld_universe & window_set)) if window_set else tuple(sorted(ld_universe))

    def in_universe(idx):
        s = _safe_idx_str(idx)
        return s is not None and s in ld_universe

    index_status = {
        "variant 1": in_universe(idx1_ref),
        "variant 2": in_universe(idx2_ref),
        "variant 3": in_universe(idx3_ref),
    }

    def extract_vec(idx_ref):
        s = _safe_idx_str(idx_ref)
        if s is None:
            return {}
        if s in rows:
            vec = ld_matrix.loc[s]
        elif s in cols:
            vec = ld_matrix[s]
        else:
            return {}
        if isinstance(vec, pd.DataFrame):
            vec = vec.iloc[0]
        vec = pd.to_numeric(vec, errors="coerce").dropna()
        vec.index = vec.index.astype(str)
        out = {k: float(v) for k, v in vec.items()}
        out[s] = 1.0
        return out

    ld_maps = {
        "r2_1": extract_vec(idx1_ref),
        "r2_2": extract_vec(idx2_ref),
        "r2_3": extract_vec(idx3_ref),
    }
    return ld_maps, index_status, ref_snps
