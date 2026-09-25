"""LD parsing, PLINK/internal-1000G LD computation, and clumping helpers.

Helpers extracted from app.py. They may run PLINK subprocesses and use
temporary files, but they never touch Streamlit, session state, or the UI.
"""

import os
import subprocess
import tempfile
import time
from io import BytesIO

import numpy as np
import pandas as pd

from locusblend_web.config import BIN_DIR
from locusblend_web.references import (
    _find_plink_exec,
    _resolve_bfile_prefix,
    normalize_chrom,
)
from locusblend_web.variants import chrom_mask, dedup_columns


def log(msg):
    """Print a timestamped message; mirrors app.py's log() to avoid a circular import."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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


def prepare_clump_candidates(df_ref, selected_chrom, center_bp, window_bp):
    """Filter to selected chromosome and BP window, keep rows with valid P
    and REF_SNP, sort by P ascending.  Returns columns: REF_SNP, DISPLAY_ID,
    CHR, BP, P."""
    df = dedup_columns(df_ref).copy()
    mask = (
        chrom_mask(df, selected_chrom)
        & (df["BP"] >= center_bp - window_bp)
        & (df["BP"] <= center_bp + window_bp)
    )
    cand = df.loc[mask].copy()
    cand["P"] = pd.to_numeric(cand["P"], errors="coerce")
    cand = cand[
        cand["P"].notna() & (cand["P"] > 0) & (cand["P"] <= 1)
    ].copy()
    if "REF_SNP" not in cand.columns:
        raise KeyError("Candidate table missing REF_SNP column — run reference matching first.")
    cand = cand[cand["REF_SNP"].notna()].copy()
    cand = cand.drop_duplicates("REF_SNP").copy()
    cand = cand.sort_values("P").reset_index(drop=True)
    if "DISPLAY_ID" not in cand.columns:
        cand["DISPLAY_ID"] = cand["REF_SNP"]
    return cand[["REF_SNP", "DISPLAY_ID", "CHR", "BP", "P"]]


def run_plink_clump_for_auto_indices(
    candidates,
    bfile_prefix,
    selected_chrom,
    start_bp,
    end_bp,
    max_indices,
    clump_r2=0.01,
    plink_path=None,
):
    """Run PLINK --clump on the candidate list, returning up to max_indices
    independent lead SNPs (ranked by PLINK clump order)."""
    if plink_path is None:
        plink_path = _find_plink_exec()
    bfile_prefix = _resolve_bfile_prefix(bfile_prefix)
    selected_chrom = normalize_chrom(selected_chrom)

    if candidates.empty:
        raise ValueError("No valid candidate SNPs for clumping.")

    window_kb = max(1, int((int(end_bp) - int(start_bp)) / 1000) + 1)

    with tempfile.TemporaryDirectory(prefix="plink_clump_") as tmpdir:
        clump_in = os.path.join(tmpdir, "clump_input.txt")
        extract_in = os.path.join(tmpdir, "extract.snplist")
        out_prefix = os.path.join(tmpdir, "clump_out")

        # Write clump input (SNP P)
        with open(clump_in, "w") as f:
            f.write("SNP P\n")
            for _, row in candidates.iterrows():
                f.write(f"{row['REF_SNP']} {row['P']:.6e}\n")

        # Write extract list (candidate REF_SNPs only)
        with open(extract_in, "w") as f:
            for snp in candidates["REF_SNP"]:
                f.write(f"{snp}\n")

        cmd = [
            plink_path,
            "--threads", "1",
            "--bfile", bfile_prefix,
            "--chr", selected_chrom,
            "--from-bp", str(int(start_bp)),
            "--to-bp", str(int(end_bp)),
            "--extract", extract_in,
            "--clump", clump_in,
            "--clump-snp-field", "SNP",
            "--clump-field", "P",
            "--clump-p1", "1",
            "--clump-p2", "1",
            "--clump-r2", str(clump_r2),
            "--clump-kb", str(window_kb),
            "--out", out_prefix,
        ]
        if selected_chrom == "X":
            cmd.extend(["--allow-extra-chr"])

        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        clumped_path = out_prefix + ".clumped"
        if res.returncode != 0:
            log(f"PLINK clump failed:\n{res.stdout}")
            if not os.path.exists(clumped_path):
                raise ValueError(f"PLINK clumping failed. Output:\n{res.stdout}")

        if not os.path.exists(clumped_path):
            raise ValueError("PLINK clumping produced no .clumped output.")

        clumped = pd.read_csv(clumped_path, sep=r"\s+")
        if clumped.empty or "SNP" not in clumped.columns:
            raise ValueError("PLINK clumping produced an empty result — no independent lead SNPs found.")

        lead_snps = clumped["SNP"].astype(str).str.strip().head(max_indices).tolist()

    # Map back to candidate metadata
    selected_rows = []
    for snp in lead_snps:
        match = candidates[candidates["REF_SNP"].astype(str) == snp]
        if not match.empty:
            selected_rows.append(match.iloc[0].to_dict())

    if not selected_rows:
        raise ValueError("None of the PLINK lead SNPs matched the candidate table.")

    out = pd.DataFrame(selected_rows)
    out["rank"] = range(1, len(out) + 1)
    return out[["rank", "DISPLAY_ID", "REF_SNP", "CHR", "BP", "P"]]


def greedy_clump_uploaded_ld_for_auto_indices(
    candidates,
    max_indices=3,
    clump_r2=0.01,
    ld_long_table=None,
    ld_matrix_table=None,
):
    """Greedy independent-variant selection using uploaded LD.

    Sorts candidates by P ascending, iteratively picks the top-significant
    SNP not yet blocked, then blocks all SNPs with r² > clump_r2 to the
    lead SNP.
    """
    ld_long = ld_long_table
    ld_matrix = ld_matrix_table

    if ld_matrix is not None:
        rows_set = {str(x) for x in ld_matrix.index}
        cols_set = {str(x) for x in ld_matrix.columns}
        ld_universe = rows_set | cols_set

        def get_partners(ref_snp):
            partners = {}
            s = str(ref_snp)
            if s in rows_set:
                vec = ld_matrix.loc[s]
                if isinstance(vec, pd.DataFrame):
                    vec = vec.iloc[0]
                vec = pd.to_numeric(vec, errors="coerce").dropna()
                partners.update({str(k): float(v) for k, v in vec.items()})
            if s in cols_set:
                # column lookup
                vec = ld_matrix[s]
                if isinstance(vec, pd.DataFrame):
                    vec = vec.iloc[:, 0]
                vec = pd.to_numeric(vec, errors="coerce").dropna()
                vec.index = ld_matrix.index.astype(str)
                for k, v in vec.items():
                    if str(k) not in partners:
                        partners[str(k)] = float(v)
            return partners
    elif ld_long is not None:
        ld_universe_set = set()
        ld_map = {}
        for _, row in ld_long.iterrows():
            a, b, r = str(row["SNP_A"]), str(row["SNP_B"]), float(row["R2"])
            ld_universe_set.update([a, b])
            ld_map.setdefault(a, {})[b] = r
            ld_map.setdefault(b, {})[a] = r
        ld_universe = ld_universe_set

        def get_partners(ref_snp):
            return ld_map.get(str(ref_snp), {})
    else:
        raise ValueError("No uploaded LD table provided for greedy clumping.")

    # Filter to candidates present in LD universe
    cand = candidates[candidates["REF_SNP"].astype(str).isin(ld_universe)].copy()
    if cand.empty:
        raise ValueError("No candidate SNPs found in the uploaded LD reference.")

    blocked = set()
    selected = []

    for _, row in cand.iterrows():
        s = str(row["REF_SNP"])
        if s in blocked:
            continue
        selected.append(row.to_dict())
        if len(selected) >= max_indices:
            break
        partners = get_partners(s)
        for partner, r in partners.items():
            if r > clump_r2:
                blocked.add(partner)

    if not selected:
        raise ValueError(
            f"No independent lead SNPs found at r²={clump_r2} "
            f"within the selected locus window."
        )

    out = pd.DataFrame(selected)
    out["rank"] = range(1, len(out) + 1)
    needed_cols = ["rank", "DISPLAY_ID", "REF_SNP", "CHR", "BP", "P"]
    for c in needed_cols:
        if c not in out.columns:
            out[c] = pd.NA
    return out[needed_cols]


def read_reference_bim(bfile_prefix):
    bfile_prefix = _resolve_bfile_prefix(bfile_prefix)
    bim_path = f"{bfile_prefix}.bim"

    log(f"Using BIM path: {bim_path}")
    log(f"cwd: {os.getcwd()}")

    bim = pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"]
    )
    bim["CHR"] = bim["CHR"].map(normalize_chrom)
    bim["BP"] = pd.to_numeric(bim["BP"], errors="coerce")
    bim["A1"] = bim["A1"].astype(str).str.upper()
    bim["A2"] = bim["A2"].astype(str).str.upper()
    bim["SNP"] = bim["SNP"].astype(str)
    bim = bim.dropna(subset=["BP"]).copy()
    return bim


def _read_plink_ld_table(ld_path, index_snp):
    if not os.path.exists(ld_path):
        return {}

    ld = pd.read_csv(ld_path, sep=r"\s+")
    if ld.empty or "SNP_A" not in ld.columns or "SNP_B" not in ld.columns or "R2" not in ld.columns:
        return {}

    ld = ld.copy()
    ld["SNP_A"] = ld["SNP_A"].astype(str)
    ld["SNP_B"] = ld["SNP_B"].astype(str)
    ld["R2"] = pd.to_numeric(ld["R2"], errors="coerce")

    sub = ld[ld["SNP_A"] == str(index_snp)][["SNP_B", "R2"]].dropna()
    out = dict(zip(sub["SNP_B"], sub["R2"]))
    out[str(index_snp)] = 1.0
    return out


def build_ld_maps_with_plink(
    bfile_prefix,
    chrom,
    start,
    end,
    window_snps,
    idx1_ref,
    idx2_ref,
    idx3_ref,
    plink_path=str(BIN_DIR / "plink"),
    bim_loader=None,
):
    """Run PLINK --r2 for each index variant and build ld_maps.

    bim_loader lets the caller inject the app's cached BIM reader; when it is
    omitted, the uncached read_reference_bim is used.
    """
    plink_exec = _find_plink_exec(plink_path)
    bfile_prefix = _resolve_bfile_prefix(bfile_prefix)
    chrom = normalize_chrom(chrom)

    load_bim = bim_loader or read_reference_bim
    bim = load_bim(bfile_prefix)
    ref_window = bim[
        (bim["CHR"] == chrom) &
        (bim["BP"] >= int(start)) &
        (bim["BP"] <= int(end))
    ].copy()

    ref_snps = set(ref_window["SNP"].tolist())
    extract_snps = sorted(set([str(x) for x in window_snps if pd.notna(x) and str(x) in ref_snps]))

    index_status = {
        "variant 1": idx1_ref is not None and str(idx1_ref) in ref_snps,
        "variant 2": idx2_ref is not None and str(idx2_ref) in ref_snps,
        "variant 3": idx3_ref is not None and str(idx3_ref) in ref_snps,
    }

    ld_maps = {"r2_1": {}, "r2_2": {}, "r2_3": {}}

    if len(extract_snps) == 0:
        return ld_maps, index_status, tuple()

    kb_span = max(1000, int((int(end) - int(start)) / 1000) + 100)

    with tempfile.TemporaryDirectory(prefix="plink_ld_") as tmpdir:
        extract_path = os.path.join(tmpdir, "extract.snplist")
        with open(extract_path, "w") as f:
            for snp in extract_snps:
                f.write(f"{snp}\n")

        index_inputs = [
            ("r2_1", idx1_ref),
            ("r2_2", idx2_ref),
            ("r2_3", idx3_ref),
        ]

        for r2_col, index_snp in index_inputs:
            if index_snp is None or str(index_snp) not in ref_snps:
                continue

            out_prefix = os.path.join(tmpdir, f"ld_{r2_col}")
            cmd = [
                plink_exec,
                "--threads", "1",
                "--bfile", bfile_prefix,
                "--chr", chrom,
                "--from-bp", str(int(start)),
                "--to-bp", str(int(end)),
                "--extract", extract_path,
                "--ld-snp", str(index_snp),
                "--r2",
                "--ld-window", "999999",
                "--ld-window-kb", str(kb_span),
                "--ld-window-r2", "0",
                "--out", out_prefix
            ]
            if chrom == "X":
                cmd.extend(["--allow-extra-chr"])

            res = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            if res.returncode != 0:
                log(f"PLINK failed for {index_snp}:\n{res.stdout}")
                continue

            ld_maps[r2_col] = _read_plink_ld_table(out_prefix + ".ld", index_snp)

    return ld_maps, index_status, tuple(sorted(ref_snps))


def build_ld_annot_for_window(window_union_df, ld_maps, ref_snps, idx1_ref, idx2_ref, idx3_ref):
    out = window_union_df[["REF_SNP", "CHR", "BP"]].drop_duplicates("REF_SNP").copy()
    ref_set = set(ref_snps)

    out["in_ref"] = out["REF_SNP"].notna() & out["REF_SNP"].isin(ref_set)
    out["r2_1"] = out["REF_SNP"].map(ld_maps.get("r2_1", {}))
    out["r2_2"] = out["REF_SNP"].map(ld_maps.get("r2_2", {}))
    out["r2_3"] = out["REF_SNP"].map(ld_maps.get("r2_3", {}))

    if idx1_ref in ref_set:
        out.loc[out["REF_SNP"] == idx1_ref, "r2_1"] = 1.0
    if idx2_ref in ref_set:
        out.loc[out["REF_SNP"] == idx2_ref, "r2_2"] = 1.0
    if idx3_ref in ref_set:
        out.loc[out["REF_SNP"] == idx3_ref, "r2_3"] = 1.0

    return out[["REF_SNP", "r2_1", "r2_2", "r2_3", "in_ref"]]


def merge_ld_annot(df, ld_annot):
    df = dedup_columns(df)
    ld_annot = dedup_columns(ld_annot)

    base = df.drop(columns=["r2_1", "r2_2", "r2_3"], errors="ignore").copy()
    out = base.merge(ld_annot, on="REF_SNP", how="left")
    out["in_ref"] = out["in_ref"].fillna(False).astype(bool)
    return dedup_columns(out)
