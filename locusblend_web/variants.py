"""Variant and summary-statistic normalization helpers.

Pure pandas transformations extracted from app.py. Nothing in this module
touches Streamlit, session state, PLINK, the filesystem, or gene/LD data.
"""

import pandas as pd

from locusblend_web.references import normalize_chrom


def dedup_columns(df):
    return df.loc[:, ~pd.Index(df.columns).duplicated()].copy()


def _clean_locus_df(df, source_name="uploaded data"):
    df = dedup_columns(df).copy()

    # ---------- normalize raw column names ----------
    df.columns = pd.Index(df.columns).map(
        lambda x: str(x).replace("\ufeff", "").strip()
    )

    # Keep the original column names for debugging.
    original_cols = list(df.columns)

    # Build a case-insensitive lookup table.
    lower_to_actual = {}
    for c in df.columns:
        cl = str(c).strip().lower()
        if cl not in lower_to_actual:
            lower_to_actual[cl] = c

    def pick(*aliases):
        """Return the first existing actual column name from aliases (case-insensitive)."""
        for a in aliases:
            key = str(a).strip().lower()
            if key in lower_to_actual:
                return lower_to_actual[key]
        return None

    rename_map = {}

    # ---------- standard required fields ----------
    c = pick("CHR", "chr", "#chr", "chrom", "chromosome")
    if c and c != "CHR":
        rename_map[c] = "CHR"

    c = pick("BP", "bp", "pos", "position", "base_pair_location")
    if c and c != "BP":
        rename_map[c] = "BP"

    c = pick("P", "p", "pval", "pvalue", "p_value", "P-value", "P_VALUE")
    if c and c != "P":
        rename_map[c] = "P"

    # ---------- rsid / snp ----------
    c = pick("rsid", "RSID", "SNP", "snp", "MarkerName", "markername", "ID", "id")
    if c and c != "rsid":
        rename_map[c] = "rsid"

    # ---------- alleles ----------
    # This dataset uses ALLELE1 / ALLELE0.
    c = pick("A1", "a1", "EA", "ea", "effect_allele", "effect allele", "ALLELE1", "allele1")
    if c and c != "A1":
        rename_map[c] = "A1"

    c = pick(
        "A2",
        "a2",
        "NEA",
        "nea",
        "other_allele",
        "non_effect_allele",
        "non effect allele",
        "ALLELE0",
        "allele0",
    )
    if c and c != "A2":
        rename_map[c] = "A2"

    # ---------- optional/common fields ----------
    c = pick("BETA", "beta", "Effect", "effect", "estimate")
    if c and c != "BETA":
        rename_map[c] = "BETA"

    c = pick("SE", "se", "StdErr", "stderr", "standard_error")
    if c and c != "SE":
        rename_map[c] = "SE"

    c = pick("A1FREQ", "a1freq", "EAF", "eaf", "MAF", "maf", "freq")
    if c and c != "A1FREQ":
        rename_map[c] = "A1FREQ"

    c = pick("N", "n", "N_incl", "n_incl", "samplesize", "sample_size")
    if c and c != "N":
        rename_map[c] = "N"

    df = df.rename(columns=rename_map)

    # ---------- clean values ----------
    if "CHR" in df.columns:
        df["CHR"] = (
            df["CHR"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.strip()
        )

    for col in ["BP", "P", "BETA", "SE", "A1FREQ", "N"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ["A1", "A2"]:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype(str)
                .str.strip()
                .str.upper()
            )

    if "rsid" in df.columns:
        df["rsid"] = df["rsid"].astype(str).str.strip()

    # ---------- required columns ----------
    required = ["CHR", "BP", "P", "A1", "A2"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"[DEBUG] {source_name} original columns: {original_cols}", flush=True)
        print(f"[DEBUG] {source_name} normalized columns: {list(df.columns)}", flush=True)
        raise ValueError(f"{source_name} missing columns: {missing}")

    # ---------- optional helper column ----------
    # Build a CHR:BP posID when it is missing.
    if "posID" not in df.columns:
        df["posID"] = df["CHR"].astype(str) + ":" + df["BP"].astype("Int64").astype(str)

    # Used later for uniqueid matching.
    if "uniqueid" not in df.columns:
        df["uniqueid"] = (
            df["CHR"].astype(str)
            + ":"
            + df["BP"].astype("Int64").astype(str)
            + ":"
            + df["A1"].astype(str)
            + ":"
            + df["A2"].astype(str)
        )

    # ---------- DISPLAY_ID (for UI / hover) ----------
    if "DISPLAY_ID" not in df.columns:
        if "rsid" in df.columns:
            df["DISPLAY_ID"] = df["rsid"].astype(str).str.strip()
        elif "SNP" in df.columns:
            df["DISPLAY_ID"] = df["SNP"].astype(str).str.strip()
        elif "uniqueid" in df.columns:
            df["DISPLAY_ID"] = df["uniqueid"].astype(str).str.strip()
        elif "posID" in df.columns:
            df["DISPLAY_ID"] = df["posID"].astype(str).str.strip()
        else:
            df["DISPLAY_ID"] = (
                df["CHR"].astype(str)
                + ":"
                + pd.to_numeric(df["BP"], errors="coerce").astype("Int64").astype(str)
            )

    return df


def chrom_mask(df, chrom):
    """Return a boolean mask matching df['CHR'] to the selected chromosome.

    Both sides are normalised (chr prefix stripped) before comparison.
    """
    target = normalize_chrom(chrom)
    if "CHR" not in df.columns:
        raise KeyError(f"DataFrame has no 'CHR' column; columns: {list(df.columns)}")
    return df["CHR"].map(normalize_chrom) == target


def attach_reference_snp_two_pass(df, bim_ref):
    """
    Match input summary stats to BIM reference in two passes:
    1) CHR:BP:A1:A2
    2) CHR:BP:A2:A1

    Returns original df plus:
      - QUERY_UID
      - REF_UID
      - REF_UID_FLIP
      - REF_MATCH
      - REF_MATCH_TYPE   ("forward", "flip", or NA)
      - NEED_FLIP        (True if matched on swapped alleles)
    """

    out = df.copy()
    ref = bim_ref.copy()

    # ---------- helper ----------
    def norm_chr(s):
        return s.map(normalize_chrom)

    def norm_allele(s):
        return (
            s.astype(str)
             .str.strip()
             .str.upper()
        )

    # ---------- make sure df has required cols ----------
    required_df = ["CHR", "BP", "A1", "A2"]
    miss_df = [c for c in required_df if c not in out.columns]
    if miss_df:
        raise ValueError(f"input df missing columns for reference matching: {miss_df}")

    out["CHR"] = norm_chr(out["CHR"])
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce")
    out["A1"] = norm_allele(out["A1"])
    out["A2"] = norm_allele(out["A2"])

    # Use Int64 so NA values do not crash the pipeline.
    out["BP"] = out["BP"].astype("Int64")

    # ---------- build QUERY_UID on df ----------
    out["QUERY_UID"] = (
        out["CHR"].astype(str) + ":"
        + out["BP"].astype(str) + ":"
        + out["A1"].astype(str) + ":"
        + out["A2"].astype(str)
    )

    out["QUERY_UID_FLIP"] = (
        out["CHR"].astype(str) + ":"
        + out["BP"].astype(str) + ":"
        + out["A2"].astype(str) + ":"
        + out["A1"].astype(str)
    )

    # ---------- normalize bim_ref ----------
    # Accept common BIM column names.
    rename_map = {}
    if "#CHROM" in ref.columns and "CHR" not in ref.columns:
        rename_map["#CHROM"] = "CHR"
    if "chr" in ref.columns and "CHR" not in ref.columns:
        rename_map["chr"] = "CHR"
    if "pos" in ref.columns and "BP" not in ref.columns:
        rename_map["pos"] = "BP"
    if "bp" in ref.columns and "BP" not in ref.columns:
        rename_map["bp"] = "BP"
    if "a1" in ref.columns and "A1" not in ref.columns:
        rename_map["a1"] = "A1"
    if "a2" in ref.columns and "A2" not in ref.columns:
        rename_map["a2"] = "A2"
    if "snp" in ref.columns and "SNP" not in ref.columns:
        rename_map["snp"] = "SNP"
    if "id" in ref.columns and "SNP" not in ref.columns:
        rename_map["id"] = "SNP"

    ref = ref.rename(columns=rename_map)

    required_ref = ["CHR", "BP", "A1", "A2"]
    miss_ref = [c for c in required_ref if c not in ref.columns]
    if miss_ref:
        raise ValueError(f"bim_ref missing columns for reference matching: {miss_ref}")

    ref["CHR"] = norm_chr(ref["CHR"])
    ref["BP"] = pd.to_numeric(ref["BP"], errors="coerce").astype("Int64")
    ref["A1"] = norm_allele(ref["A1"])
    ref["A2"] = norm_allele(ref["A2"])

    # SNP name from the BIM file.
    if "SNP" not in ref.columns:
        ref["SNP"] = (
            ref["CHR"].astype(str) + ":"
            + ref["BP"].astype(str) + ":"
            + ref["A1"].astype(str) + ":"
            + ref["A2"].astype(str)
        )

    # Forward UID.
    ref["REF_UID"] = (
        ref["CHR"].astype(str) + ":"
        + ref["BP"].astype(str) + ":"
        + ref["A1"].astype(str) + ":"
        + ref["A2"].astype(str)
    )

    # Reverse UID.
    ref["REF_UID_FLIP"] = (
        ref["CHR"].astype(str) + ":"
        + ref["BP"].astype(str) + ":"
        + ref["A2"].astype(str) + ":"
        + ref["A1"].astype(str)
    )

    # ---------- first pass: forward ----------
    m1 = (
        ref[["REF_UID", "SNP"]]
        .drop_duplicates("REF_UID")
        .rename(columns={
            "REF_UID": "QUERY_UID",
            "SNP": "REF_MATCH_FWD",
        })
    )

    out = out.merge(m1, on="QUERY_UID", how="left")

    # ---------- second pass: allele flipped ----------
    m2 = (
        ref[["REF_UID_FLIP", "SNP"]]
        .drop_duplicates("REF_UID_FLIP")
        .rename(columns={
            "REF_UID_FLIP": "QUERY_UID",
            "SNP": "REF_MATCH_REV",
        })
    )

    out = out.merge(
        m2,
        left_on="QUERY_UID",
        right_on="QUERY_UID",
        how="left",
    )

    # ---------- combine ----------
    out["REF_MATCH"] = out["REF_MATCH_FWD"]
    out.loc[out["REF_MATCH"].isna(), "REF_MATCH"] = out.loc[
        out["REF_MATCH"].isna(), "REF_MATCH_REV"
    ]

    out["REF_MATCH_TYPE"] = pd.NA
    out.loc[out["REF_MATCH_FWD"].notna(), "REF_MATCH_TYPE"] = "forward"
    out.loc[
        out["REF_MATCH_FWD"].isna() & out["REF_MATCH_REV"].notna(),
        "REF_MATCH_TYPE",
    ] = "flip"

    out["NEED_FLIP"] = out["REF_MATCH_TYPE"].eq("flip")

    # backward compatibility
    out["REF_SNP"] = out["REF_MATCH"]

    # keep REF_UID columns for downstream use/debugging
    out["REF_UID"] = out["QUERY_UID"]
    out["REF_UID_FLIP"] = out["QUERY_UID_FLIP"]

    return out


def resolve_index_variant_from_input(df_top_ref, df_bottom_ref, user_text):
    """
    Resolve user-entered SNP text against dataset 1 (top) and dataset 2 (bottom).

    Matching priority:
      1) DISPLAY_ID
      2) rsid
      3) SNP
      4) REF_MATCH
      5) uniqueid
      6) QUERY_UID

    Returns one matched row (as Series).
    """

    user_text = str(user_text).strip()
    if user_text == "":
        raise ValueError("Empty index SNP input.")

    def prep(df, source_label):
        x = df.copy()

        # Cast common columns to string when present.
        for col in ["DISPLAY_ID", "rsid", "SNP", "REF_MATCH", "uniqueid", "QUERY_UID"]:
            if col in x.columns:
                x[col] = x[col].astype(str).str.strip()

        # Build a uniqueid when it is missing.
        if "uniqueid" not in x.columns and all(c in x.columns for c in ["CHR", "BP", "A1", "A2"]):
            x["uniqueid"] = (
                x["CHR"].map(normalize_chrom)
                + ":"
                + pd.to_numeric(x["BP"], errors="coerce").astype("Int64").astype(str)
                + ":"
                + x["A1"].astype(str).str.strip().str.upper()
                + ":"
                + x["A2"].astype(str).str.strip().str.upper()
            )

        # Generate DISPLAY_ID automatically.
        if "DISPLAY_ID" not in x.columns:
            if "rsid" in x.columns:
                x["DISPLAY_ID"] = x["rsid"]
            elif "SNP" in x.columns:
                x["DISPLAY_ID"] = x["SNP"]
            elif "REF_MATCH" in x.columns:
                x["DISPLAY_ID"] = x["REF_MATCH"]
            elif "uniqueid" in x.columns:
                x["DISPLAY_ID"] = x["uniqueid"]
            elif "QUERY_UID" in x.columns:
                x["DISPLAY_ID"] = x["QUERY_UID"]
            else:
                x["DISPLAY_ID"] = pd.NA

        x["SOURCE"] = source_label
        return x

    top = prep(df_top_ref, "top")
    bottom = prep(df_bottom_ref, "bottom")
    merged = pd.concat([top, bottom], axis=0, ignore_index=True, sort=False)

    # Strip whitespace consistently.
    for col in ["DISPLAY_ID", "rsid", "SNP", "REF_MATCH", "uniqueid", "QUERY_UID"]:
        if col in merged.columns:
            merged[col] = merged[col].astype(str).str.strip()

    # Try each candidate column in order.
    search_cols = ["DISPLAY_ID", "rsid", "SNP", "REF_MATCH", "uniqueid", "QUERY_UID"]

    for col in search_cols:
        if col in merged.columns:
            hit = merged[merged[col].astype(str) == user_text].copy()
            if len(hit) > 0:
                # Prefer top, then bottom; this rule can be changed.
                hit["_priority"] = hit["SOURCE"].map({"top": 0, "bottom": 1}).fillna(9)
                hit = hit.sort_values(["_priority"]).drop(columns=["_priority"])
                return hit.iloc[0]

    # Then match case-insensitively again (mainly for rsid / SNP).
    user_upper = user_text.upper()
    for col in search_cols:
        if col in merged.columns:
            hit = merged[merged[col].astype(str).str.upper() == user_upper].copy()
            if len(hit) > 0:
                hit["_priority"] = hit["SOURCE"].map({"top": 0, "bottom": 1}).fillna(9)
                hit = hit.sort_values(["_priority"]).drop(columns=["_priority"])
                return hit.iloc[0]

    available = [c for c in search_cols if c in merged.columns]
    raise ValueError(
        f"Could not resolve index SNP '{user_text}'. "
        f"Searched columns: {available}"
    )


def attach_uploaded_ld_keys(df, ld_universe):
    """
    Uploaded-LD equivalent of attach_reference_snp_two_pass: set REF_SNP
    on each row to whichever identifier column first appears in the
    uploaded LD universe. Does not require a 1000G BIM.

    Match priority: DISPLAY_ID, rsid, SNP, uniqueid, posID, QUERY_UID,
    QUERY_UID_FLIP. Rows with no identifier in ld_universe get
    REF_SNP=NaN and render as grey/missing downstream.
    """
    out = dedup_columns(df.copy())

    if all(c in out.columns for c in ["CHR", "BP", "A1", "A2"]):
        chr_s = out["CHR"].map(normalize_chrom)
        bp_s = pd.to_numeric(out["BP"], errors="coerce").astype("Int64").astype(str)
        a1_s = out["A1"].astype(str).str.strip().str.upper()
        a2_s = out["A2"].astype(str).str.strip().str.upper()
        uid_fwd = chr_s + ":" + bp_s + ":" + a1_s + ":" + a2_s
        uid_rev = chr_s + ":" + bp_s + ":" + a2_s + ":" + a1_s
        if "uniqueid" not in out.columns:
            out["uniqueid"] = uid_fwd
        if "QUERY_UID" not in out.columns:
            out["QUERY_UID"] = uid_fwd
        if "QUERY_UID_FLIP" not in out.columns:
            out["QUERY_UID_FLIP"] = uid_rev

    universe = {str(x) for x in (ld_universe or ())}

    candidate_cols = ["DISPLAY_ID", "rsid", "SNP", "uniqueid", "posID", "QUERY_UID", "QUERY_UID_FLIP"]
    present_cols = [c for c in candidate_cols if c in out.columns]

    ref_snp = pd.Series([pd.NA] * len(out), index=out.index, dtype=object)
    for col in present_cols:
        col_str = out[col].astype(str).str.strip()
        mask = ref_snp.isna() & col_str.isin(universe)
        if mask.any():
            ref_snp.loc[mask] = col_str.loc[mask]

    out["REF_SNP"] = ref_snp
    return out
