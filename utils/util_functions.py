import polars as pl
import pandas as pd



def create_max_spliceai_vep_ploras(df):
    df = df.with_columns(
        max_spliceai_score = pl.when((pl.col("SpliceAI_pred").is_not_null()) & (pl.col("SpliceAI_pred") != "-") )
        .then(
            pl.col("SpliceAI_pred")
            .str.split("|")
            .list.slice(1, 4) 
            .list.eval(pl.element().cast(pl.Float64))
            .list.max()
        )
        .otherwise(pl.lit(0.0)) # Assign 0.0 if missing, so the filter can handle it
    )
    return df


def unique_vep_res_per_gene(df):
    impact_map = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    df['impact_rank'] = df['IMPACT'].map(impact_map)
    
    df = df[df["Gene"] != "-"]
    df = df.sort_values(
        by=['Gene', '#Uploaded_variation', 'impact_rank', 'max_spliceai_score'], 
        ascending=[True, True, False, False]
    )
    if 'sampleID' in df.columns:
        df_final = df.drop_duplicates(
            subset=['sampleID', 'Gene', '#Uploaded_variation']
        ).drop(columns="impact_rank")
    else:
        df_final = df.drop_duplicates(
            subset=['Gene', '#Uploaded_variation']
        ).drop(columns="impact_rank")
    return df_final


def left_align_trim_vcf(row):
    ref, alt = row["REF"], row["ALT"]

    # 1. Find and strip common prefix
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    if i > 0:
        ref = ref[i:]
        alt = alt[i:]

    # 2. Find and strip common suffix
    j = 0
    while (
        j < len(ref)
        and j < len(alt)
        and ref[len(ref) - 1 - j] == alt[len(alt) - 1 - j]
    ):
        j += 1
    if j > 0:
        ref = ref[: len(ref) - j]
        alt = alt[: len(alt) - j]

    # 3. Represent empty strings as '-'
    ref = ref if ref != "" else "-"
    alt = alt if alt != "" else "-"

    return pd.Series([ref, alt])

def compute_vep_location(row):
    chrom = row["#CHROM"]
    pos = int(row["POS"])
    ref = row["REF"]
    alt = row["ALT"]

    # VEP-style insertion
    if ref == "-" and alt != "-":
        return f"{chrom}:{pos}-{pos+1}"

    # VEP-style deletion
    if alt == "-" and ref != "-":
        if len(ref) == 1:
            return f"{chrom}:{pos}"
        return f"{chrom}:{pos+1}-{pos + 1 + len(ref) - 1}"

    # SNV
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:{pos}"

    # fallback
    span = max(len(ref), len(alt))
    end = pos + span - 1
    return f"{chrom}:{pos}" if end == pos else f"{chrom}:{pos}-{end}"



def get_midpoint(pos_str):
    # Ensure it's a string to avoid errors with potential NaNs
    pos_str = str(pos_str)
    
    if '-' in pos_str:
        # Split by the hyphen and calculate the average
        parts = pos_str.split('-')
        start, end = int(parts[0]), int(parts[1])
        return (start + end) // 2
    else:
        # If it's just a single number, return it as a float/int
        return float(pos_str)

def to_vep_alleles(row):
    ref, alt = str(row["Ref"]), str(row["Alt"])

    # Find the length of the common prefix
    common_len = 0
    for r, a in zip(ref, alt):
        if r == a:
            common_len += 1
        else:
            break

    # If they share a prefix, slice it off
    if common_len > 0:
        ref_trimmed = ref[common_len:]
        alt_trimmed = alt[common_len:]

        # Replace empty strings with VEP's hyphen '-'
        ref = ref_trimmed if ref_trimmed else "-"
        alt = alt_trimmed if alt_trimmed else "-"

    return pd.Series([ref, alt])


def vep_normalize_row(row):
    """Trims padding bases from REF and ALT to mimic VEP normalization."""
    ref = row["VEP_Ref"]
    alt = row["VEP_Alt"]

    # 1. Strip identical trailing (right) bases
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # 2. Strip identical leading (left) bases
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]

    # 3. Replace empty strings with VEP's '-'
    return pd.Series({"VEP_Ref": ref if ref else "-", "VEP_Alt": alt if alt else "-"})