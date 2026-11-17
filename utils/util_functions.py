import polars as pl
import pandas as pd
import re



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



def calculate_proportions(df, expression_col="CNV", variant_type="germline", omic_type = "gene"):
    # Work on a copy to prevent SettingWithCopyWarning and index fragmentation
    df = df.copy()
    
    for expression_direction in ["overexpression", "underexpression"]:
        vus_col = f"VUS_{expression_direction}_{expression_col}"
        df[vus_col] = False

        if expression_direction == "underexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=True) 
        elif expression_direction == "overexpression":
            if (df["Method"] == "sf-normalized zScores").any():
                df = df.sort_values("zScore", ascending=False)

        # 1. Define VUS based on CNV/VEP type
        if expression_col == "CNV":
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True) 
                if variant_type == "germline":
                    df.loc[df[expression_col].str.contains("Germline_DEL", na=False), vus_col] = True
                else:
                    df.loc[df[expression_col].str.contains("Somatic_HDEL", na=False), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df[expression_col].str.contains("AMP|DUP", na=False), vus_col] = True
                
        elif expression_col == "VEP":
            df.loc[(df[f"{variant_type}_IMPACT"].str.contains("HIGH")), vus_col] = True
            if omic_type == "protein":
                df.loc[(df[f"{variant_type}_Consequence"].str.contains("missense")), vus_col] = True
            
        elif expression_col == "CNV_VEP":
            # Combined logic: either high impact VEP or relevant CNV
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True)
                df.loc[df["CNV"].str.contains("DEL", na=False) | (df["IMPACT"] == "HIGH"), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df["CNV"].str.contains("AMP|DUP", na=False) | (df["IMPACT"] == "HIGH"), vus_col] = True
        
        else: # Fallback for other custom expression_cols
            if expression_direction == "underexpression":
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=True)   
                df.loc[df[expression_col].str.contains("DEL", na=False), vus_col] = True
            else:
                if (df["Method"] == "sf-normalized zScores").any():
                    df = df.sort_values("zScore", ascending=False)
                df.loc[df[expression_col].str.contains("AMP|DUP", na=False), vus_col] = True
                df.loc[(df["IMPACT"] == "HIGH"), vus_col] = True

        # 2. Subset based on zScore direction for proportion calculation
        if expression_direction == "underexpression":
            mask = df['zScore'] < 0
        else:
            mask = df['zScore'] > 0
            
        # Create a temporary subset to calculate running statistics
        # This ensures we don't have indexing conflicts with the main DF
        subset = df.loc[mask].copy()
        
        if len(subset) > 0:
            # Calculate ranks and cumulative counts
            u_rank = np.arange(1, len(subset) + 1)
            u_vus = (subset[vus_col] == True).to_numpy()
            u_x = np.cumsum(u_vus)
            
            # Vectorized Clopper–Pearson interval calculation
            alpha = 0.05
            # beta.ppf handles the statistics for the CI
            ci_min = beta.ppf(alpha/2, u_x, u_rank - u_x + 1)
            ci_max = beta.ppf(1 - alpha/2, u_x + 1, u_rank - u_x)

            # Map the calculated values back to the original 'df' using the mask
            df.loc[mask, f"{expression_direction}_rank"] = u_rank
            df.loc[mask, f"{expression_direction}_proportions_{expression_col}_{variant_type}"] = u_x / u_rank
            df.loc[mask, f"{expression_direction}_ci_min_{expression_col}_{variant_type}"] = ci_min
            df.loc[mask, f"{expression_direction}_ci_max_{expression_col}_{variant_type}"] = ci_max

    return df





non_coding_consequence_list = ["consequence_intron_variant", "consequence_upstream_gene_variant", 
                    "consequence_downstream_gene_variant", "consequence_3_prime_utr_variant"]


coding_consequence_list = [
    # HIGH impacts
    'transcript_ablation',
    'stop_gained',
    'stop_lost',
    'missense_variant',
    'frameshift_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'transcript_amplification',
    'start_lost',
    'feature_elongation',
    'feature_truncation',
    
    # Other coding variants
    'inframe_deletion',
    'inframe_insertion',
    'synonymous_variant',
    'coding_sequence_variant',
    'stop_retained_variant',
    'incomplete_terminal_codon_variant',
]

# Create a regex pattern: "transcript_ablation|stop_gained|stop_lost|..."
coding_pattern = "|".join(map(re.escape, set(coding_consequence_list)))
    

 
 # 'consequence_protein_altering_variant',

 # 'consequence_splice_donor_5th_base_variant',
 # 'consequence_splice_donor_region_variant',
 # 'consequence_splice_polypyrimidine_tract_variant',
 # 'consequence_splice_region_variant',
 # 'consequence_start_retained_variant',





def _false_discovery_control(ps, *, axis=0, method='bh'):
    # Input Validation and Special Cases
    ps = np.asarray(ps)

    ps_in_range = (np.issubdtype(ps.dtype, np.number)
                   and np.all(ps == np.clip(ps, 0, 1)))
    if not ps_in_range:
        raise ValueError("`ps` must include only numbers between 0 and 1.")

    methods = {'bh', 'by'}
    if method.lower() not in methods:
        raise ValueError(f"Unrecognized `method` '{method}'."
                         f"Method must be one of {methods}.")
    method = method.lower()

    if axis is None:
        axis = 0
        ps = ps.ravel()

    axis = np.asarray(axis)[()]
    if not np.issubdtype(axis.dtype, np.integer) or axis.size != 1:
        raise ValueError("`axis` must be an integer or `None`")

    if ps.size <= 1 or ps.shape[axis] <= 1:
        return ps[()]

    ps = np.moveaxis(ps, axis, -1)
    m = ps.shape[-1]

    # Main Algorithm
    # Equivalent to the ideas of [1] and [2], except that this adjusts the
    # p-values as described in [3]. The results are similar to those produced
    # by R's p.adjust.

    # "Let [ps] be the ordered observed p-values..."
    order = np.argsort(ps, axis=-1)
    ps = np.take_along_axis(ps, order, axis=-1)  # this copies ps

    # Equation 1 of [1] rearranged to reject when p is less than specified q
    i = np.arange(1, m + 1)
    ps *= m / i

    # Theorem 1.3 of [2]
    if method == 'by':
        ps *= np.sum(1 / i)

    # accounts for rejecting all null hypotheses i for i < k, where k is
    # defined in Eq. 1 of either [1] or [2]. See [3]. Starting with the index j
    # of the second to last element, we replace element j with element j+1 if
    # the latter is smaller.
    np.minimum.accumulate(ps[..., ::-1], out=ps[..., ::-1], axis=-1)

    # Restore original order of axes and data
    np.put_along_axis(ps, order, values=ps.copy(), axis=-1)
    ps = np.moveaxis(ps, -1, axis)

    return np.clip(ps, 0, 1)


