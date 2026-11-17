## some outliers + VUS were missing acmg annotations. there were pulled manually from FRANKLIN and are updated using this 

import pandas as pd



def add_missing_acmgs(df):

    # 1. Raw variant lines
    raw_data = """
    # chr9:111642368:G>A : Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr9-111642368-G-A
    # chr1-155208082-G-A  Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr1-155208082-G-A 
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27019512-G-A
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr11-94170372-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr9-111651708-T-A
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr7-75684236-G-GC
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr22-41565354-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr4-84384070-A-T
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr10-99258202-C-T
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr10-99220447-C-T
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr9-111685190-C-T
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr11-94180454-G-A
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr7-75684236-G-GC
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr15-89864214-G-A
    # Likely Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr1-27684931-G-A
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr9-111665142-C-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr12-13062413-G-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr16-72825380-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27019512-G-A
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr7-75684236-G-GC
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr7-75693391-T-G
    # Likely Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr7-75694276-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr7-75693391-T-G
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr6-100016213-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr20-5946059-A-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27683248-C-T
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr15-89862456-T-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27690404-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-23219372-G-A
    # Likely Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr1-27686731-C-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr15-89861200-A-G
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27684507-CACA-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr12-13044800-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27683351-C-T
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr22-19245125-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr20-5946059-A-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27683084-T-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr15-89861200-A-G
    # Likely Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr1-27686731-C-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr17-8054145-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr6-100016213-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr9-22008798-T-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27687322-T-C
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr14-68806097-G-T
    # VUS (Likly benign) https://franklin.genoox.com/clinical-db/variant/snp/chr1-27683216-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr12-13044800-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr7-75693391-T-G
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27682647-C-G
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27683351-C-T
    # Pathogenic https://franklin.genoox.com/clinical-db/variant/snp/chr15-89862456-T-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr1-27687628-G-A
    # VUS https://franklin.genoox.com/clinical-db/variant/snp/chr12-21624056-T-C
    """

    # 2. Parse lines into a mapping DataFrame
    records = []
    for line in raw_data.strip().split("\n"):
        if not line.strip():
            continue
    
        # Extract variant from the end of the URL
        url_part = line.split()[-1]
        variant_str = url_part.split("/")[-1]  # e.g., 'chr9-111642368-G-A'
        chrom, pos, ref, alt = variant_str.split("-")
        chrom = chrom.split("chr")[1]
    
        # Map classification label
        if "Likely Pathogenic" in line or "Likely Pathogenic" in line:
            classification = "Likely Pathogenic"
        elif "Pathogenic" in line:
            classification = "Pathogenic"
        elif "VUS (LB)" in line:
            classification = "Likely benign"
        elif "VUS" in line:
            classification = "Uncertain Significance"
        else:
            classification = None

        df.loc[(df["#CHROM"] == chrom) & (df["POS"] == int(pos)) & (df["REF"] == ref) & (df["ALT"] == alt), "germline_ACMG_Classification"] = classification 

    return df
