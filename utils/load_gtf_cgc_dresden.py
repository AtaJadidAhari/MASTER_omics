import pandas as pd

# --- Load annotation table ---
gene_annot_dt = pd.read_csv(
    "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv",
    sep=","
)

# Create geneID_short (strip version suffix)
gene_annot_dt["geneID_short"] = gene_annot_dt["gene_id"].str.replace(r"\..*", "", regex=True)

# --- Load CGC table ---
cgc = pd.read_csv(
    "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/resources/Cosmic_CancerGeneCensus_v102_GRCh37.tsv",
    sep="\t"
)

# Inner join by gene_name = GENE_SYMBOL
cgc = gene_annot_dt.merge(cgc, left_on="gene_name", right_on="GENE_SYMBOL", how="inner")

# Rename gene_id → geneID
cgc = cgc.rename(columns={"gene_id": "geneID"})

# Filter TSG and oncogenes
cgc_tsg = cgc[cgc["ROLE_IN_CANCER"].str.contains("TSG", na=False)]
cgc_oncogene = cgc[cgc["ROLE_IN_CANCER"].str.contains("oncogene", na=False)]

# --- Dresden list ---
dresden_list = [
    "AIP", "AKT1", "ALK", "ANKRD26", "APC", "ATM", "BAP1", "BARD1", "BLM", "BMPR1A",
    "BRCA1", "BRCA2", "BRIP1", "BUB1B", "CBL", "CDC73", "CDH1", "CDK4", "CDKN1B",
    "CDKN1C", "CDKN2A", "CEBPA", "CHEK2", "COL7A1", "CREBBP", "CYLD", "DDB2",
    "DDX41", "DICER1", "DIS3L2", "DKC1", "DOCK8", "EGFR", "ELANE", "EPCAM",
    "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ETV6", "EXT1", "EXT2", "EZH2",
    "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI",
    "FANCL", "FAS", "FH", "FLCN", "GATA2", "GJB2", "GPC3", "HFE", "HMBS", "HPS1",
    "HRAS", "IDH1", "KIT", "KRAS", "MAX", "MEN1", "MET", "MITF", "MLH1", "MSH2",
    "MSH6", "MUTYH", "NBN", "NF1", "NF2", "NRAS", "NSD1", "NTHL1", "PALB2",
    "PARN", "PDGFRA", "PHOX2B", "PIK3CA", "PIK3R1", "PMS2", "POLD1", "POLE",
    "POLH", "POT1", "PRKAR1A", "PRSS1", "PTCH1", "PTEN", "PTPN11", "RAD51C",
    "RAD51D", "RB1", "RECQL4", "RET", "RHBDF2", "RMRP", "RPS19", "RPS24", "RPS26",
    "RTEL1", "RUNX1", "SAMD9", "SBDS", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD",
    "SETBP1", "SLX4", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SOS1", "SPINK1",
    "STAT3", "STK11", "SUFU", "TERC", "TERT", "TINF2", "TMEM127", "TP53", "TRIM37",
    "TSC1", "TSC2", "VHL", "WAS", "WRAP53", "WRN", "WT1", "XPA", "XPC", "XRCC2"
]

# Create DataFrame and merge with annotation
dresden_dt = pd.DataFrame({"gene_name": dresden_list})
dresden_dt["predisposition_gene"] = True
dresden_dt = dresden_dt.merge(gene_annot_dt, on="gene_name", how="left")

dresden_dt["geneID"] = dresden_dt["gene_id"]

# Add CGC roles
dresden_dt_cgc = dresden_dt.merge(cgc[["geneID", "ROLE_IN_CANCER"]], on="geneID", how="left")

dresden_dt_tsg = dresden_dt_cgc[dresden_dt_cgc["ROLE_IN_CANCER"].str.contains("TSG", na=False)]
dresden_dt_oncogene = dresden_dt_cgc[dresden_dt_cgc["ROLE_IN_CANCER"].str.contains("oncogene", na=False)]

# --- Extended Dresden list ---
extended_dresden_list = [
    "ABRAXAS1", "ACD", "AIP", "AKT1", "AKT2", "ALK", "ANKRD26", "APC", "AR", "ARMC5", "ATG12", "ATM",
  "ATR", "ATRIP", "ATRX", "AURKA", "AXIN2", "BAP1", "BARD1", "BIK", "BLM", "BMPR1A", "BPTF", "BRAF",
  "BRCA1", "BRCA2", "BRIP1", "BUB1", "BUB1B", "CASR", "CBL", "CDC73", "CDH1", "CDK4", "CDKN1B",
  "CDKN1C", "CDKN2A", "CEBPA", "CEP57", "CHEK2", "CMTR2", "COL7A1", "CREBBP", "CTC1", "CTLA4",
  "CTNNA1", "CTNNB1", "CTR9", "CYLD", "CYP2D6", "DDB2", "DDX41", "DHCR24", "DHCR7", "DICER1",
  "DIS3L2", "DKC1", "DLST", "DOCK8", "DROSHA", "EGFR", "EGLN1", "EGLN2", "ELANE", "ELP1",
  "EPAS1", "EPCAM", "ERBB2", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6L2", "ETV6", "EXT1",
  "EXT2", "EZH2", "FAN1", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG",
  "FANCI", "FANCL", "FANCM", "FAS", "FH", "FLCN", "FLT3", "FOXE1", "G6PC1", "GALNT12",
  "GATA1", "GATA2", "GBA1", "GJB2", "GNA11", "GNAQ", "GPC3", "GPR101", "GPR161", "GREM1",
  "HAVCR2", "HAX1", "HMBS", "HOXB13", "HPS1", "HRAS", "IDH1", "IKZF1", "ITK", "KIF1B",
  "KIT", "KMT2D", "KRAS", "LEMD3", "LIG4", "LZTR1", "MAD1L1", "MAP2K2", "MAP3K1", "MAX",
  "MBD4", "MCM8", "MECOM", "MEN1", "MET", "MITF", "MLH1", "MLH3", "MNX1", "MPL",
  "MRE11", "MSH2", "MSH3", "MSH6", "MTAP", "MUTYH", "NBN", "NF1", "NF2", "NOP10",
  "NOTCH3", "NPAT", "NPC1", "NRAS", "NSD1", "NTHL1", "PALB2", "PALLD", "PARN", "PAX5",
  "PDGFRA", "PDGFRB", "PHOX2B", "PIK3CA", "PIK3R1", "PLA2G2A", "PMS2", "POLD1", "POLE", "POLH",
  "POT1", "POU6F2", "PPM1D", "PPP1CB", "PRF1", "PRKAR1A", "PRSS1", "PTCH1", "PTCH2", "PTEN",
  "PTPN11", "PTPRJ", "RABL3", "RAD50", "RAD51", "RAD51C", "RAD51D", "RB1", "RECQL", "RECQL4",
  "RECQL5", "REST", "RET", "RHBDF2", "RINT1", "RMRP", "RNASEL", "RNF43", "ROS1", "RPS19",
  "RPS20", "RPS24", "RPS26", "RTEL1", "RUNX1", "SAMD9", "SAMD9L", "SASH1", "SBDS", "SDHA",
  "SDHAF2", "SDHB", "SDHC", "SDHD", "SEC23B", "SEMA4A", "SETBP1", "SH2B3", "SH2D1A", "SHOC2",
  "SLC25A11", "SLX4", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SOS1", "SPINK1", "SPRED1", "SPRTN",
  "STAT3", "STK11", "SUFU", "TBXT", "TERC", "TERT", "TG", "TINF2", "TMEM127", "TP53",
  "TREX1", "TRIM28", "TRIM37", "TRIP13", "TSC1", "TSC2", "TSHR", "UBA1", "UBE2T", "VHL",
  "WAS", "WRAP53", "WRN", "WT1", "XPA", "XPC", "XRCC2"
]



# Create DataFrame and merge with annotation
extended_dresden_dt = pd.DataFrame({"gene_name": extended_dresden_list})
extended_dresden_dt["extended_predisposition_gene"] = True
extended_dresden_dt = extended_dresden_dt.merge(gene_annot_dt, on="gene_name", how="left")

extended_dresden_dt["geneID"] = extended_dresden_dt["gene_id"]
extended_dresden_dt_cgc = extended_dresden_dt.merge(cgc[["geneID", "ROLE_IN_CANCER"]], on="geneID", how="left")

extended_dresden_dt_cgc["Predisposition"] = "Extended list"
extended_dresden_dt_cgc.loc[extended_dresden_dt_cgc["gene_name"].isin(dresden_dt["gene_name"]), "Predisposition"] = "Predisposition"



cgc["Predisposition"] = False
cgc.loc[cgc["geneID_short"].isin(extended_dresden_dt["geneID_short"]), "Predisposition"] = "Extended predisp"
cgc.loc[cgc["geneID_short"].isin(dresden_dt["geneID_short"]), "Predisposition"] = "Predisposition"


