library(UniProt.ws)
library(data.table)
library(biomaRt)


up <- UniProt.ws(taxId = 9606)

protein_intensities <- fread("/omics/odcf/analysis/hipo/hipo_021/proteomics_TUM/MASTER_preprocessed_fp_with_annot_scores_298.csv")

View(protein_intensities[,c("Fasta headers", "fasta_gene_name", "Majority protein IDs")])

extract_uniprot_entry <- function(headers) {
  h <- trimws(sub("^>", "", as.character(headers)))
  out <- rep(NA_character_, length(h))

  sptr <- grepl("^(?:sp|tr)\\|", h)
  # From sp|...|... or tr|...|..., take the 3rd field up to a space or '|'
  out[sptr] <- sub("^(?:sp|tr)\\|[^|]+\\|([^\\s|]+).*", "\\1", h[sptr], perl = TRUE)

  # Otherwise, look for a NAME_ORG token anywhere in the line
  m <- regexpr("[A-Z0-9]+_[A-Z0-9]+", h[!sptr], perl = TRUE)
  out[!sptr][m > 0] <- regmatches(h[!sptr], m)

  out
}

protein_intensities[, uniprotID := extract_uniprot_entry(`Fasta headers`)]
protein_intensities[, fasta_gene_name := sub(".*GN=([^ ]+).*", "\\1", `Fasta headers`)]


map_symbols_to_ensembl <- function(symbols, mart) {
  orig <- data.table(input_symbol = toupper(trimws(as.character(symbols))),
                     idx = seq_along(symbols))
  syms <- unique(orig$input_symbol[nzchar(orig$input_symbol)])
  if (!length(syms)) {
    return(data.table(input_symbol = character(), ensembl_gene_id = character(),
                      gene_biotype = character(), chromosome_name = character()))
  }
  
  # attrs (gene-level only)
  attrs <- c("hgnc_symbol","external_gene_name","ensembl_gene_id",
             "gene_biotype","chromosome_name")
  
  # helper: collapse duplicates per query_symbol with deterministic priorities
  collapse_stage <- function(dt) {
    if (!nrow(dt)) return(dt)
    # normalize chromosome (handles chrMT/chrX/etc.)
    dt[, chromosome_name := as.character(chromosome_name)]
    dt[, chrom_clean := sub("^chr", "", chromosome_name, ignore.case = TRUE)]
    # chromosome rank: 1..22 < X < MT < Y < others
    dt[, chrom_rank := fifelse(chrom_clean %chin% as.character(1:22), as.integer(chrom_clean),
                               fifelse(chrom_clean=="X", 23L,
                                       fifelse(chrom_clean=="MT",24L,
                                               fifelse(chrom_clean=="Y", 90L, 99L))))]
    # biotype preference: protein_coding (others allowed but lower)
    prio <- c("protein_coding","IG_C_gene","TR_C_gene")
    dt[, gene_biotype_rank := fifelse(gene_biotype %chin% prio,
                                      match(gene_biotype, prio), 999L)]
    setorder(dt, query_symbol, gene_biotype_rank, chrom_rank, ensembl_gene_id)
    out <- dt[, .SD[1], by = query_symbol]
    out[, chromosome_name := chrom_clean][, chrom_clean := NULL]
    out[]
  }
  
  # ----------------
  # Stage 1: exact symbol/display name
  # ----------------
  s1 <- as.data.table(getBM(attributes = attrs,
                            filters    = "external_gene_name",
                            values     = syms,
                            mart       = mart))
  if (nrow(s1)) s1[, query_symbol := toupper(external_gene_name)]
  s1 <- s1[query_symbol %chin% syms]
  s1c <- collapse_stage(s1)
  
  mapped1 <- unique(s1c$query_symbol)
  remaining1 <- setdiff(syms, mapped1)
  
  # ----------------
  # Stage 2: HGNC symbol (only for remaining)
  # ----------------
  s2 <- if (length(remaining1)) as.data.table(getBM(
    attributes = attrs, filters = "hgnc_symbol", values = remaining1, mart = mart
  )) else data.table()
  if (nrow(s2)){
    s2[, query_symbol := toupper(hgnc_symbol)]
    s2 <- s2[query_symbol %chin% remaining1]
    s2c <- collapse_stage(s2)
  } else{
    cols <- colnames(s1c)
    s2c <- dt <- as.data.table(setNames(replicate(length(cols), character(), simplify = FALSE), cols))
  }
  
  mapped2 <- unique(s2c$query_symbol)
  remaining2 <- setdiff(remaining1, mapped2)
  
  # ----------------
  # Stage 3: synonyms (only for still-unmapped)
  # ----------------
  attrs_syn <- c(attrs, "external_synonym")
  s3 <- if (length(remaining2)) as.data.table(getBM(
    attributes = attrs_syn, filters = "external_synonym", values = remaining2, mart = mart
  )) else data.table()
  if (nrow(s3)) s3[, query_symbol := toupper(external_synonym)]
  s3 <- s3[query_symbol %chin% remaining2]
  s3c <- collapse_stage(s3)
  
  s3c[query_symbol == "SMAP", ensembl_gene_id := "ENSG00000110696"] ## from genecards for  C11orf58
  s3c[query_symbol == "TAZ", ensembl_gene_id := "ENSG00000102125"] 
  s3c[query_symbol == "KIAA1107", ensembl_gene_id := "ENSG00000069712"]

  # combine, preserving stage priority (s1 > s2 > s3)
  combined <- rbindlist(list(
    s1c[, .(query_symbol, ensembl_gene_id, gene_biotype, chromosome_name)],
    s2c[, .(query_symbol, ensembl_gene_id, gene_biotype, chromosome_name)],
    s3c[, .(query_symbol, ensembl_gene_id, gene_biotype, chromosome_name)]
  ), use.names = TRUE, fill = TRUE)
  
  # keep first occurrence per query_symbol (earlier stage wins)
  setkey(combined, query_symbol)
  combined <- combined[!duplicated(query_symbol)]
  
  # stitch back to original order; unmapped stay NA
  final <- merge(orig[, .(input_symbol, idx)],
               combined[, .(input_symbol = query_symbol,
                            ensembl_gene_id, gene_biotype, chromosome_name)],
               by = "input_symbol", all.x = TRUE, sort = FALSE)[order(idx)]

  
  manual_naming <- data.table::data.table(
    gene_name = c(
      "CDKN2A-P16","L1RE1","LOC102724159","IGHV1-8","AKAP2","IGHV3-9","TRA","HCG_1984214",
      "GAGE6","CNK3/IPCEF1","CDKN2A-P14","IGHV3-38-3","IGHV4-38-2","IGLV5-39","CLDN18-2",
      "IGHV4-30-2","ERVS71-1","LINC00696","LOC100996750","TSPY26P","PALM2","OVOS2","FAM243B"
    ),
    ensembl_gene_id = c(
      "ENSG00000147889_16", NA, NA, NA, "ENSG00000241978", "ENSG00000211940", NA, NA,
      NA, NA, "ENSG00000147889_14", NA, NA, NA, NA,
      NA, NA, NA, NA, "ENSG00000235217", "ENSG00000243444", "ENSG00000177359", "ENSG00000277277"
    )
  )
  
  
  
  final <- left_join(final, manual_naming, by=c("input_symbol"="gene_name"))
  final[is.na(ensembl_gene_id.x), ensembl_gene_id.x := ensembl_gene_id.y]
  setnames(final, "ensembl_gene_id.x", "ensembl_gene_id")
  return (final)
}






protein_intensities[, unique_protein_id := sub("\\s*;.*$", "", `Majority protein IDs`)]
protein_intensities[, unique_protein_id := toupper(trimws(as.character(unique_protein_id)))]
protein_intensities[, `Majority protein IDs` := toupper(trimws(as.character(`Majority protein IDs`)))]
protein_intensities[, fasta_gene_name := toupper(trimws(as.character(fasta_gene_name)))]

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")  # pick a fixed version if needed

dt <- map_symbols_to_ensembl(protein_intensities[, fasta_gene_name], ensembl)
dt <- dt[!is.na(ensembl_gene_id)]

dt_res <- merge(protein_intensities, dt, by.x = "fasta_gene_name", by.y = "input_symbol", all.x=T)
dt_res <- dt_res[!is.na(ensembl_gene_id)]



dt_res[, ensembl_gene_id.y := NULL]
dt_res[, idx := NULL]

View(dt_res[,c("Fasta headers", "fasta_gene_name", "Majority protein IDs", "ensembl_gene_id")])

fwrite(dt_res, "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protein_intensities_ensembl_majority.tsv", sep="\t")


gene_annot_dt <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_data/preprocess/v19/gene_name_mapping_v19.tsv")
gene_annot_dt <- gene_annot_dt %>%
  mutate(geneID_short = sub("\\..*", "", gene_id))
old_res_names <- fread("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/sample_data/protrider_intensities.tsv")
old_res_names <- left_join(old_res_names, gene_annot_dt, by=c("protein_ID"= "gene_name")) 
