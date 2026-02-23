
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/Esophagus_Stomach/results_per_junction.tsv'),
    output = list('/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/results/v19/fraser/Esophagus_Stomach/results_gene_all.tsv'),
    params = list('Esophagus_Stomach', "cohort" = 'Esophagus_Stomach'),
    wildcards = list('Esophagus_Stomach', "cohort" = 'Esophagus_Stomach'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 70000, "mem_mib" = 66758, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/local/a379i/47076705/cluster_tmp'),
    config = list(),
    rule = 'splicing_dna_first',
    bench_iteration = as.numeric(NA),
    scriptdir = '/home/a379i/Scripts/AberrantSplicing',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(FRASER)
library(dplyr)
library(BiocParallel)

register(MulticoreParam(8))

make_subsets_for_FDR_from_list <- function(genes, sampleIDs,
                                           subset_name = "Genes_to_test_on_all_samples") {

    stopifnot(is.character(genes))
    stopifnot(is.character(sampleIDs))

    sub_ls <- rep(list(genes), length(sampleIDs))
    names(sub_ls) <- sampleIDs

    subsets_expanded <- list()
    subsets_expanded[[subset_name]] <- sub_ls

    return(subsets_expanded)
}

source("~/Scripts/utils/load_gtf_cgc_dresden.R")
source("~/Scripts/utils/drop_utils.R")


genes_to_subset <- dresden_list

cohort <- "ACC"
cohort <- snakemake@params[["cohort"]]

fds <- loadFraserDataSet(dir="/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_splicing/datasets", name=paste(cohort, "v19", sep = '--'))

sample_ids <- colData(fds)$sampleID

subsets <- make_subsets_for_FDR_from_list(
    genes = dresden_list,
    sampleIDs = sample_ids
)
subsets <- subsets$Genes_to_test_on_all_samples
s <- subsets[1]

names(s) <- "Genes_to_test_on_all_samples"

type = "jaccard"
fds <- calculatePadjValues(fds, type=type, subsets=s)

fds <- saveFraserDataSet(fds)



snakemake <- readRDS(paste0("/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/.drop/tmp/AS/", cohort, "--v19/08_results.Rds"))



source(snakemake@input$setup, echo=FALSE)
source(snakemake@input$add_HPO_cols)
library(AnnotationDbi)

annotation    <- snakemake@wildcards$annotation
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load fds and create a new one
fds <- loadFraserDataSet(dir=workingDir, name=paste(dataset, annotation, sep = '--'))

# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=snakemake@params$padjCutoff,
                    deltaPsiCutoff=snakemake@params$deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
all_junc_res  <- as.data.table(results(fds, all=TRUE))
print('Results per junction extracted')

# Add features
if(nrow(res_junc_dt) > 0){
    # number of samples per gene and variant
    res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
    res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
    res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
} else{
    warning("The aberrant splicing pipeline gave 0 intron-level results for the ", dataset, " dataset.")
}

# Extract full results by gene
res_gene <- results(fds, psiType=psiTypes,
                    aggregate=TRUE, collapse=FALSE,
                    all=TRUE)
res_genes_dt   <- as.data.table(res_gene)
print('Results per gene extracted')
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_full)

# Subset gene results to aberrant
padj_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
res_genes_dt <- res_genes_dt[do.call(pmin, c(res_genes_dt[,padj_cols, with=FALSE], 
                                                list(na.rm = TRUE))) <= snakemake@params$padjCutoff &
                                    abs(deltaPsi) >= snakemake@params$deltaPsiCutoff & 
                                    totalCounts >= 5,]

if(nrow(res_genes_dt) > 0){
    # add HPO overlap information
    sa <- fread(snakemake@config$sampleAnnotation, 
                colClasses = c(RNA_ID = 'character', DNA_ID = 'character'))
    if(!is.null(sa$HPO_TERMS)){
        if(!all(is.na(sa$HPO_TERMS)) & ! all(sa$HPO_TERMS == '')){
            res_genes_dt <- add_HPO_cols(res_genes_dt, hpo_file = snakemake@params$hpoFile)
        }
    }
} else{
    warning("The aberrant splicing pipeline gave 0 gene-level results for the ", dataset, " dataset.")
}

# Annotate results with spliceEventType and blacklist region overlap
txdb <- loadDb(snakemake@input$txdb)
    
# annotate the type of splice event and UTR overlap
if(nrow(res_junc_dt) > 0){
    res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
}
if(nrow(res_genes_dt) > 0){
    res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)
}
    
# set genome assembly version to load correct blacklist region BED file (hg19 or hg38)
assemblyVersion <- snakemake@config$genomeAssembly
if(grepl("grch37", assemblyVersion, ignore.case=TRUE)) assemblyVersion <- "hg19"
if(grepl("grch38", assemblyVersion, ignore.case=TRUE)) assemblyVersion <- "hg38"

# annotate overlap with blacklist regions
if(assemblyVersion %in% c("hg19", "hg38")){
    if(nrow(res_junc_dt) > 0){
        res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, 
                                        assemblyVersion=assemblyVersion)
    }
    if(nrow(res_genes_dt) > 0){
        res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, 
                                         assemblyVersion=assemblyVersion)
    }
} else{
    message(date(), ": cannot annotate blacklist regions as no blacklist region\n", 
            "BED file is available for genome assembly version ", assemblyVersion, 
            " as part of FRASER.")
}

# Results
options(scipen=999)
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene_aberrant)
fwrite(all_junc_res, sep='\t', file=snakemake@output$resultTableJunc_full)
