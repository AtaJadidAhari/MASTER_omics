library(OUTRIDER)


res_dir <- "/omics/odcf/analysis/hipo/hipo_021/outlier_analysis/drop_runs/drop_master_202502_allGenes/processed_results/aberrant_expression/v19/outrider"
folders <- list.dirs(res_dir, full.names = TRUE, recursive = FALSE)


robust_z <- function(x) {
        med <- median(x, na.rm = TRUE)
        mad_val <- mad(x, na.rm = TRUE)  # constant=1 = unscaled MAD
        if (mad_val == 0) return(rep(0.0005, length(x)))   # avoid division by zero
        return(0.6745 *(x - med) / mad_val)
}

calculate_zScores <- function(ods){
    drop_group_name <- tail(unlist(strsplit(drop_group, "/")), 1)
    ods <- readRDS(paste0(drop_group, "/ods.Rds"))
    #fwrite(zScore_res_sorted, paste0(drop_group, "/zScores.tsv"))
    log2_fpkms <- log2(fpkm(ods) + 1)
    log2_fpkm_zcores <- as.data.table(t(apply(log2_fpkms, 1, scale)))

    log2_fpkm_robust_zcores <- as.data.table(t(apply(log2_fpkms, 1, robust_z))) 
    setnames(log2_fpkm_robust_zcores, colnames(log2_fpkms))
    
    rownames(log2_fpkm_zcores) <- rownames(log2_fpkms)
    setnames(log2_fpkm_zcores, colnames(log2_fpkms))
    
    size_factors <- sizeFactors(ods)
    log2_fpkm_zcores_sizeFactor_normalized <- sweep(log2_fpkm_zcores, 2, sizeFactors(ods), FUN = "/")
    
    
    log2_fpkm_zcores_sizeFactor_normalized <- as.data.table(log2_fpkm_zcores_sizeFactor_normalized)
    log2_fpkm_zcores_sizeFactor_normalized[, "geneID" := rownames(log2_fpkms)]
    log2_fpkm_robust_zcores[, geneID := rownames(log2_fpkms)]
    zScore_res <- melt(log2_fpkm_zcores_sizeFactor_normalized, value.name = "zScore", variable.name = "sampleID")

    log2_fpkm_robust_zcores <- as.data.table(log2_fpkm_robust_zcores)
    robust_res <- melt(log2_fpkm_robust_zcores, id.vars = "geneID",
                       value.name = "robust_zScore", variable.name = "sampleID")
    
    # Merge both by geneID and sampleID
    zScore_res <- merge(zScore_res, robust_res, by = c("geneID", "sampleID"))
    
    zScore_res[, abs_zScore := abs(zScore)]
    zScore_res_sorted <- zScore_res[order(zScore_res$abs_zScore, decreasing = TRUE), ]

    zScore_res_sorted[, DROP_GROUP:= drop_group_name]
    
    p_values <- 2 * pnorm(-abs(zScore_res_sorted$abs_zScore))
    zScore_res_sorted$pvalue <- p_values
    
    p_adj <- p.adjust(p_values, method = "BY")
    zScore_res_sorted$padjust <- p_adj

    return (zScore_res_sorted)

}


for (drop_group in folders){
    drop_group_name <- tail(unlist(strsplit(drop_group, "/")), 1)
    ods <- readRDS(paste0(drop_group, "/ods.Rds"))
    zScore_res_sorted <- calculate_zScores(ods)
    fwrite(zScore_res_sorted, paste0(drop_group, "/zScores.tsv"))
}






