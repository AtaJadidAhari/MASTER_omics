library(MultiAssayExperiment)
library(RaggedExperiment)


load("/omics/odcf/analysis/OE0246_projects/datamaster/dataMASTER_release/object/dataMASTER.RData")

# Subset all data types by samples which have RNA data.
idx = !is.na(colData(dataMASTER)$RNASample)
myobj = dataMASTER[,idx,]
# 



master_sa <- as.data.table(colData(myobj), keep.rownames = TRUE)

names(master_sa)

master_sa <- master_sa[, c("PatientID", "TumorID", "ControlID", "TumorCellContent", "SufficientTumorCellContent", "DNASeq", 
              "RNASample", "Sex", "Comment", "TissueMaterial", "RNASampleID"
              )]
master_sa[, sampleID := rownames(colData(myobj))]

master_sa
