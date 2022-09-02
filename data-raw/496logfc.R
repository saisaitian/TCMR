load("F:/Github/TCMR/logfc20210824.Rdata")

logfc[is.na(logfc)]=0
rownames(logfc) <- logfc$gene_id
TCM_RNAseq <- logfc[,-c(1:2)]

usethis::use_data(TCM_RNAseq)
