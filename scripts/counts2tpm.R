library(IOBR)
counts <- readRDS("data/exp_counts.rds") %>% as.data.frame()
rownames(counts) <- counts$Gene_ID
counts$Gene_ID <- NULL
counts <- as.matrix(counts)
tpm <- count2tpm(countMat = counts, source = "local", idType = "Ensembl")
saveRDS(tpm,file = "data/exp_tpm.rds")