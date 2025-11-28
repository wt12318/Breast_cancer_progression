library(dplyr)
# .libPaths(c("/data/sdd/wt/miniconda3/envs/deconvolution/lib/R/library/",
#             "/home/wt/R/x86_64-pc-linux-gnu-library/4.2",
#             "/data/sdd/wt/miniconda3/envs/R_pack/lib/R/library/"))
library(immunedeconv)
tpm <- readRDS("data/exp_tpm.rds")
dt <- tpm
dt$genes <- rownames(dt)
dt <- dt %>% select(genes,everything())
write.table(dt,"data/cibersort_input.txt",quote = F,row.names = F,sep = "\t")

est <- deconvolute_estimate(tpm)
est$cell_type <- paste0("estimate~",rownames(est))
est <- est %>% select(cell_type,everything())

t1 <- deconvolute(tpm, "quantiseq")t1 <- decoestnvolute(tpm, "quantiseq")
t2 <- deconvolute(tpm, "timer",indications=rep("BRCA",ncol(tpm)))
t3 <- deconvolute(tpm, "consensus_tme",indications=rep("BRCA",ncol(tpm))) 
t4 <- deconvolute(tpm, "mcp_counter")
t5 <- xCell::xCellAnalysis(tpm)
t5 <- as.data.frame(t5)
t5$cell_type <- rownames(t5)
t5 <- t5 %>% select(cell_type,everything())
t6 <- deconvolute(tpm, "epic")
t7 <- deconvolute(tpm, "abis")

all_immune <- bind_rows(
  est,
  t1 %>% mutate(cell_type = paste0("quantiseq~",cell_type)),
  t2 %>% mutate(cell_type = paste0("timer~",cell_type)),
  t3 %>% mutate(cell_type = paste0("consensus_tme~",cell_type)),
  t4 %>% mutate(cell_type = paste0("mcp_counter~",cell_type)),
  t5 %>% mutate(cell_type = paste0("xCell~",cell_type)),
  t6 %>% mutate(cell_type = paste0("epic~",cell_type)),
  t7 %>% mutate(cell_type = paste0("abis~",cell_type)),
)
rownames(all_immune) <- 1:nrow(all_immune)
###cibersort
cibersort <- data.table::fread("data/cibersort_res.txt",data.table = F)
cibersort <- cibersort %>% select(-c(24,25,26))
rownames(cibersort) <- cibersort$Mixture
cibersort$Mixture <- NULL
cibersort <- as.data.frame(t(cibersort))
cibersort$cell_type <- rownames(cibersort)
cibersort <- cibersort %>% select(cell_type, everything()) %>% 
  mutate(cell_type = paste0("cibersort~",cell_type))

cibersort_abs <- data.table::fread("data/cibersort_res_absolute.txt",data.table = F)
cibersort_abs <- cibersort_abs %>% select(-c(24,25,26,27))
rownames(cibersort_abs) <- cibersort_abs$Mixture
cibersort_abs$Mixture <- NULL
cibersort_abs <- as.data.frame(t(cibersort_abs))
cibersort_abs$cell_type <- rownames(cibersort_abs)
cibersort_abs <- cibersort_abs %>% select(cell_type, everything()) %>% 
  mutate(cell_type = paste0("cibersort_abs~",cell_type))

all_immune <- bind_rows(
  all_immune,
  cibersort,
  cibersort_abs
)
rownames(all_immune) <- 1:nrow(all_immune)
saveRDS(all_immune, "data/all_immune_8.rds")