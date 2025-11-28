library(dplyr)
tcga_exp <- data.table::fread("data/TCGA-BRCA.star_tpm.tsv.gz",data.table = F)
mapping <- data.table::fread("data/gencode.v36.annotation.gtf.gene.probemap",data.table = F)
tcga_exp <- left_join(tcga_exp %>% rename(id = Ensembl_ID),
                      mapping %>% select(id,gene)) %>% 
  select(gene,everything()) %>% select(-id) %>% 
  distinct(gene,.keep_all = T)

tcga_sur <- data.table::fread("data/TCGA-BRCA.survival.tsv.gz",data.table = F)
need_sample <- colnames(tcga_exp)[2:ncol(tcga_exp)]
need_sample <- need_sample[which(as.numeric(substr(need_sample,14,15)) < 11)]  
need_sample <- need_sample[which(substr(need_sample,1,12) %in% tcga_sur$`_PATIENT`)]
exp_filter <- tcga_exp %>% 
  select(gene, need_sample)
sample_info <- left_join(
  data.frame(samples = colnames(exp_filter)[2:ncol(exp_filter)]) %>% 
    mutate(patient = substr(samples,1,12)),
  tcga_sur %>% select(`_PATIENT`,OS.time,OS) %>% rename(patient = `_PATIENT`) %>% distinct_all()
)
saveRDS(exp_filter,"data/tcga_brca_filter_exp.rds")
saveRDS(sample_info,"data/tcga_brca_sampleinfo.rds")