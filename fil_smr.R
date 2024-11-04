library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

datas = fread('./eqtl_fil_os.msmr') %>%
  filter(p_SMR_multi < 0.05 & p_HEIDI > 0.05)

gene = datas$probeID 

gene = bitr(gene,fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db) 
gene = mutate(gene,probeID = ENSEMBL)


datas = inner_join(datas,gene) 
sigg = gene$SYMBOL
sige = gene$ENSEMBL

write.table(sigg,file = '../sig_gene.txt',quote = F,row.names = F,col.names = F)
write.table(sige,file = '../sig_ensg.txt',quote = F,row.names = F,col.names = F)
