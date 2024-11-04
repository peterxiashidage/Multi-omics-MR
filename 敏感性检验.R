library(TwoSampleMR)
library(dplyr)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)

dats = list.files('./eqtl_split/') %>%
  {.[str_detect(.,'f_')]}
gene = str_split(dats, "_", simplify = FALSE) %>%
  {unlist(sapply(., `[`, 2))}

total = list()
diymr = function(dats) {
  files = fread(paste0('./eqtl_split/',dats,sep = ''))
  mrResult=mr(files)
  mrTab=generate_odds_ratios(mrResult)
  all = inner_join(mrResult,mrTab)
  return(all)
}  

for (i in seq_along(dats)) {
  total[[gene[i]]] = diymr(dats[i])
}


df_list <- lapply(names(total), function(name) {
  df <- total[[name]]
  df$name <- name
  return(df)
})
final_df <- do.call(rbind, df_list)

final_df = final_df%>%
  filter(pval < 0.05)

ids = bitr(final_df$name,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db) %>%
  mutate(name = ENSEMBL)
final_df = inner_join(final_df,ids)

sig_mr_ensg = unique(final_df$ENSEMBL)
sig_mr_gene = unique(final_df$SYMBOL)

write.table(sig_mr_ensg,file = './sig_mr_ensg.txt',quote = F,row.names = F,col.names = F)
write.table(sig_mr_gene,file = './sig_mr_gene.txt',quote = F,row.names = F,col.names = F)
