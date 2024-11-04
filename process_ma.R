library(dplyr)
library(data.table)

datas = fread('./BEurope-Bmd-As-C-Gwas-SumStats.txt')
datas = datas %>%
  select(-SNPID,-CHR,-BP,-INFO) %>%
  rename(A1 = ALLELE1,
         A2 = ALLELE0,
         freq = A1FREQ,
         b = BETA,
         p = P,
         n = N) %>%
  filter(!grepl("^\\.", SNP) & SNP != ".")

write.table(datas,'Osteoporosis_GWAS.ma',row.names = F,quote = F)
