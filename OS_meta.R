library(GEOmirror)
library(GEOquery)
library(limma)
library(stringr)
library(AnnoProbe)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(metafor)
library(data.table)

options(timeout = 100000)

# 下载并提取数据集
gse1 <- getGEO("GSE56815", GSEMatrix = TRUE)
gse2 <- getGEO("GSE7158", GSEMatrix = TRUE)
gse3 <- getGEO("GSE56116", GSEMatrix = TRUE)
gse4 <- getGEO("GSE56814", GSEMatrix = TRUE)
gse5 <- getGEO("GSE35958", GSEMatrix = TRUE)
gse6 <- getGEO("GSE35956", GSEMatrix = TRUE)

# 假设每个数据集中只有一个表达矩阵
exprs1 <- exprs(gse1[[1]])
exprs2 <- exprs(gse2[[1]])
exprs3 <- exprs(gse3[[1]])
exprs4 <- exprs(gse4[[1]])
exprs5 <- exprs(gse5[[1]])
exprs6 <- exprs(gse6[[1]])

# 标准化（如果需要）
exprs1 <- normalizeBetweenArrays(exprs1)
exprs2 <- normalizeBetweenArrays(exprs2)
exprs3 <- normalizeBetweenArrays(exprs3)
exprs4 <- normalizeBetweenArrays(exprs4)
exprs5 <- normalizeBetweenArrays(exprs5)
exprs6 <- normalizeBetweenArrays(exprs6)

#limma
group1 = c(rep('control',40),rep('osteoporosis',40))
group1 = factor(group1,levels = c('control','osteoporosis'))
design1 <- model.matrix(~ group1) # 根据实际分组调整
fit1 <- lmFit(exprs1, design1)
fit1 <- eBayes(fit1)
topTable1 <- topTable(fit1, adjust="fdr", number=nrow(exprs1))

group2 = c(rep('control',14),rep('osteoporosis',12))
group2 = factor(group2,levels = c('control','osteoporosis'))
design2 <- model.matrix(~ group2) # 根据实际分组调整
fit2 <- lmFit(exprs2, design2)
fit2 <- eBayes(fit2)
topTable2 <- topTable(fit2, adjust="fdr", number=nrow(exprs2))

pd3 = pData(gse3[[1]])
group3 = c(rep('osteoporosis',4),'control',rep('osteoporosis',6),rep('control',2))
group3 = factor(group3,levels = c('control','osteoporosis'))
design3 <- model.matrix(~ group3) # 根据实际分组调整
fit3 <- lmFit(exprs3, design3)
fit3 <- eBayes(fit3)
topTable3 <- topTable(fit3, adjust="fdr", number=nrow(exprs3))

pd4 = pData(gse4[[1]]) %>%
  mutate(group = if_else(str_detect(pd4$`bone mineral density:ch1`,'high'),'control','osteoporosis'))
group4 = pd4$group
group4 = factor(group4,levels = c('control','osteoporosis'))
design4 <- model.matrix(~ group4) # 根据实际分组调整
fit4 <- lmFit(exprs4, design4)
fit4 <- eBayes(fit4)
topTable4 <- topTable(fit4, adjust="fdr", number=nrow(exprs4))

group5 = c(rep('control',4),rep('osteoporosis',5))
group5 = factor(group5,levels = c('control','osteoporosis'))
design5 <- model.matrix(~ group5) # 根据实际分组调整
fit5 <- lmFit(exprs5, design5)
fit5 <- eBayes(fit5)
topTable5 <- topTable(fit5, adjust="fdr", number=nrow(exprs5))

group6 = c(rep('control',5),rep('osteoporosis',5))
group6 = factor(group6,levels = c('control','osteoporosis'))
design6 <- model.matrix(~ group6) # 根据实际分组调整
fit6 <- lmFit(exprs6, design6)
fit6 <- eBayes(fit6)
topTable6 <- topTable(fit6, adjust="fdr", number=nrow(exprs6))
#准备meta分析的数据
# annotation
ids1 = AnnoProbe::idmap('GPL96')
ids1 = ids1[!duplicated(ids1$symbol),]

ids2 = AnnoProbe::idmap('GPL570')
ids2 = ids2[!duplicated(ids2$symbol),]

ids3 = fread('./GPL4133-12599.txt')%>%
  select(ID,GENE_SYMBOL)
ids3 = ids3[!duplicated(ids3$GENE_SYMBOL),]
rownames(ids3) = ids3$ID
ids3$ID = as.character(ids3$ID)

ids4 = fread('./GPL5175-3188.txt')%>%
  select(ID,gene_assignment)
temp = str_split(ids4$gene_assignment,'//')
temp <- lapply(temp, `[`, 2)
temp <- do.call(rbind, temp)
temp = trimws(temp)
ids4$gene_assignment = temp
ids4 = ids4[!duplicated(ids4$gene_assignment),]
rownames(ids4) = ids4$ID
ids4$ID = as.character(ids4$ID)

ids5 = AnnoProbe::idmap('GPL570')
ids5 = ids5[!duplicated(ids5$symbol),]

ids6 = AnnoProbe::idmap('GPL570')
ids6 = ids6[!duplicated(ids6$symbol),]
#选出共同基因
topTable1 = mutate(topTable1,probe_id = rownames(topTable1)) %>%
  inner_join(ids1)
topTable2 = mutate(topTable2,probe_id = rownames(topTable2)) %>%
  inner_join(ids2)
topTable3 = mutate(topTable3,ID = rownames(topTable3)) %>%
  inner_join(ids3)
topTable4 = mutate(topTable4,ID = rownames(topTable4)) %>%
  inner_join(ids4)
topTable5 = mutate(topTable5,probe_id = rownames(topTable5)) %>%
  inner_join(ids5)
topTable6 = mutate(topTable6,probe_id = rownames(topTable6)) %>%
  inner_join(ids6)


rownames(topTable1) = topTable1$symbol
rownames(topTable2) = topTable2$symbol
rownames(topTable3) = topTable3$GENE_SYMBOL
rownames(topTable4) = topTable4$gene_assignment
rownames(topTable5) = topTable5$symbol
rownames(topTable6) = topTable6$symbol

genes = unique(c(rownames(topTable1),rownames(topTable2),rownames(topTable3),rownames(topTable4),rownames(topTable5),rownames(topTable6)))

# 提取效应大小和标准误
effect1 <- topTable1$logFC
se1 <- topTable1$logFC / topTable1$t
effect2 <- topTable2$logFC
se2 <- topTable2$logFC / topTable2$t
effect3 <- topTable3$logFC
se3 <- topTable3$logFC / topTable3$t
effect4 <- topTable4$logFC
se4 <- topTable4$logFC / topTable4$t
effect5 <- topTable5$logFC
se5 <- topTable5$logFC / topTable5$t
effect6 <- topTable6$logFC
se6 <- topTable6$logFC / topTable6$t

# 合并成一个数据框
meta_data <- data.frame(
  gene = genes,
  effect1 = NA,
  se1 = NA,
  effect2 = NA,
  se2 = NA,
  effect3 = NA,
  se3 = NA,
  effect4 = NA,
  se4 = NA,
  effect5 = NA,
  se5 = NA,
  effect6 = NA,
  se6 = NA
)

meta_data$effect1[match(rownames(topTable1), genes)] <- effect1
meta_data$se1[match(rownames(topTable1), genes)] <- se1
meta_data$effect2[match(rownames(topTable2), genes)] <- effect2
meta_data$se2[match(rownames(topTable2), genes)] <- se2
meta_data$effect3[match(rownames(topTable3), genes)] <- effect3
meta_data$se3[match(rownames(topTable3), genes)] <- se3
meta_data$effect4[match(rownames(topTable4), genes)] <- effect4
meta_data$se4[match(rownames(topTable4), genes)] <- se4
meta_data$effect5[match(rownames(topTable5), genes)] <- effect5
meta_data$se5[match(rownames(topTable5), genes)] <- se5
meta_data$effect6[match(rownames(topTable6), genes)] <- effect6
meta_data$se6[match(rownames(topTable6), genes)] <- se6

# 进行meta分析
results <- data.frame()
for(i in 1:nrow(meta_data)) {
  yi = c(meta_data$effect1[i],meta_data$effect2[i],meta_data$effect3[i],meta_data$effect4[i],meta_data$effect5[i],meta_data$effect6[i])
  sei = c(meta_data$se1[i],meta_data$se2[i],meta_data$se3[i],meta_data$se4[i],meta_data$se5[i],meta_data$se6[i])
  valid <- !is.na(yi) & !is.na(sei)
  if(sum(valid) > 1) { # 至少需要两个有效值
    res <- rma(yi = yi[valid], sei = sei[valid], method = "REML")
    results <- rbind(results, data.frame(gene = meta_data$gene[i], 
                                         est = res$beta, 
                                         se = res$se, 
                                         zval = res$zval, 
                                         pval = res$pval, 
                                         ci.lb = res$ci.lb, 
                                         ci.ub = res$ci.ub))
  }
}

results$fdr <- p.adjust(results$pval, method = "fdr")
write.csv(results,file = './meta_result.csv')

#森林图
OS = fread('./apoptosis.csv') %>%
  filter(`Relevance score` > 0.5)
OS = OS$`Gene Symbol`
results_fil = results[results$gene %in% OS,]
results_fil = results_fil[results_fil$fdr < 0.05,]
write.csv(results_fil,file = './meta_apoptosis.csv')
os_sig = results_fil$gene
os_ensg = bitr(os_sig, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
os_ensg = os_ensg$ENSEMBL
write.table(os_ensg,file = 'os_ensg.txt',row.names = F,quote = F,col.names = F)

results_fil = results_fil %>%
  mutate(abs = abs(est))%>%
  top_n(20,wt =abs)

forest(
  x = results_fil$est,
  sei = results_fil$se,
  slab = results_fil$gene,
  xlab = "Effect Size",
  xlim = c(-1, 1),
  at = c(-1, -0.5, 0, 0.5, 1),
  ylim = c(-1, nrow(results_fil) + 5),
  alim = c(-1, 1),
  cex = 1,
  psize = 1,
  main = "Top 20 Effect Size Forest Plot of Meta-Analysis Results",
  col = "black"
)

