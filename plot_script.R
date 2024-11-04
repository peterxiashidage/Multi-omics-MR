library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(karyoploteR)
library(ggplotify)
library(showtext)
source('./plot/plot_SMR.r',regular="path/to/Times New Roman.ttf")

font_add(family = 'Times New Roman',regular = "./Times New Roman.ttf")

fread('../mqtl/bl_mqtl_lite_chr2.epi') %>%
  .$V2 %>% unique() %>% write.table('./m_chr2_cg',sep ='\t',quote = F,row.names = F,col.names = F)

#做一个SMRdata对象
fread('./LRRFIP1.txt') %>%
  dplyr::select(SNP,A1,A2,Freq,b,SE,p) %>%
  {colnames(.)=c('SNP','A1','A2','freq','b','se','p');.} %>%
  mutate(n=142487) %>%
  write.table('./LRRFIP1.ma',quote = F,sep = '\t',row.names = F,col.names = T)

#画eqtl-op
SMRData = ReadSMRData("./LRRFIP1_op.ENSG00000124831.txt")
pdf('./LRRFIP1_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="LRRFIP1_Osteoporosis",
              cisWindow = 1000,transWindow = 1000)
dev.off()



SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)

#画met-op
SMRData = ReadSMRData("./plot.cg14458575.txt")
pdf('./cg14458575_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="cg14458575_Osteoporosis",
              xlabel = 'mQTL effect sizes',
              ylabel = 'GWAS effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()

#画met-eqtl
SMRData = ReadSMRData("./cg14458575_LRRFIP1.cg14458575.txt")
pdf('./cg14458575_LRRFIP1.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="cg14458575_LRRFIP1",
              xlabel = 'mQTL effect sizes',
              ylabel = 'eQTL effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()



#smr locus
file_eqtl=read.table("./LRRFIP1.txt",header = T)
file_gwas=read.table("../Osteoporosis_GWAS.ma",header = T)
file_meqtl=read.table("./chr2_mqtl.txt",header = T,stringsAsFactors = F)


file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000124831",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
eqtl = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'cis-eQTL \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, face = "plain"),  
        axis.title = element_text(size = 24, face = "bold"))
ggsave("./LRRFIP1.eqtl.zoom.pdf",width = 20,height = 5)

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000124831"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
write.table(file_plot,file = "./tmp.zoom.file.txt",quote = F,row.names = F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
gwas = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'GWAS \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, face = "plain"),  
        axis.title = element_text(size = 24, face = "bold"))
ggsave("./LRRFIP1.gwas.zoom.pdf",width = 20,height = 5)


file_plot=file_meqtl
file_plot=file_plot[file_plot$Probe=="cg14458575",]
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000124831"],]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
file_plot$p=as.numeric(as.character(file_plot$p))
file_plot$BP=as.numeric(as.character(file_plot$BP))
mqtl = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#009988")+theme_bw()+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'cis-mQTL \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, face = "plain"),  
        axis.title = element_text(size = 24, face = "bold"))
ggsave("./LRRFIP1.cg14458575.zoom.pdf",width = 20,height = 5)


pdf(file = 'chr_locus.pdf',width = 20,height = 5)
# 创建空白的染色体图
kp <- plotKaryotype(genome = "hg19", chromosomes = "chr2")

# 添加染色体带

kpAddCytobandsAsLine(kp)

# 添加基因位置框
gene_start <- 54300434 #  起始位置
gene_end <- 54422763   #  结束位置
kpRect(kp, chr = "chr2", x0 = gene_start, x1 = gene_end, y0 = 0, y1 = 0.2, col = NA, border = "red", lwd = 2)

# 添加基因名称
kpText(kp, chr = "chr2", x = (gene_start + gene_end) / 2, y = 0.3, labels = "LRRFIP1", col = "red", cex = 1, pos = 4)


# 保存图表
dev.off()

locus = gwas/eqtl/mqtl
ggsave("./smr_locus.pdf",width = 20,height = 12)







#画gtex-op
SMRData = ReadSMRData("./plot.ENSG00000074603.txt")
pdf('./DPP8_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="DPP8_Osteoporosis",
              xlabel = 'eQTL effect sizes',
              ylabel = 'GWAS effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()

SMRData = ReadSMRData("./ENSG00000105738_op.ENSG00000105738.txt")
pdf('./SIPA1L3_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="SIPA1L3_Osteoporosis",
              xlabel = 'eQTL effect sizes',
              ylabel = 'GWAS effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()

SMRData = ReadSMRData("./ENSG00000168036_op.ENSG00000168036.txt")
pdf('./CTNNB1_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="CTNNB1_Osteoporosis",
              xlabel = 'eQTL effect sizes',
              ylabel = 'GWAS effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()

SMRData = ReadSMRData("./ENSG00000196923_op.ENSG00000196923.txt")
pdf('./PDLIM7_op.pdf',width = 10,height = 10)
SMREffectPlot(data=SMRData,
              trait_name="PDLIM7_Osteoporosis",
              xlabel = 'eQTL effect sizes',
              ylabel = 'GWAS effect sizes',
              cisWindow = 1000,transWindow = 1000)
dev.off()





#gtex smr locusplot
file_eqtl=read.table("./Muscle_Skeletal_fil.txt",header = T)
file_gwas=read.table("../Osteoporosis_GWAS.ma",header = T)



file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000074603",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
DPP8 = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'DPP8 eqtl \n-log10(p)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000105738",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
SIPA1L3 = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'SIPA1L3 eqtl \n-log10(p)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000168036",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
CTNNB1 = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'CTNNB1 eqtl \n-log10(p)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_eqtl
file_plot=file_plot[file_plot$Probe=="ENSG00000196923",]
file_plot$shape=21
file_plot$shape=as.factor(file_plot$shape)
PDLIM7 = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#DDAA33")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_fill_gradient2(low="navy", mid = "lightblue",high="red2",midpoint = 0.4)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'PDLIM7 eqtl \n-log10(p)') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000074603"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
ddp8gwas = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'GWAS \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000105738"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
sipagwas = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'GWAS \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000168036"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
ctnnbgwas = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'GWAS \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_plot=file_gwas
file_plot=file_plot[file_plot$SNP %in% file_eqtl$SNP[file_eqtl$Probe=="ENSG00000196923"],]
file_plot=merge(file_plot,file_eqtl[,c("SNP","BP")],by="SNP",all=F)
file_plot$shape=NA
file_plot$shape[file_plot$SNP=="rs744166"]=23
file_plot$shape[file_plot$SNP!="rs744166"]=21
file_plot$shape=as.factor(file_plot$shape)
pdlgwas = ggplot(file_plot, aes(BP, -log10(p))) +
  geom_point(aes(shape=shape), color = "black", size = 2,fill="#4477AA")+theme_bw()+
  #xlim(4500000,6000000)+
  scale_shape_manual(values = c(21,23))+guides(shape=F)+
  labs(y = 'GWAS \n-log10(p)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



pdf(file = 'PDLIM7_chr_locus.pdf',width = 20,height = 5)
# 创建空白的染色体图
kp <- plotKaryotype(genome = "hg19", chromosomes = "chr5")

# 添加染色体带

kpAddCytobandsAsLine(kp)

# 添加基因位置框
gene_start <- 176910395 #  起始位置
gene_end <- 176924605   #  结束位置
kpRect(kp, chr = "chr5", x0 = gene_start, x1 = gene_end, y0 = 0, y1 = 0.2, col = NA, border = "red", lwd = 2)

# 添加基因名称
kpText(kp, chr = "chr5", x = (gene_start + gene_end) / 2, y = 0.3, labels = "PDLIM7", col = "red", cex = 1, pos = 4)


# 保存图表
dev.off()

locus = pdlgwas/PDLIM7
ggsave("./pdlim7_smr_locus.pdf",width = 20,height = 12)
