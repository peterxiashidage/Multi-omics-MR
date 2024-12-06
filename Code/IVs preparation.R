library(tidyverse)
library(data.table)
library(ivpack)
library(meta)
library(devtools)
library(pacman)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(phenoscanner)
library(LDlinkR)
library(mr.raps)
library(MRPRESSO)
library(extrafont)
library(anchors)
library(plinkbinr)
library(friendly2MR)

# IVs preparation function
MR_prep <- function(exp, out) {
  
  # Extracting instruments
  exp_dat <- read_exposure_data(
    filename = paste0('./eqtl_split/',exp,sep = ''),
    clump = FALSE, #是否去除连锁不平衡
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b", #效应值大小
    se_col = "SE", #beta值的标准差
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "p",
    min_pval = 1e-200, 
    log_pval = FALSE, 
    chr_col = "Chr",
    eaf_col = 'eaf',
    pos_col = "BP"
  )
  
  # Clumping instruments
  exp_dat <- exp_dat %>% 
    rename(
      rsid = SNP,
      pval = pval.exposure
    )
  
  #本地clump
  exp_dat_clumped <- ld_clump( #ld_clump是ieugwasr包里面利用PLINK方法去除连锁不平衡的函数
    dat = exp_dat, #输入暴露矩阵
    clump_kb = 10000,  #去除的窗口大小
    clump_r2 = 0.001,  #两个SNP之间的关联性，越接近1越强
    clump_p = 5e-8, #p值
    plink_bin = get_plink_exe(), #这个plinkbinr包也要自己安装一下
    bfile = '../../bfile/g1000_eur' #这几个.bed .bim .fam文件一定要解压出来和Rproject在同一文件夹
  )
  
  exp_dat_clumped <- exp_dat_clumped %>% 
    rename(
      SNP = rsid,
      pval.exposure = pval
    )
  
  
  # Printing number of IVs for exposure
  print(paste0("Number of IVs: ", as.character(length(exp_dat_clumped$SNP))))
  
  # Extracting instruments from outcome GWAS
  out_dat <- read_outcome_data(
    filename = out, 
    snps = exp_dat_clumped$SNP,
    sep = "\t", 
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "P",
    eaf_col = 'A1FREQ',
    min_pval = 1e-200, 
    log_pval = FALSE, 
    chr_col = "CHR",
    pos_col = "BP"
  )
  
  
  # Identifying & printing exposure instruments missing from outcome GWAS
  missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
  print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
  print("List of IVs missing from outcome GWAS:")
  for (i in 1:length(missing_IVs)) {
    print(paste0(missing_IVs[i]))
  }
  
  # Replacing missing instruments from outcome GWAS with proxies
  if(length(missing_IVs) == 0) {
    
    print("All exposure IVs found in outcome GWAS.")
    
  } else {
    
    print("Some exposure IVs missing from outcome GWAS.")
    out_full <- fread(paste0(deparse(substitute(out)), "_outcome.txt", sep = "")) #读取文件，类似于read.table()
    
    for (i in 1:length(missing_IVs)) {
      
      proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = "54ef0a1dd56f", file = FALSE) 
      # LDproxy：它用于在人类基因组中查找与给定变异存在相关性的代理变异。这些代理变异是与给定变异在连锁不平衡（LD）中相关的变异。
      # pop:人口；token：需要用户自己注册
      proxies <- proxies[proxies$R2 > 0.8, ]
      proxy_present = FALSE
      
      if(length(proxies$RS_Number) == 0){
        
        print(paste0("No proxy SNP available for ", missing_IVs[i]))
        
      } else {
        
        for (j in 1:length(proxies$RS_Number)) {
          
          proxy_present <- proxies$RS_Number[j] %in% out_full$SNP
          
          if (proxy_present) {
            proxy_SNP = proxies$RS_Number[j]
            proxy_SNP_allele_1 = str_sub(proxies$Alleles[j], 2, 2)
            proxy_SNP_allele_2 = str_sub(proxies$Alleles[j], 4, 4)
            original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
            original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
            break
          }
        }
      }
      
      if(proxy_present == TRUE) {
        print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
        proxy_row <- out_dat[1, ]
        proxy_row$SNP = missing_IVs[i]
        proxy_row$beta.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
        proxy_row$se.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
        if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
        if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
        if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
        if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
        proxy_row$pval.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "P"])
        proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N"])
        if("N_case" %in% colnames(out_full)) proxy_row$ncase.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_case"])
        if("N_control" %in% colnames(out_full))proxy_row$ncontrol.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_control"])
        proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "chr.exposure"])
        proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "pos.exposure"])
        if("AF1" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "AF1"])
        out_dat <- rbind(out_dat, proxy_row)
      }
      
      if(proxy_present == FALSE) {
        print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
      }
    }
    
  }
  
  # Harmonising exposure and outcome datasets
  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped, 
    outcome_dat = out_dat, 
    action = 2
  )
  
  write.table(dat,file = paste0('./eqtl_split/','f_',exp,sep = ''),quote = F,row.names = F)
}



#切分文件
dir.create('./eqtl_split')
data <- read.table("blood_sig_exposure.txt", header = TRUE, sep = "\t")
split_data <- split(data, data$Probe)

# 使用 lapply 来对分割后的每一部分应用 write.csv
lapply(names(split_data), function(probe) {
  write.table(split_data[[probe]], file = paste0("./eqtl_split/", probe, "_exposure.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
})


gene = list.files('./eqtl_split/')
for (i in gene) {
  exp = i
  out = './OP_outcome.txt'
  MR_prep(exp,out)
}
