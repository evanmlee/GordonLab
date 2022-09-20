#!/usr/bin/Rscript
rm(list=ls())


library(edgeR)
library(ggplot2)

###Edited from original Hao-Wei script
DE_dir = "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE"
input_dir = "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/input/edgeR_filt"
output_dir = "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/output"

setwd(DE_dir)
# count_fpath = sprintf("%s/bpm_filt_count.csv",input_dir)
count_fpath = sprintf("%s/full_merged_count_int.csv",input_dir)
count_data = read.csv(count_fpath)

###End Edited 
dim(count_data)
# [1] 48390   31 #updated

count_data[1:3, 1:2]
#        target_id Pup_1.cecal_contents_53_1C_Pup_1
# 1 ANCJAENF_00001                            31593
# 2 ANCJAENF_00002                           614048
# 3 ANCJAENF_00003                            31593

count_data.data_only = count_data[, 2:ncol(count_data)]
dim(count_data.data_only)
row.names(count_data.data_only) <- count_data$target_id
dim(count_data.data_only)
# [1] 2604   15


## prepare sample info table
source_str = sapply(strsplit(colnames(count_data.data_only), "\\."), "[[", 2)
sample_info.df = data.frame(sample_info = colnames(count_data.data_only),
                            mouse_id = sapply(strsplit(source_str, "_"), "[[", 6),
                            dam_or_pup = sapply(strsplit(source_str, "_"), "[[", 5),
                            treatment = sapply(strsplit(source_str, "_"), "[[", 4),
                            postnatal_day_when_sac = sapply(strsplit(source_str, "_"), "[[", 3),
                            sample_from = sapply(strsplit(source_str, "_"), "[[", 1))
head(sample_info.df)

# sanity check: check to make sure there is no dam
unique(sample_info.df$dam_or_pup)
# [1] "Pup
# pass sanity check


## get CAZyme annotation information                                        ## TODO: need to remove Rto from "CAZyme_annotation_summary.txt"
CAZy_fpath = sprintf("%s/CAZyme_annotation_summary.txt",input_dir)
CAZy = read.csv(CAZy_fpath, sep = '\t')
colnames(CAZy) = c("bacterial_name", "proteinId", "description", "modelnotes", "defline") 


## subset data to only include 1B, 1C, 2B
# sample_info.df.only1B1C2B <- sample_info.df[which(sample_info.df$treatment %in% c("1B", "1C", "2B")), ]
sample_info.df.only1B1C2B <- sample_info.df[which(sample_info.df$treatment %in% c("1C", "2B")), ]
nrow(sample_info.df.only1B1C2B)
# [1] 30
count_data.1B1C2Bdata_only <- count_data.data_only[ ,which(colnames(count_data.data_only) %in% sample_info.df.only1B1C2B$sample_info)]
dim(count_data.1B1C2Bdata_only)
# [1] 48390    30

sample_info.df.only1B1C2B$treatment_sample_from <- paste(sample_info.df.only1B1C2B$treatment, sample_info.df.only1B1C2B$sample_from, sep = "")

design_ready.1B1C2B = data.frame(file = sample_info.df.only1B1C2B$sample_info,
                                 condition = sample_info.df.only1B1C2B$treatment_sample_from)
count_data.1B1C2Bdata_only <- round(count_data.1B1C2Bdata_only, digits = 0)


###########################################################################################################
#                             bacteria info and sample group design matrix                               #
###########################################################################################################

out_dir = "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/output"

bacteria_ID <- c("NJCFFJJN", "ONMCJBAG", "OOAPABDJ", "NBCBLOMG", "HIAFFLEM", "ANCJAENF", "NNMNEBCK",
                 "AGCELHME", "BACHIBIB", "EMCKBHPA", "AKLJOBCP", "HMLHAAEP", "CGOBPECO", "BOMJNLPF",
                 "GJBELKAJ", "CGCJBKJN", "JKAMACKM", "LDOIJNDB")
bacteria_name <- c("Pco", "Mmu", "Dlo", "Pst", "Rob", "Bbr", "Bca", "Blu", "Lru", "Dfo", "Eav", "Eco",
                   "Fpr", "Lga4B6", "Rgn", "Rto", "Sga", "Spa")

design_ready.1B1C2B = data.frame(file = sample_info.df.only1B1C2B$sample_info,
                                 condition = sample_info.df.only1B1C2B$treatment_sample_from)
design_ready.1B1C2B.cecal = design_ready.1B1C2B[which(grepl("cecal", design_ready.1B1C2B$condition) == T), ]
design_ready.1B1C2B.ileal = design_ready.1B1C2B[which(grepl("ileal", design_ready.1B1C2B$condition) == T), ]


bacteria_info <- data.frame(bacteria_ID = bacteria_ID,
                            bacteria_name = bacteria_name,
                            up_p1Cn2B = rep(NA, length(bacteria_ID)),
                            down_p1Cn2B = rep(NA, length(bacteria_ID)))

###########################################################################################################
#                             Filtering by expression - edgeR or external tpm filtered                    #
###########################################################################################################


count_data.1B1C2Bdata_only.cecal <- count_data.1B1C2Bdata_only[which(grepl("cecal", colnames(count_data.1B1C2Bdata_only)) == T)]
dgeFull.cecal <- DGEList(counts = count_data.1B1C2Bdata_only.cecal, group = design_ready.1B1C2B.cecal$condition)  

print(dgeFull.cecal$samples$lib.size)
keep.cecal <- filterByExpr(dgeFull.cecal)
# keep.cecal <- filterByExpr(dgeFull.cecal,min.count=3) #non-default edgeR filtering 
length(keep.cecal[keep.cecal])

count.cecal.keep <- count_data.1B1C2Bdata_only.cecal[keep.cecal,]
edger_filtered_fpath <- sprintf("%s/edgeR_filtered_count_cecal.csv",out_dir)
write.csv(count.cecal.keep, edger_filtered_fpath,quote=FALSE)

#Using tpm-filtered dataset to filter count_data using same criteria as SVD/PCA 
tpm_prev_filt_rlog_fpath = sprintf("%s/tpm3_prev0.25_rlog.csv",input_dir)
tpm_prev_filt_rlog = read.csv(tpm_prev_filt_rlog_fpath,row.names=1)

DATA_FILTER = "edgeR"

if (DATA_FILTER == "tpm") {
  tpm_filt_counts.1C2Bonly <- count_data.1B1C2Bdata_only.cecal[row.names(tpm_prev_filt_rlog),]  
} else if (DATA_FILTER == "edgeR") {
  tpm_filt_counts.1C2Bonly <- count_data.1B1C2Bdata_only.cecal[keep.cecal,]  
}

design_ready.1B1C2B.cecal = design_ready.1B1C2B[which(grepl("cecal", design_ready.1B1C2B$condition) == T), ]
dgeFull.cecal <- DGEList(counts = tpm_filt_counts.1C2Bonly, group = design_ready.1B1C2B.cecal$condition)  
dge.cecal.norm <- calcNormFactors(dgeFull.cecal)

#Whole dataset library size and norm factors 
print(dge.cecal.norm$samples$lib.size)
print(dge.cecal.norm$samples$norm.factors)

#Deprecated - filter to bpm annotated genes before running DE -> bpm_filt.dge; bpm_filt.dge is no longer used
mcseed_fpath = sprintf("%s/19isolates_mcseed_pathwaycomplete.csv",input_dir)
mcseed.df <- read.csv(mcseed_fpath)
bpm_filt.dge <- dge.cecal.norm[row.names(dge.cecal.norm) %in% mcseed.df$Locus.tag,,keep.lib.sizes=TRUE]
#====#

###########################################################################################################
#                             Iterative DGE by bacteria                                                   #
###########################################################################################################

#Read cecal abundance data 
cecal_abundance_fpath = sprintf("%s/cecal_abundance_E6.csv",input_dir)
# cecal_abundance_fpath = sprintf("%s/cecal_abundance_E6_floored.csv",input_dir)
cecal_abundance.df <- read.csv(cecal_abundance_fpath,row.names=1)
dropped_bacteria <- c("Pco","Pst","Dfo","Blu","Fpr")

strain_level_norm = TRUE #TRUE will provide better separation of DE vs not DE within each strain
ABUNDANCE_MODEL = TRUE #Include abundance of each strain as a covariate in the design matrix/ GLM
BPM_FILT = FALSE #BPM_filtering before GLM/DE - FALSE, include unannotated genes
#for GLM fitting 
if (BPM_FILT) {
  dge.cecal.norm <- dge.cecal.norm[row.names(dge.cecal.norm) %in% mcseed.df$Locus.tag,,keep.lib.sizes=TRUE]
} else {
  dge.cecal.norm <- dge.cecal.norm
}

for(i in 1:nrow(bacteria_info)){
  
  #Skip strains if no reads mapped or no reads survived filtering 
  if (bacteria_info$bacteria_name[i] %in% c("Pco","Pst","Dfo","Blu","Fpr")) next 
  
  bact_name <- bacteria_info$bacteria_name[i]
  print(bact_name)
  bacteria_ID <- bacteria_info$bacteria_ID[i]
  
  #Filter down to genes from strain X; recalc norm factors if strain_level_norm
  if (strain_level_norm) {
    dge.cecal.norm.X <- dge.cecal.norm[which(grepl(bacteria_info$bacteria_ID[i], 
                                               row.names(dge.cecal.norm)) == T),,
                                   keep.lib.sizes=FALSE]
    dge.cecal.norm.X <- calcNormFactors(dge.cecal.norm.X) 
  }
  else {
    dge.cecal.norm.X <- dge.cecal.norm[which(grepl(bacteria_info$bacteria_ID[i], 
                                               row.names(dge.cecal.norm)) == T),,
                                   keep.lib.sizes=TRUE]
  }
  #Sanity check output to double check that correct number of genes being used for DE 
  print(length(row.names(dge.cecal.norm.X$counts)))
  #Add cecal abundance to design matrix 
  cecal_abundance.X <- cecal_abundance.df[bact_name]
  design_mat.X <- data.frame(design_ready.1B1C2B.cecal)
  design_mat.X["abundance_X"] <- cecal_abundance.X
  ## GLM estimates of dispersion
  #From Hao-Wei's original comment:
  # Fitting a model in edgeR takes several steps. First, you must fit the common dispersion. Then you need 
  # to fit a trended model (if you do not fit a trend, the default is to use the common dispersion as a trend). 
  # Then you can fit the tagwise dispersion which is a function of this model.
  if (ABUNDANCE_MODEL) {
    design.mat.X.cecal.glm <- model.matrix(~ 0 + design_mat.X$condition + 
                                             design_mat.X$abundance_X)  # 0: no intersect in the lm
    colnames(design.mat.X.cecal.glm) <- c("1Ccecal","2Bcecal","abundance_X")
  }
  else {
    design.mat.X.cecal.glm <- model.matrix(~ 0 + design_mat.X$condition)  # 0: no intersect in the lm
    colnames(design.mat.X.cecal.glm) <- c("1Ccecal","2Bcecal")
  }
  
  dge.cecal.norm.X.GLMCommonDisp <- estimateGLMCommonDisp(dge.cecal.norm.X, design.mat.X.cecal.glm)
  dge.cecal.norm.X.GLMCommonDisp <- estimateGLMTrendedDisp(dge.cecal.norm.X.GLMCommonDisp, design.mat.X.cecal.glm)
  dge.cecal.norm.X.GLMCommonDisp <- estimateGLMTagwiseDisp(dge.cecal.norm.X.GLMCommonDisp, design.mat.X.cecal.glm)
  
  ## GLM testing for differential expression
  fit.X.cecal <- glmFit(dge.cecal.norm.X.GLMCommonDisp, design.mat.X.cecal.glm)
  if (ABUNDANCE_MODEL) {
    lrt12.X.cecal <- glmLRT(fit.X.cecal,contrast=c(1,-1,0))
  } else {
    lrt12.X.cecal <- glmLRT(fit.X.cecal,contrast=c(1,-1))
  }
  
  topTags(lrt12.X.cecal, n=10, adjust.method="BH", sort.by="PValue", p.value=1)
  lrt12.X.cecal.out <- lrt12.X.cecal$table
  # #Bpm filter results
  # lrt12.X.cecal.out <- lrt12.X.cecal.out[row.names(lrt12.X.cecal.out) %in% mcseed.df$Locus.tag,]
  
  #Some different FDR correction options - eventually settled on FDR for all genes within each strain 
  #i.e. instead of FDR for all genes across all strains
  lrt12.X.cecal.out$FDR <- p.adjust(lrt12.X.cecal.out$PValue, method = 'BH', n = nrow(dge.cecal.norm.X))
  # lrt12.X.cecal.out$FDR <- p.adjust(lrt12.X.cecal.out$PValue, method = 'BH', n = nrow(dge.cecal.norm))
  # lrt12.X.cecal.out$FDR <- p.adjust(lrt12.X.cecal.out$PValue, method = 'bonferroni', n = nrow(tpm_filt_counts.1C2Bonly))
  lrt12.X.cecal.out <- lrt12.X.cecal.out[order(lrt12.X.cecal.out$PValue), ]
  
  #BPM filtering and BPM and FDR significant filtering 
  lrt12.X.cecal.bpm <- lrt12.X.cecal.out[row.names(lrt12.X.cecal.out) %in% mcseed.df$Locus.tag,]
  lrt12.X.cecal.sig <- lrt12.X.cecal.bpm[which(lrt12.X.cecal.bpm$FDR < 0.05), ]
  
  lrt12.X.cecal.sig$direction <- ifelse(lrt12.X.cecal.sig$logFC > 0, "pos", "neg")
  lrt12.X.cecal.sig_res <- as.data.frame(table(lrt12.X.cecal.sig$direction))
  
  lrt12.pos <- ifelse((c("pos") %in% lrt12.X.cecal.sig_res$Var1), lrt12.X.cecal.sig_res$Freq[which(lrt12.X.cecal.sig_res$Var1 == "pos")], 0)
  lrt12.neg <- ifelse((c("neg") %in% lrt12.X.cecal.sig_res$Var1), lrt12.X.cecal.sig_res$Freq[which(lrt12.X.cecal.sig_res$Var1 == "neg")], 0)
  bacteria_info$up_p1Cn1B[i] = lrt12.pos
  bacteria_info$down_p1Cn1B[i] = lrt12.neg
  
  DE_outdir = sprintf("%s/DE_LRT_cecal",out_dir)
  dir.create(DE_outdir,showWarnings = FALSE)
  all_file_name = paste(bacteria_info$bacteria_name[i], "_cecal_p1Cn2B_all.csv", sep = "")
  bpm_file_name = paste(bacteria_info$bacteria_name[i], "_cecal_p1Cn2B_bpm.csv", sep = "")
  sig_file_name = paste(bacteria_info$bacteria_name[i], "_cecal_p1Cn2B_bpm_sig.csv", sep = "")
  # strain_out_dir = sprintf("%s/DE/%s",out_dir,bacteria_name)
  # dir.create(strain_out_dir)
  lrt_all_fpath = sprintf("%s/%s",DE_outdir,all_file_name)
  write.csv(lrt12.X.cecal.out,lrt_all_fpath,quote=FALSE)
  lrt_bpm_fpath = sprintf("%s/%s",DE_outdir,bpm_file_name)
  write.csv(lrt12.X.cecal.bpm,lrt_bpm_fpath,quote=FALSE)
  lrt_sig_fpath = sprintf("%s/%s",DE_outdir,sig_file_name)
  write.csv(lrt12.X.cecal.sig,lrt_sig_fpath,quote=FALSE)
}

###########################################################################################################
###                                   ILEAL CONTENTS MICROBIAL RNASEQ                                  ####
###########################################################################################################

###########################################################################################################
#                             Filtering by expression - edgeR or external tpm filtered                    #
###########################################################################################################


count_data.1B1C2Bdata_only.ileal <- count_data.1B1C2Bdata_only[which(grepl("ileal", colnames(count_data.1B1C2Bdata_only)) == T)]
dgeFull.ileal <- DGEList(counts = count_data.1B1C2Bdata_only.ileal, group = design_ready.1B1C2B.ileal$condition)  

print(dgeFull.ileal$samples$lib.size)
keep.ileal <- filterByExpr(dgeFull.ileal)
# keep.ileal <- filterByExpr(dgeFull.ileal,min.count=3) #non-default edgeR filtering 
length(keep.ileal[keep.ileal])

count.ileal.keep <- count_data.1B1C2Bdata_only.ileal[keep.ileal,]
edger_filtered_fpath <- sprintf("%s/edgeR_filtered_count_ileal.csv",out_dir)
write.csv(count.ileal.keep, edger_filtered_fpath,quote=FALSE)

#Using tpm-filtered dataset to filter count_data using same criteria as SVD/PCA 
tpm_prev_filt_rlog_fpath = sprintf("%s/tpm3_prev0.25_rlog.csv",input_dir)
tpm_prev_filt_rlog = read.csv(tpm_prev_filt_rlog_fpath,row.names=1)

DATA_FILTER = "edgeR"

if (DATA_FILTER == "tpm") {
  tpm_filt_counts.1C2Bonly <- count_data.1B1C2Bdata_only.ileal[row.names(tpm_prev_filt_rlog),]  
} else if (DATA_FILTER == "edgeR") {
  tpm_filt_counts.1C2Bonly <- count_data.1B1C2Bdata_only.ileal[keep.ileal,]  
}

design_ready.1B1C2B.ileal = design_ready.1B1C2B[which(grepl("ileal", design_ready.1B1C2B$condition) == T), ]
dgeFull.ileal <- DGEList(counts = tpm_filt_counts.1C2Bonly, group = design_ready.1B1C2B.ileal$condition)  
dge.ileal.norm <- calcNormFactors(dgeFull.ileal)

#Whole dataset library size and norm factors 
print(dge.ileal.norm$samples$lib.size)
print(dge.ileal.norm$samples$norm.factors)

#Deprecated - filter to bpm annotated genes before running DE -> bpm_filt.dge; bpm_filt.dge is no longer used
mcseed_fpath = sprintf("%s/19isolates_mcseed_pathwaycomplete.csv",input_dir)
mcseed.df <- read.csv(mcseed_fpath)
bpm_filt.dge <- dge.ileal.norm[row.names(dge.ileal.norm) %in% mcseed.df$Locus.tag,,keep.lib.sizes=TRUE]
#====#

###########################################################################################################
#                             Iterative DGE by bacteria                                                   #
###########################################################################################################

#Read ileal abundance data 
ileal_abundance_fpath = sprintf("%s/ileal_abundance_E6.csv",input_dir)
# ileal_abundance_fpath = sprintf("%s/ileal_abundance_E6_floored.csv",input_dir)
ileal_abundance.df <- read.csv(ileal_abundance_fpath,row.names=1)
dropped_bacteria <- c("Pco","Pst","Dfo","Blu","Fpr")

strain_level_norm = TRUE #TRUE will provide better separation of DE vs not DE within each strain
ABUNDANCE_MODEL = TRUE #Include abundance of each strain as a covariate in the design matrix/ GLM
BPM_FILT = FALSE #BPM_filtering before GLM/DE - FALSE, include unannotated genes
#for GLM fitting 
if (BPM_FILT) {
  dge.ileal.norm <- dge.ileal.norm[row.names(dge.ileal.norm) %in% mcseed.df$Locus.tag,,keep.lib.sizes=TRUE]
} else {
  dge.ileal.norm <- dge.ileal.norm
}

for(i in 1:nrow(bacteria_info)){
  
  #Skip strains if no reads mapped or no reads survived filtering 
  if (bacteria_info$bacteria_name[i] %in% c("Pco","Pst","Dfo","Blu","Fpr")) next 
  
  bact_name <- bacteria_info$bacteria_name[i]
  print(bact_name)
  bacteria_ID <- bacteria_info$bacteria_ID[i]
  
  #Filter down to genes from strain X; recalc norm factors if strain_level_norm
  if (strain_level_norm) {
    dge.ileal.norm.X <- dge.ileal.norm[which(grepl(bacteria_info$bacteria_ID[i], 
                                                   row.names(dge.ileal.norm)) == T),,
                                       keep.lib.sizes=FALSE]
    dge.ileal.norm.X <- calcNormFactors(dge.ileal.norm.X) 
  }
  else {
    dge.ileal.norm.X <- dge.ileal.norm[which(grepl(bacteria_info$bacteria_ID[i], 
                                                   row.names(dge.ileal.norm)) == T),,
                                       keep.lib.sizes=TRUE]
  }
  #Sanity check output to double check that correct number of genes being used for DE 
  print(length(row.names(dge.ileal.norm.X$counts)))
  #Add ileal abundance to design matrix 
  ileal_abundance.X <- ileal_abundance.df[bact_name]
  design_mat.X <- data.frame(design_ready.1B1C2B.ileal)
  design_mat.X["abundance_X"] <- ileal_abundance.X
  ## GLM estimates of dispersion
  #From Hao-Wei's original comment:
  # Fitting a model in edgeR takes several steps. First, you must fit the common dispersion. Then you need 
  # to fit a trended model (if you do not fit a trend, the default is to use the common dispersion as a trend). 
  # Then you can fit the tagwise dispersion which is a function of this model.
  if (ABUNDANCE_MODEL) {
    design.mat.X.ileal.glm <- model.matrix(~ 0 + design_mat.X$condition + 
                                             design_mat.X$abundance_X)  # 0: no intersect in the lm
    colnames(design.mat.X.ileal.glm) <- c("1Cileal","2Bileal","abundance_X")
  }
  else {
    design.mat.X.ileal.glm <- model.matrix(~ 0 + design_mat.X$condition)  # 0: no intersect in the lm
    colnames(design.mat.X.ileal.glm) <- c("1Cileal","2Bileal")
  }
  
  dge.ileal.norm.X.GLMCommonDisp <- estimateGLMCommonDisp(dge.ileal.norm.X, design.mat.X.ileal.glm)
  dge.ileal.norm.X.GLMCommonDisp <- estimateGLMTrendedDisp(dge.ileal.norm.X.GLMCommonDisp, design.mat.X.ileal.glm)
  dge.ileal.norm.X.GLMCommonDisp <- estimateGLMTagwiseDisp(dge.ileal.norm.X.GLMCommonDisp, design.mat.X.ileal.glm)
  
  ## GLM testing for differential expression
  fit.X.ileal <- glmFit(dge.ileal.norm.X.GLMCommonDisp, design.mat.X.ileal.glm)
  if (ABUNDANCE_MODEL) {
    lrt12.X.ileal <- glmLRT(fit.X.ileal,contrast=c(1,-1,0))
  } else {
    lrt12.X.ileal <- glmLRT(fit.X.ileal,contrast=c(1,-1))
  }
  
  topTags(lrt12.X.ileal, n=10, adjust.method="BH", sort.by="PValue", p.value=1)
  lrt12.X.ileal.out <- lrt12.X.ileal$table
  # #Bpm filter results
  # lrt12.X.ileal.out <- lrt12.X.ileal.out[row.names(lrt12.X.ileal.out) %in% mcseed.df$Locus.tag,]
  
  #Some different FDR correction options - eventually settled on FDR for all genes within each strain 
  #i.e. instead of FDR for all genes across all strains
  lrt12.X.ileal.out$FDR <- p.adjust(lrt12.X.ileal.out$PValue, method = 'BH', n = nrow(dge.ileal.norm.X))
  # lrt12.X.ileal.out$FDR <- p.adjust(lrt12.X.ileal.out$PValue, method = 'BH', n = nrow(dge.ileal.norm))
  # lrt12.X.ileal.out$FDR <- p.adjust(lrt12.X.ileal.out$PValue, method = 'bonferroni', n = nrow(tpm_filt_counts.1C2Bonly))
  lrt12.X.ileal.out <- lrt12.X.ileal.out[order(lrt12.X.ileal.out$PValue), ]
  
  #BPM filtering and BPM and FDR significant filtering 
  lrt12.X.ileal.bpm <- lrt12.X.ileal.out[row.names(lrt12.X.ileal.out) %in% mcseed.df$Locus.tag,]
  lrt12.X.ileal.sig <- lrt12.X.ileal.bpm[which(lrt12.X.ileal.bpm$FDR < 0.05), ]
  
  lrt12.X.ileal.sig$direction <- ifelse(lrt12.X.ileal.sig$logFC > 0, "pos", "neg")
  lrt12.X.ileal.sig_res <- as.data.frame(table(lrt12.X.ileal.sig$direction))
  
  lrt12.pos <- ifelse((c("pos") %in% lrt12.X.ileal.sig_res$Var1), lrt12.X.ileal.sig_res$Freq[which(lrt12.X.ileal.sig_res$Var1 == "pos")], 0)
  lrt12.neg <- ifelse((c("neg") %in% lrt12.X.ileal.sig_res$Var1), lrt12.X.ileal.sig_res$Freq[which(lrt12.X.ileal.sig_res$Var1 == "neg")], 0)
  bacteria_info$up_p1Cn1B[i] = lrt12.pos
  bacteria_info$down_p1Cn1B[i] = lrt12.neg
  
  DE_outdir = sprintf("%s/DE_LRT_ileal",out_dir)
  dir.create(DE_outdir,showWarnings = FALSE)
  all_file_name = paste(bacteria_info$bacteria_name[i], "_ileal_p1Cn2B_all.csv", sep = "")
  bpm_file_name = paste(bacteria_info$bacteria_name[i], "_ileal_p1Cn2B_bpm.csv", sep = "")
  sig_file_name = paste(bacteria_info$bacteria_name[i], "_ileal_p1Cn2B_bpm_sig.csv", sep = "")
  # strain_out_dir = sprintf("%s/DE/%s",out_dir,bacteria_name)
  # dir.create(strain_out_dir)
  lrt_all_fpath = sprintf("%s/%s",DE_outdir,all_file_name)
  write.csv(lrt12.X.ileal.out,lrt_all_fpath,quote=FALSE)
  lrt_bpm_fpath = sprintf("%s/%s",DE_outdir,bpm_file_name)
  write.csv(lrt12.X.ileal.bpm,lrt_bpm_fpath,quote=FALSE)
  lrt_sig_fpath = sprintf("%s/%s",DE_outdir,sig_file_name)
  write.csv(lrt12.X.ileal.sig,lrt_sig_fpath,quote=FALSE)
}


#=========================================================================================================#
#                                interate through all comparsion - cecal                                  #
#=========================================================================================================#
# 
# mcSEED_dir <- "/scratch/jglab/hwchang/MG10_microbial_RNAseq/map_DE_genes_to_mcSEED/add_mcseed_cecal"
# setwd(mcSEED_dir)
# 
# head(bacteria_info)
# #  bacteria_ID bacteria_name up_p1Cn1B down_p1Cn1B up_p1Bn2B down_p1Bn2B up_p1Cn2B down_p1Cn2B
# #1    NJCFFJJN           Pco        11          20       107           3       105           3
# #2    ONMCJBAG           Mmu        54          16         2          35        22          14
# #3    OOAPABDJ           Dlo      2190          25         0        2182        15          26
# #4    NBCBLOMG           Pst        32          49         0           3         0           3
# #5  ROSSTS7063           Rob        27          15         0           0        26          15
# #6    ANCJAENF           Bbr         0           6         8          14        34          49
# 
# header_list <- c("logFC", "logCPM", "LR", "PValue", 
#                  "FDR", "prokka_gene_id", "RAST_id", 
#                  "percent_id", "evalue", "functional_role", 
#                  "gene_symbol", "subsystem_1", "subsystem_2", 
#                  "subsystem_3", "subsystem_4", "from_X_bacteria")
# 
# bacteria_info.minus <- bacteria_info[-which(bacteria_info$bacteria_name == "Fpr"), ]
# 
# for (p in 1:nrow(bacteria_info.minus)){
#   #p = 13
#   
#   ##==1==## p1Cn1B
#   print(bacteria_info.minus$bacteria_name[p])
#   file <- paste(bacteria_info.minus$bacteria_name[p], "_cecal_p1Cn1B_w_mcseed.txt", sep = "")
#   w_mcseed <- read.csv(file, header = F, row.names = 1, sep = "\t")
#   colnames(w_mcseed) <- header_list
#   w_mcseed.c <- w_mcseed[, -which(colnames(w_mcseed) %in% c("prokka_gene_id", "best_mcseed_hit", "percent_id", "evalue", "from_X_bacteria"))]
#   
#   w_mcseed.c.sign <- w_mcseed.c[which(w_mcseed.c$FDR < 0.05), ]
#   bacteria_info.minus$up_p1Cn1B[p] <- nrow(w_mcseed.c.sign[which(w_mcseed.c.sign$logFC > 1) ,])
#   bacteria_info.minus$down_p1Cn1B[p] <- nrow(w_mcseed.c.sign[which(w_mcseed.c.sign$logFC < 1) ,])
#   
#   ##==2==## p1Bn2B
#   file.p1Bn2B <- paste(bacteria_info.minus$bacteria_name[p], "_cecal_p1Bn2B_w_mcseed.txt", sep = "")
#   w_mcseed.p1Bn2B <- read.csv(file.p1Bn2B, header = F, row.names = 1, sep = "\t")
#   colnames(w_mcseed.p1Bn2B) <- header_list
#   w_mcseed.p1Bn2B.c <- w_mcseed.p1Bn2B[, -which(colnames(w_mcseed.p1Bn2B) %in% c("prokka_gene_id", "best_mcseed_hit", "percent_id", "evalue", "from_X_bacteria"))]
#   
#   w_mcseed.p1Bn2B.c.sign <- w_mcseed.p1Bn2B.c[which(w_mcseed.p1Bn2B.c$FDR < 0.05), ]
#   bacteria_info.minus$up_p1Bn2B[p] <- nrow(w_mcseed.p1Bn2B.c.sign[which(w_mcseed.p1Bn2B.c.sign$logFC > 1) ,])
#   bacteria_info.minus$down_p1Bn2B[p] <- nrow(w_mcseed.p1Bn2B.c.sign[which(w_mcseed.p1Bn2B.c.sign$logFC < 1) ,])
#   
#   ##==3==## p1Cn2B
#   file.p1Cn2B <- paste(bacteria_info.minus$bacteria_name[p], "_cecal_p1Cn2B_w_mcseed.txt", sep = "")
#   w_mcseed.p1Cn2B <- read.csv(file.p1Cn2B, header = F, row.names = 1, sep = "\t")
#   colnames(w_mcseed.p1Cn2B) <- header_list
#   w_mcseed.p1Cn2B.c <- w_mcseed.p1Cn2B[, -which(colnames(w_mcseed.p1Cn2B) %in% c("prokka_gene_id", "best_mcseed_hit", "percent_id", "evalue", "from_X_bacteria"))]
#   
#   w_mcseed.p1Cn2B.c.sign <- w_mcseed.p1Cn2B.c[which(w_mcseed.p1Cn2B.c$FDR < 0.05), ]
#   bacteria_info.minus$up_p1Cn2B[p] <- nrow(w_mcseed.p1Cn2B.c.sign[which(w_mcseed.p1Cn2B.c.sign$logFC > 1) ,])
#   bacteria_info.minus$down_p1Cn2B[p] <- nrow(w_mcseed.p1Cn2B.c.sign[which(w_mcseed.p1Cn2B.c.sign$logFC < 1) ,])
# }
# 
# setwd(out_dir)
# write.csv(bacteria_info.minus, "20210602_cecal_mcSEED_DE_gene_summary.csv")
# 
