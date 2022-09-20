#!/usr/bin/Rscript
rm(list=ls())

## Objective
# apply fgsea package to edgeR data
library(fgsea)
library(tibble)
library(readr)
## edgeR data generated from
# 20220831_microbial_RNAseq_edgeR.R
cecal_LRT_dir <- "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/output/DE_LRT_edgeR_filt_E6nf"
ileal_LRT_dir <- "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/output/DE_LRT_edgeR_filt_E6nf"
gsea_input_dir <- "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/gsea_input"
gsea_output_dir <- "~/mstp/Gordon/experiments/MG00_HWDC_og/HOSVD/DE/gsea_output"
dir.create(gsea_output_dir,showWarnings = FALSE)

#Read CAZYme table
CAZy_fpath = sprintf("%s/CAZyme_annotation_summary.txt",gsea_input_dir)
CAZy = read.csv(CAZy_fpath, sep = '\t')
colnames(CAZy) = c("bacterial_name", "proteinId", "description", "modelnotes", "defline") 
# mcSEED_fpath = sprintf("%s/19isolates_mcseed_pathwaycomplete.csv",gsea_input_dir) #semicolon separated lists of Pht/Functional category
mcSEED_fpath = sprintf("%s/mcSEED_single_pht.csv",gsea_input_dir) #All entries with single phenotype or functional category, (more) duplicate locus tags 
mcSEED_df = read.csv(mcSEED_fpath)


## mcSEED annotation
# mcSEED_anno_dir <- "/scratch/jglab/hwchang/mcSEED_annotation/mcseed_anno_in_table_at_module_level_2"

## edgeR output
# edgeR_dir <- "/scratch/jglab/hwchang/MG10_microbial_RNAseq/from_R/edgeR"

## fgsea 
# output <- "/scratch/jglab/hwchang/MG10_microbial_RNAseq/from_R/fgsea_at_module_level_2"

#=========================================================================================================#
#                                     all microbes in cecal contents                                      #
#=========================================================================================================#

## mapping file
# mapping_file <- "/scratch/jglab/hwchang/MG10_microbial_RNAseq/map_DE_genes_to_mcSEED/mapping_file_for_merging_mod.txt"
# map <- read.csv(mapping_file, header = F)

#Replace map with bacteria_info analog to edgeR analysis 
bacteria_ID <- c("NJCFFJJN", "ONMCJBAG", "OOAPABDJ", "NBCBLOMG", "HIAFFLEM", "ANCJAENF", "NNMNEBCK",
                 "AGCELHME", "BACHIBIB", "EMCKBHPA", "AKLJOBCP", "HMLHAAEP", "CGOBPECO", "BOMJNLPF",
                 "GJBELKAJ", "CGCJBKJN", "JKAMACKM", "LDOIJNDB")
bacteria_name <- c("Pco", "Mmu", "Dlo", "Pst", "Rob", "Bbr", "Bca", "Blu", "Lru", "Dfo", "Eav", "Eco",
                   "Fpr", "Lga4B6", "Rgn", "Rto", "Sga", "Spa")
bacteria_info <- data.frame(bacteria_ID = bacteria_ID,
                            bacteria_name = bacteria_name)
comparison <- "p1Cn2B"
num_of_gene_in_pathway = c()
MIN_PATHWAY_SIZE = 2

# for (b in 1:nrow(map)){ #equivalent to iterating over three char strain abbreviations 
  #b = 1
  # prefix <- map$V1[b]
for(i in 1:nrow(bacteria_info)){
# for(i in 1:2){ #Mmu
# for(i in 5:6){ #Rob
# for(i in 7:8){ #Bca
# for(i in 12:13){ #Eco
  #Skip strains if no reads mapped or no reads survived filtering 
  prefix <- bacteria_info$bacteria_name[i]
  locus_prefix <- bacteria_info$bacteria_ID[i]
  if (prefix %in% c("Pco","Pst","Dfo","Blu","Fpr")) next 
  ## input files
  edgeR.f <- paste(cecal_LRT_dir,"/",prefix, "_cecal_", comparison, "_all.csv", sep = "") 
  #mcSEED originally was divided into files by strain which are equivalent to subsets 
  #of the 19 isolates table (+- weird old mcSEED module structure)
  # mcSEED <- read.csv(mcSEED_f, sep = "\t") 
  mcSEED <- mcSEED_df[which(grepl(locus_prefix,mcSEED_df$Locus.tag,fixed=TRUE)),]
  
  
  df = as.data.frame(table(mcSEED$Phenotype))
  #Note this filters based on pathway size in genome, not based on 
  #pathway size among filtered expressed genes
  gsea_path_df <- df[which(df$Freq >= MIN_PATHWAY_SIZE), ] 
  
  num_of_gene_in_pathway <- c(num_of_gene_in_pathway, df$Freq)
  ## build pathway list
  key = gsea_path_df$Var1
  path_list = vector(mode = "list", length = length(key))
  names(path_list) = key
  
  

  if (dim(gsea_path_df)[1][1] > 0){
    for(i in 1:nrow(gsea_path_df)){
      pathway_name = gsea_path_df$Var1[i]
      # path_list[[i]] <- mcSEED$gene_id[which(mcSEED$subsystem == pathway_name)]
      path_list[[i]] <- mcSEED$Locus.tag[which(mcSEED$Phenotype == pathway_name)]
    }
    

    ## prepare rank information
    # setwd(edgeR_dir)
    # edgeR.f <- paste(prefix, "_cecal_", comparison, ".csv", sep = "")
    rank_raw = read.csv(edgeR.f, row.names = 1)

    ## build gene set
    gene_list <- unique(mcSEED$Locus.tag)

    mcSEED_gene_loc = which(row.names(rank_raw) %in% gene_list)
    mcSEED_rank <- rank_raw$logFC[mcSEED_gene_loc]
    names(mcSEED_rank) <- row.names(rank_raw)[mcSEED_gene_loc]
    
    ## fgsea
    if (length(unique(mcSEED_rank)) > 0){
      if (length(unique(mcSEED_rank > 0)) > 1){   
        fgseaRes <- fgsea(pathways = path_list, 
                          stats    = mcSEED_rank,
                          eps = 0.0,
                          minSize  = MIN_PATHWAY_SIZE,
                          maxSize  = 500)     
      }else if(length(unique(mcSEED_rank > 0)) == 1 && unique(mcSEED_rank > 0) == TRUE){
        fgseaRes <- fgsea(pathways = path_list, 
                      stats    = mcSEED_rank,
                      eps = 0.0,
                      minSize  = MIN_PATHWAY_SIZE,
                      maxSize  = 500,
                      scoreType = "pos")
      }else if(length(unique(mcSEED_rank > 0)) == 1 && unique(mcSEED_rank > 0) == FALSE){
        fgseaRes <- fgsea(pathways = path_list, 
                        stats    = mcSEED_rank,
                        eps = 0.0,
                        minSize  = MIN_PATHWAY_SIZE,
                        maxSize  = 500,
                        scoreType = "neg")
      }else{"ERROR!!!"} # (unique(mcSEED_rank > 0) == TRUE)
    
      fgseaRes <- fgseaRes[order(padj), ]
      fgseaRes_sign <- fgseaRes[fgseaRes$padj < 0.05, ]
    
      if (nrow(fgseaRes_sign) > 0){
        fgseaRes_sign_res <- fgseaRes_sign[, 1:7]
      
        #setwd(output)
        report_f <- paste(prefix, "_cecal_", comparison, "_fgsea.csv", sep = "")
        #write.csv(as.data.frame(fgseaRes_sign_res), report_f)
      }else{
        print(paste(prefix, comparison, "-- no significant enriched pathways"))
      } # if (nrow(fgseaRes_sign) > 0)
    
    }else{
      print(paste(prefix, comparison, "-- no gene match targeted mcSEED pathways"))      
    
    } # if (unique(mcSEED_rank) > 0)
  }else{print(paste(prefix, comparison, "-- no mcSEED pathways that has more than 9 genes")) } # (dim(gsea_path_df)[1][1] > 0)
  cecal_fgsea_outdir = paste(gsea_output_dir,"/cecal",sep="")
  dir.create(cecal_fgsea_outdir,showWarnings = FALSE)
  fgsea_fpath = paste(cecal_fgsea_outdir,"/",prefix,"_cecal_",comparison,"_fgsea.csv",sep="")
  fgsea_sig_fpath = paste(cecal_fgsea_outdir,"/",prefix,"_cecal_",comparison,"_fgsea_sig.csv",sep="")
  fgseaRes.tibble <- as_tibble(fgseaRes)
  # fgseaRes.tibble$leadingEdgePaste <- paste(fgseaRes.tibble$leadingEdge,collapse=",")
  fgseaRes_sign.tibble <- as_tibble(fgseaRes_sign)
  write_csv(fgseaRes.tibble,fgsea_fpath)
  write_csv(fgseaRes_sign.tibble,fgsea_sig_fpath)
}

# 11 warnings
# In preparePathwaysAndStats(pathways, stats, minSize, maxSize,  ... :
# All values in the stats vector are greater than zero and scoreType is "std", maybe you should switch to scoreType = "pos".

# run loop one at a time to find out the source of those warning msg
# b = 1; pass
# b = 2; pass
# b = 3; warning msg when c = 2
# problem solved
# re-run

## 14 warnings
# In if (unique(mcSEED_rank > 0) == TRUE) { ... :
# the condition has length > 1 and only the first element will be used
# problem solved
# re-run

table(num_of_gene_in_pathway)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   27   47   48 
#1602  831  693  543  345  267  219  171  129  123  105   63   45   42   30   15    3   21   18    9    6    6    3    9    3    3

num_of_gene_in_pathway.df <- data.frame(table(num_of_gene_in_pathway))

setwd(output)
write.csv(num_of_gene_in_pathway.df, "num_of_gene_in_pathway_cecal.csv")


#=========================================================================================================#
#                                     all microbes in ileal contents                                      #
#=========================================================================================================#

## mapping file
mapping_file <- "/scratch/jglab/hwchang/MG10_microbial_RNAseq/map_DE_genes_to_mcSEED/mapping_file_for_merging_ileal.txt"
map <- read.csv(mapping_file, header = F)
comparison_list <- c("p1Cn1B", "p1Bn2B", "p1Cn2B")

num_of_gene_in_pathway_ileal <- c()

for (b in 1:nrow(map)){
  prefix <- map$V1[b]
  
  for (c in 1:length(comparison_list)){
    comparison <- comparison_list[c]
    
    ## input files
    f <- paste(prefix, "_ileal_", comparison, ".csv", sep = "") 
    mcSEED_f <- paste(prefix, "_in_table.txt", sep = "")
    
    setwd(mcSEED_anno_dir)
    mcSEED <- read.csv(mcSEED_f, sep = "\t")
    df = as.data.frame(table(mcSEED$subsystem))
    gsea_path_df <- df[which(df$Freq > 9), ]
    num_of_gene_in_pathway_ileal <- c(num_of_gene_in_pathway_ileal, df$Freq)
    
    ## build pathway list
    key = gsea_path_df$Var1
    path_list = vector(mode = "list", length = length(key))
    names(path_list) = key
    
    if (dim(gsea_path_df)[1][1] > 0){    
      for(i in 1:nrow(gsea_path_df)){
        pathway_name = gsea_path_df$Var1[i]
        path_list[[i]] <- mcSEED$gene_id[which(mcSEED$subsystem == pathway_name)]
      }
    
      ## prepare rank information
      setwd(edgeR_dir)
      edgeR.f <- paste(prefix, "_ileal_", comparison, ".csv", sep = "")
      rank_raw = read.csv(edgeR.f, row.names = 1)
    
      ## build gene set
      gene_list <- unique(mcSEED$gene_id)
    
      mcSEED_gene_loc = which(row.names(rank_raw) %in% gene_list)
      mcSEED_rank <- rank_raw$logFC[mcSEED_gene_loc]
      names(mcSEED_rank) <- row.names(rank_raw)[mcSEED_gene_loc]
    
      ## fgsea
      if (length(unique(mcSEED_rank)) > 0){
        if (length(unique(mcSEED_rank > 0)) > 1){   
          fgseaRes <- fgsea(pathways = path_list, 
                          stats    = mcSEED_rank,
                          eps = 0.0,
                          minSize  = 10,
                          maxSize  = 500,
                          nPermSimple = 10000)     
        }else if(length(unique(mcSEED_rank > 0)) == 1 && unique(mcSEED_rank > 0) == TRUE){
          fgseaRes <- fgsea(pathways = path_list, 
                          stats    = mcSEED_rank,
                          eps = 0.0,
                          minSize  = 10,
                          maxSize  = 500,
                          nPermSimple = 10000,
                          scoreType = "pos")
        }else if(length(unique(mcSEED_rank > 0)) == 1 && unique(mcSEED_rank > 0) == FALSE){
          fgseaRes <- fgsea(pathways = path_list, 
                          stats    = mcSEED_rank,
                          eps = 0.0,
                          minSize  = 10,
                          maxSize  = 500,
                          nPermSimple = 10000,
                          scoreType = "neg")
        }else{"ERROR!!!"} # (unique(mcSEED_rank > 0) == TRUE)
      
        fgseaRes <- fgseaRes[order(padj), ]
        fgseaRes_sign <- fgseaRes[fgseaRes$padj < 0.05, ]
      
        if (nrow(fgseaRes_sign) > 0){
          fgseaRes_sign_res <- fgseaRes_sign[, 1:7]
        
          setwd(output)
          report_f <- paste(prefix, "_ileal_", comparison, "_fgsea.csv", sep = "")
          write.csv(as.data.frame(fgseaRes_sign_res), report_f)
        }else{
          print(paste(prefix, comparison, "-- no significant enriched pathways"))
          #print()
        } # if (nrow(fgseaRes_sign) > 0)
      
      }else{
        print(paste(prefix, comparison, "-- no gene match targeted mcSEED pathways"))      
      
      } # if (unique(mcSEED_rank) > 0)
    }else{print(paste(prefix, comparison, "-- no mcSEED pathways that has more than 14 genes")) } # (dim(gsea_path_df)[1][1] > 0)
  }
}

# Warning messages:
#1: In fgseaMultilevel(...) :
#  There were 1 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
#2: In fgseaMultilevel(...) :
#  There were 1 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)

# update code and re-run

# fixed! no completed with warning code

table(num_of_gene_in_pathway_ileal)

num_of_gene_in_pathway_ileal.df <- data.frame(table(num_of_gene_in_pathway_ileal))

setwd(output)
write.csv(num_of_gene_in_pathway_ileal.df, "num_of_gene_in_pathway_ileal.csv")

