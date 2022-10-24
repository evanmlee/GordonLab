import pandas as pd
import numpy as np
import os
import re 
import warnings
"""Utility functions for loading standard data tables and expression data for microbial RNAseq"""

# global STRAIN_ABBREVS, CAZYME_STRAIN_ABBREVS, STRAIN_TAXA_NAMES, MCSEED_ISOLATE_NAMES,MCSEED_TAXA_NAMES,STRAIN_ABBREVS_TO_ISOLATE_NAMES
STRAIN_ABBREVS = ["Bbr","Bca","Bli2D9","Blu","Rob","Dfo","Dlo","Eav","Eco","FprB","Lga4B6","Lru","Mmu","Pco","Pst",
              "Rgn","Rto","Sga","Spa"]
CAZYME_STRAIN_ABBREVS = ["Bbr","Bca","Bli2D9","Bli463","Blu","Bwe","Dfo","Dlo","Eav","Eco","FprB","Lga4B6","Mmu","Pco","Pst",
                  "Rgn","Rob","Rto","Sga","Spa"]
STRAIN_TAXA_NAMES = [] 
MCSEED_ISOLATE_NAMES = []
MCSEED_TAXA_NAMES = []
STRAIN_ABBREVS_TO_ISOLATE_NAMES = {}

############################ Counts Loading and Simple Filtering ############################
def load_kallisto_pseudocounts(counts_fpath,col_stripping="",col_filtering_pats=[],dropped_locus_tags=[]):
    """Load kallisto_pseudocounts table from counts_fpath. 
    @param counts_fpath str: file path to counts_table. Can be .tsv or .csv, which will use separators r'\s+' or default 
    read_csv options. The column name for locus tags should be 'target_id'. 
    @param col_stripping str, {'est_counts', 'tpm'}: For use with raw kallisto outputs only. 
    If provided, will strip out corresponding column tags (such as .est_counts) from column names (default "")
    @param col_filtering_pats, list or array-like: if provided, filters columns in returned counts_df to only (default []) 
    @param dropped_locus_tags (list or array-like): If provided, removes rows corresponding to locus tags containing entries in 
    dropped_locus_tags. This can be desirable if reads were mapped to wrong genomes which should be excluded from analysis 
    (default []). 
    @return counts_df: DataFrame of index of locus tags, columns of sample names, and values of kallisto data 
    (either pseudocounts or tpm). 

    Usages:
    For raw kallisto outputs: provide col_stripping (i.e. one of {'est_counts','tpm'}), col_filtering_pats and 
    dropped_locus tags if desired. 
    For already processed kallisto_outputs: Providing col_stripping when the column names don't contain the corresponding 
    tags ('.est_counts' or .'tpm') will cause errors. 
    """
    COUNTS_STRIP_PAT = r'(.*)\.est_counts'
    TPM_STRIP_PAT = r'(.*)\.tpm'

    #Flexible whitespace parsing - tpm first column is tab separated then space separated 
    if re.search(r'\.csv',counts_fpath):
        counts_df = pd.read_csv(counts_fpath)
    elif re.search(r'\.tsv',counts_fpath):
        counts_df = pd.read_csv(counts_fpath,sep=r"\s+")
    counts_df.set_index("target_id",inplace=True)

    #Remove .tpm and .est_counts tags from column names in tpm and count dfs 
    if col_stripping:
        if col_stripping == "est_counts":
            counts_df.columns = counts_df.columns.str.extract(COUNTS_STRIP_PAT,expand=False)
        elif col_stripping == "tpm":
            counts_df.columns = counts_df.columns.str.extract(TPM_STRIP_PAT,expand=False)
    #filter columns down to columns which contain entries in col_filtering_pats
    if len(col_filtering_pats)>0:
        for pat in col_filtering_pats:
            counts_df = counts_df.loc[:,counts_df.columns.str.contains(pat)]
    #Drop all loci matching entries in dropped_locus_tags (i.e. due to errors mapping to non-target genomes)
    if len(dropped_locus_tags) > 0:
        counts_df = filter_out_locus_tags(counts_df,dropped_locus_tags)
    return counts_df

def load_multiple_kallisto_pseudocounts(counts_fpaths=[],col_stripping="",col_filtering_pats=[],dropped_locus_tags=[]):
    return pd.concat([load_kallisto_pseudocounts(fpath,col_stripping,col_filtering_pats,dropped_locus_tags) for fpath in counts_fpaths])

def filter_out_locus_tags(counts_df,dropped_locus_tags=[]):
    dropped_loci_re_pat = "|".join(dropped_locus_tags)
    return counts_df.loc[~counts_df.index.str.contains(dropped_loci_re_pat)]

def int_round_counts(counts_df,how="ceil"):
    counts_int = counts_df.copy()
    if how == "ceil":
        counts_int.iloc[:,:] = np.ceil(counts_int.iloc[:,:]).astype(int)
    elif how == "floor":
        counts_int.iloc[:,:] = np.ceil(counts_int.iloc[:,:]).astype(int)
    else:
        raise ValueError("Unknown option for how; use {'ceil','floor'}.")
    return counts_int

############################ Abundance Table Loading and Transformations ############################
def load_cecal_abundance_table(cecal_abundance_fpath,abundance_strain_abbrevs=[]):
    cecal_abundance_df = pd.read_csv(cecal_abundance_path,sep=",",index_col="MouseID")
    if len(abundance_strain_abbrevs) > 0: 
        pass #TODO finish me 


def log10_to_raw_abundance(abundance_df,abundance_scaling_factor=10**6,abundance_floor=0,drop_samples_re="1B",
                            abundance_strain_abbrevs=[]):
    if len(abundance_strain_abbrevs) == 0:
        abundance_strain_abbrevs = STRAIN_ABBREVS
    # STRAIN_ABBREVS = ["Bbr","Bca","Bli2D9","Bli463","Blu","Rob","Dfo","Dlo","Eav","Eco","FprB","Lga4B6","Lru","Mmu","Pco","Pst",
    #               "Rgn","Rto","Sga","Spa"]
    abundance_converted  = abundance_df.copy()
    if abundance_strain_abbrevs[0] not in abundance_converted.columns:
        bact_columns = [col.strip() for col in abundance_converted.columns[3:-1].values]
        abbrevs_map = dict(zip(bact_columns,abundance_strain_abbrevs))
        abundance_converted.rename(columns=abbrevs_map,inplace=True)
    abundance_converted = abundance_converted.loc[:,abundance_strain_abbrevs]
    
    abundance_converted = abundance_converted.replace("ND",0)
    abundance_converted.loc[:,abundance_strain_abbrevs] = abundance_converted.loc[:,abundance_strain_abbrevs].astype(float)
    
    for col in abundance_strain_abbrevs:
        abundance_converted.loc[:,col] = np.power([10]*len(abundance_converted),
                                                        abundance_converted.loc[:,col].values)
        if abundance_floor: #For absolute abundance < abundance floor, replace with abundance_floor 
            abundance_converted.loc[(abundance_converted.loc[:,col] < abundance_floor) & 
                                          (abundance_converted.loc[:,col] != 1),col] = abundance_floor
        abundance_converted.loc[:,col] = abundance_converted.loc[:,col].replace(1,0)
        #Replace 1 (from 10^0 abundance) with 0 
        abundance_converted.loc[:,col] = abundance_converted.loc[:,col]/(abundance_scaling_factor) #Scale by millions of genome equivalents
        
    #Drop 1B samples from abundance 
    if drop_samples_re and abundance_converted.columns.str.contains(drop_samples_re).any():
        sample_arm_re=r'[\w_]+[\w_]+_([12][ABC])_[\w+]'
        abundance_converted.loc[:,"Arm"] = abundance_converted.index.str.extract(sample_arm_re,expand=False)
        abundance_converted = abundance_converted.loc[~(abundance_converted["Arm"]=="1B"),:]
    return abundance_converted

############################ Annotation Loading (mcSEED, CAZyme) ############################
def load_mcSEED(mcSEED_fpath,list_cols=False):
    """
    Read mcSEED annotations from the file at mcSEED_fpath 
    @param: mcSEED_fpath, str: file path to mcSEED table
    @param: list_cols, boolean: if True, semicolon separated strings are converted to lists for 
    the columns ["Functional category","Functional pathway", "Phenotype"]

    @return mcSEED_df, DataFrame: DataFrame with columns from original file, with appropriate 
    columns as list instead of ';'-separated strings if list_cols=True. Add column 
    ["Strain"] for strain abbreviation
    @return pht_pathway_df, DataFrame: DataFrame mapping pathway abbreviations (i.e. Mal) to
    Functional pathway and Functional category (i.e. Maltose Utilization, Carbohydrate utilization)
    """
    #Modified from miRNA_SVD.ipynb
    mcSEED_df = pd.read_csv(mcSEED_fpath)
    mcSEED_df.set_index("Locus tag",inplace=True)
    mcSEED_df.loc[:,"Strain"] = mcSEED_df.loc[:,"Isolate name"].map(dict(zip(mcSEED_df["Isolate name"].unique(),STRAIN_ABBREVS)))
    global MCSEED_ISOLATE_NAMES, STRAIN_ABBREVS_TO_ISOLATE_NAMES
    MCSEED_ISOLATE_NAMES = mcSEED_df["Isolate name"].unique().tolist()
    STRAIN_ABBREVS_TO_ISOLATE_NAMES = dict(zip(STRAIN_ABBREVS,MCSEED_ISOLATE_NAMES))
    #Attempt to process isolate names to unique taxa names (Genus species)
    TAXA_NAMES_MATCH_PAT = r'([A-Z][a-z]+[ _][a-z]+)[ _][\w_]+' #Flexible space or underscore separation, extract genus and species
    #Note that TAXA_NAMES_MATCH_PAT will result in redundant taxa names if two strains of same species are present
    taxa_names =[re.search(TAXA_NAMES_MATCH_PAT,iso_name).group(1) for iso_name in MCSEED_ISOLATE_NAMES]
    taxa_names = [re.sub(r'[a-z]+','.',tn,count=1) for tn in taxa_names]
    taxa_names = [re.sub(r'_',' ',tn,count=1) for tn in taxa_names]
    if len(set(taxa_names)) < len(set(MCSEED_ISOLATE_NAMES)):
        warnings.warn("Isolate names result in redundant taxa names (i.e. two strains of same species are present).")
        warnings.warn("Please use config file to specify the taxa names with unique identifiers.")
    else: 
        global MCSEED_TAXA_NAMES
        MCSEED_TAXA_NAMES = taxa_names
    pht_pathway_df = pd.DataFrame(columns=["Functional pathway","Functional category"])
    for pht_str in mcSEED_df["Phenotype"].unique():
        pht_str_matches = mcSEED_df.loc[mcSEED_df["Phenotype"]==pht_str]
        assert(len(pht_str_matches["Functional pathway"].unique())==1)
        functional_pathway = pht_str_matches.iloc[0]["Functional pathway"]
        functional_category = pht_str_matches.iloc[0]["Functional category"]
        split_phts = [pht.strip() for pht in pht_str.split(";")]
        split_paths = [cat.strip().title() for cat in functional_pathway.split(";")]
        for pht, cat in zip(split_phts,split_paths):
            pht_pathway_df.loc[pht,"Functional pathway"] = cat
        if len(functional_category.split(";")) == 1:
            pht_pathway_df.loc[split_phts,"Functional category"] = functional_category

    if list_cols:
        for col in ["Functional category","Functional pathway", "Phenotype"]:
            mcSEED_df.loc[:,col] = mcSEED_df.loc[:,col].str.split(";")

    
    return mcSEED_df, pht_pathway_df

def load_BPM_df(BPM_fpath):
    bpm_df = pd.read_csv(BPM_fpath)
    bpm_df.set_index("Isolate name",inplace=True)
    bpm_df = bpm_df.transpose().drop(index="# functions")
    return bpm_df

def load_CAZyme(CAZyme_fpath):
    """Load CAZyme annotations from CAZyme_fpath into a DataFrame.
    """
    cazyme = pd.read_csv(CAZyme_fpath,sep='\t',names=["Isolate","Locus tag","GH","Model notes","Defline"])
    cazyme.loc[:,"Strain"] = cazyme.loc[:,"Isolate"].map(dict(zip(cazyme["Isolate"].unique(),CAZYME_STRAIN_ABBREVS)))
    return cazyme 


def vc_strain_locus_tag(expr_df,mcseed_df,drop_loci_tags=[]):
    """Returns a DataFrame indexed on locus tags for each strain, containing columns: 
        "Filtered Loci" - number of loci corresponding to that locus tag (ie loci per strain)
        "Strain" - full strain name corresponding to locus tag 
        :param pd.DataFrame expr_df: DataFrame indexed by ORFs containing (transformed) expression data, columns
        are samples 
        :param pd.DataFrame mcseed_df: DataFrame containing mcSEED annotations for loci, not necessarily for 
        all loci in expr_df 
    """
    vc_by_strain_locus_tag = expr_df.index.str.extract(r'(\w+)_\d+',expand=False).value_counts()
    if len(drop_loci_tags) > 0: 
        filt_drop_loci_tags = [lt for lt in drop_loci_tags if lt in vc_by_strain_locus_tag.index]
        vc_by_strain_locus_tag.drop(filt_drop_loci_tags,inplace=True)
#     vc_by_strain_locus_tag.drop("ROSSTS7063_a2",inplace=True) #2nd3rd_trial specific error 
    locus_tag_isos = [mcseed_df.loc[mcseed_df.index.str.contains(lt),"Isolate name"].values[0] 
                             for lt in vc_by_strain_locus_tag.index] 
    locus_tag_strains = [dict(zip(MCSEED_ISOLATE_NAMES,STRAIN_ABBREVS))[iso] for iso in locus_tag_isos]
    locus_vc_df = pd.DataFrame(index=vc_by_strain_locus_tag.index,
        columns=["Filtered Loci","Isolate name","Strain"])
    locus_vc_df.loc[:,"Filtered Loci"] = vc_by_strain_locus_tag
    locus_vc_df.loc[:,"Isolate name"] = locus_tag_isos
    locus_vc_df.loc[:,"Strain"] = locus_tag_strains

    return locus_vc_df

# def main():
    

# if __name__ == "__main__":
#     main()

