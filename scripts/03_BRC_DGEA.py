# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#got 2 attempts with this
#the original data is with no number
#second attemt data folders and manifests are marked by number 2
#for second attempt all data was manually selected to be for Ductal and Lobular Neoplasms
#This script is for attempt 1, gotta change filenames for attempt 2

import pandas as pd
import os

#set the working directory
os.chdir("/Users/smallparty/Desktop/Breast Carcinoma Practice")

#make the list for DataFrames
def DF_list_from_dir(path_to_dir, sep, head, skip):
    import pandas as pd
    import os
    
    names = os.listdir(path_to_dir) #get the list of file names
    
    #make the list of data frames
    data_lst = []
    for i in names:
        try:
            #read each DF from a directory
            data = pd.read_csv(path_to_dir + i, 
                               sep = sep, 
                               header = head, 
                               skiprows = skip)
            #add it to the list
            data_lst.append(data)
        except UnicodeDecodeError:
            continue
    return data_lst

#get the normal data list
norm_data_lst = DF_list_from_dir(path_to_dir = "normal_data/", sep = "\t", head = 0, skip = 1)

#get the tumor data list
tum_data_lst = DF_list_from_dir(path_to_dir = "tumor_data/", sep = "\t", head = 0, skip = 1)


#now we merge all the data frames into two (only leaving columns gene_id, gene_name and unstranded)
normal_data = norm_data_lst[0][["gene_id", "unstranded"]]
tumor_data = tum_data_lst[0][["gene_id", "unstranded"]]

for i in range(1,len(norm_data_lst)):
    normal_data = pd.merge(normal_data, norm_data_lst[i][["gene_id", "unstranded"]], 
                           on = "gene_id", 
                           how = "inner",
                           suffixes=('', f'_{i}'))
    tumor_data = pd.merge(tumor_data, tum_data_lst[i][["gene_id", "unstranded"]], 
                           on = "gene_id", 
                           how = "inner",
                           suffixes=('', f'_{i}'))


normal_data.columns = ['gene_id'] + [f'wt_{i+1}' for i in range(len(norm_data_lst))]
tumor_data.columns = ['gene_id'] + [f'tm_{i+1}' for i in range(len(tum_data_lst))]

#now we merge the two count data frames into one and it will be the data frame to work with
data = pd.merge(normal_data, tumor_data, on = "gene_id", how = "inner")
data = data.set_index("gene_id")

#filter out genes with low counts
filtered = data[data.sum(axis=1) > 10]

#=====================================================================================================================
#Doing the DE analysis
from pydeseq2.dds import DeseqDataSet

#create the meta data (a df that says which sample is which)
n_lst = ["norm" for _ in range(30)]
t_lst = ["tumor" for _ in range(30)]
meta = pd.DataFrame(n_lst + t_lst, columns = ["condition"], index = filtered.columns)

#make a DeseqDataSet (aka get log2FC values with GLM)
ds = DeseqDataSet(counts = filtered.T, #give it a df where rows = samples, cols = genes
                  metadata = meta, #give it the metadata
                  design = "~condition", #compare based on the condition column from metadata
                  refit_cooks = True) #filter out outliers

ds.deseq2() #fit LFCs and other statistical stuff (i'm interested in LFCs)

#now we get the results
from pydeseq2.ds import DeseqStats

#perform statistical analysis (make sure that calculated LFCs are not zero!)
res = DeseqStats(ds,
                 contrast = ["condition", "tumor", "norm"],
                 alpha = 0.05)
res.summary() #add the summary to the res DataSet
res.results_df #The resulting table is here

#subset the statistically significant LFCs that are 4 and bigger
res.results_df[(res.results_df["log2FoldChange"] >= 4) & (res.results_df["padj"] <= 0.05)]

#LFC shrinkage
ds.obsm["design_matrix"] #we take the condition[T.tumor] as coeff
res.lfc_shrink(coeff = "condition[T.tumor]")

#=====================================================================================================================
#Visualisation. VolcanoPlot
DE_res = res.results_df.copy() #get a full copy of the result data frame

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#mark the over/underexpressed genes with pandas subsetting and addint a new Bool column to the results df
DE_res["significant"] = (DE_res["padj"] < 0.001) & (abs(DE_res["log2FoldChange"]) >= 2)

DE_res["padj"] = DE_res["padj"].replace(0, np.nan) #all zero p-adjusted become NAs
DE_res = DE_res.dropna(subset = ["padj"]) #get rid of NAs



plt.figure(figsize = (15, 9))
sns.scatterplot(
    data = DE_res, #set the data frame
    x = "log2FoldChange", #what column is for X-axis
    y = -np.log10(DE_res["padj"]), #what colump is for Y-axis (negative log10 of p-adjusted)
    hue = "significant", #what to color the dots by
    palette = {True: "red", False: "grey"}, #how to color the dots
    edgecolor = None,
    alpha = 0.5,
    size = 0.5)

plt.xlabel("log2 Fold Change") #add X-axis label
plt.ylabel("-log10(adjusted p-value)") #add Y-axis label
plt.title("Volcano Plot") #add title
plt.axhline(-np.log10(0.001), color='black', linestyle='--', linewidth=1)  # p-value threshold
plt.axvline(-2, color='blue', linestyle='--', linewidth=1)  # LFC threshold
plt.axvline(2, color='blue', linestyle='--', linewidth=1)  # LFC threshold

plt.show() #view the plot


#=====================================================================================================================
#Gene Onthology Analysis (a failure)
#I need genes that are associated with plasma membrane or cell surface
import mygene as myg 
from gprofiler import GProfiler

#filtering out upregulated genes
upregulated = DE_res[(DE_res["padj"] < 0.05) & (DE_res["log2FoldChange"] >= 1)]

#change IDs from ENSG*******.** to ENSG******* and save them to a list
ensembl_ids = upregulated.index.str.replace(r'\.\d+', '', regex = True).tolist()

#get gene names from mygene
mg = myg.MyGeneInfo() #get an object with methods to get gene names
gene_info = mg.querymany(ensembl_ids, scopes = "ensembl.gene", fields = "symbol", species = "human") #get a list of gene names

#convert to DataFrame and merge with upregulated DataFrame
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()
upregulated = upregulated.reset_index() #turn gene_ids into a column
upregulated["gene_id"] = upregulated["gene_id"].str.replace(r'\.\d+', '', regex = True) #change IDs for the merge
upregulated = upregulated.merge(gene_map, left_on = "gene_id", right_on = "query") #do the merge


#Doing the GO Enrichment analysis via Gprofiler (it doesn't work)
gp = GProfiler(return_dataframe = True) #get an object/database with methods to do a GO analysis
GO_results = gp.profile(organism = "hsapiens", 
                        query = upregulated['symbol'].tolist(),
                        user_threshold = 0.05,
                        sources = ["GO:CC"])


#Doing the GO Analysis via Enrichr from GSEApy
from gseapy import enrichr
GO_results = enrichr(
    gene_list = upregulated['symbol'].dropna().tolist(),
    gene_sets = 'GO_Cellular_Component_2021',
    organism = "Human",
    cutoff = 0.05)
    
#Now I filter the results
GO_df = GO_results.res2d.copy()

#maybe something wrong with statistics, I get 4-7k of upregulated genes (that's too much)
#gotta recheck the DE analysis part
#GO Analysis does not give any meaningful results, I will try to do it in a straightforward way
#Also gotta check it for other data

#=====================================================================================================================
#Downloaded surfaceome table of proteins from https://wlab.ethz.ch/surfaceome/
surfaceome = pd.read_csv("surfaceome_Wollscheid_lab.txt", 
                   sep = "\t", 
                   header = 0, 
                   skiprows = 1,
                   encoding = "latin1")

#filtering out upregulated genes
upregulated = DE_res[(DE_res["padj"] < 0.05) & (DE_res["log2FoldChange"] >= 1)]

#change IDs from ENSG*******.** to ENSG******* and save them to a list
ensembl_ids = upregulated.index.str.replace(r'\.\d+', '', regex = True).tolist()

#get gene names from mygene
mg = myg.MyGeneInfo() #get an object with methods to get gene names
gene_info = mg.querymany(ensembl_ids, scopes = "ensembl.gene", fields = "symbol", species = "human") #get a list of gene names

#convert to DataFrame and merge with upregulated DataFrame
gene_map = pd.DataFrame(gene_info)[['query', 'symbol']].dropna()
upregulated = upregulated.reset_index() #turn gene_ids into a column
upregulated["gene_id"] = upregulated["gene_id"].str.replace(r'\.\d+', '', regex = True) #change IDs for the merge
upregulated = upregulated.merge(gene_map, left_on = "gene_id", right_on = "query") #do the merge

#Now I have a DF with upregulated genes and their symbols and IDs
#Intersecting those with surfaceome
surfaceome["UniProt name"] = surfaceome["UniProt name"].str.slice(0, -6) 
surface_targets = upregulated[upregulated["symbol"].isin(surfaceome["UniProt name"])]

#Sorting by LFCs
surface_targets = surface_targets.sort_values("log2FoldChange", ascending = False)






