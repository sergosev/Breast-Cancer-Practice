# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os


#set the working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))


print("="*120)
print("🟢 Preparing gene count matrix...")

#make the list for DataFrames
def DF_list_from_dir(path_to_dir, sep, head, skip):
    import pandas as pd
    import os
    
    #make the list of data frames
    data_lst = []
    
    for file in os.listdir(path_to_dir): #a cycle that takes each file name from a given directory
            #gett a full path to a file
            full_path = os.path.join(path_to_dir, file)
            try:
                data = pd.read_csv(full_path, 
                                   sep = sep, 
                                   header = head, 
                                   skiprows = skip,
                                   encoding = "utf-8")
                if data.shape[0] == 0 or data.shape[1] < 2:
                    print(f"Skipping {file}: Empty or malformed")
                    continue
                #add it to the list
                data_lst.append(data)
                
            #there may be a .DS_Store file, ignore it    
            except pd.errors.EmptyDataError:
                print(f"Skipping empty or invalid file: {file}")
            except Exception as e:
                print(f"Error reading {file}: {e}")
    return data_lst

#get the normal data list
norm_data_lst = DF_list_from_dir(path_to_dir = "../data/normal_data/", sep = "\t", head = 0, skip = 1)

#get the tumor data list
tum_data_lst = DF_list_from_dir(path_to_dir = "../data/tumor_data/", sep = "\t", head = 0, skip = 1)


#now we merge all the data frames into two (only leaving columns gene_id, gene_name and unstranded)
normal_data = norm_data_lst[0][["gene_id", "unstranded"]]
tumor_data = tum_data_lst[0][["gene_id", "unstranded"]]

for i in range(1,len(norm_data_lst)):
    normal_data = pd.merge(normal_data, norm_data_lst[i][["gene_id", "unstranded"]], 
                           on = "gene_id", 
                           how = "inner",
                           suffixes=('', f'_{i}'))

for i in range(1,len(tum_data_lst)):
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

normal_data.to_csv("../results/unstranded_normal_counts.csv")
tumor_data.to_csv("../results/unstranded_tumor_counts.csv")
data.to_csv("../results/full_count_matrix.csv")
filtered.to_csv("../results/filtered_count_matrix.csv")

#=====================================================================================================================
#Doing the DE analysis
from pydeseq2.dds import DeseqDataSet

#create the meta data (a df that says which sample is which)
n_lst = ["norm"]*30
t_lst = ["tumor"]*30
meta = pd.DataFrame(n_lst + t_lst, columns = ["condition"], index = filtered.columns)

print("="*120)
print("🟢 Making the DeseqDataSet...")


#make a DeseqDataSet (aka get log2FC values with GLM)
ds = DeseqDataSet(counts = filtered.T, #give it a df where rows = samples, cols = genes
                  metadata = meta, #give it the metadata
                  design = "~condition", #compare based on the condition column from metadata
                  refit_cooks = True) #filter out outliers

print("="*120)
print("🟢 Fitting glm parameters...")

ds.deseq2() #fit LFCs and other statistical stuff (i'm interested in LFCs)

#now we get the results
from pydeseq2.ds import DeseqStats

print("="*120)
print("🟢 Doing the DGEA...")

#perform statistical analysis (make sure that calculated LFCs are not zero!)
res = DeseqStats(ds,
                 contrast = ["condition", "tumor", "norm"],
                 alpha = 0.05)
res.summary() #add the summary to the res DataSet

res.results_df.to_csv("../results/DE_results.csv") #The resulting table is here

print("="*120)
print("🟢 LFC shrinkage...")

#LFC shrinkage
ds.obsm["design_matrix"] #we take the condition[T.tumor] as coeff
res.lfc_shrink(coeff = "condition[T.tumor]")
res.results_df.to_csv("../results/DE_res_LFC_shrank.csv")
#=====================================================================================================================
#Visualisation. VolcanoPlot

print("="*120)
print("🟢 Drawing VolcanoPlot...")


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

print("="*120)
print("🟢 Performing GO Enrichment Analysis...")

#filtering out upregulated genes
upregulated = DE_res[(DE_res["padj"] < 0.05) & (DE_res["log2FoldChange"] >= 2.5)]

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


#Doing the GO Analysis via Enrichr from GSEApy
from gseapy import enrichr
GO_results = enrichr(
    gene_list = upregulated['symbol'].dropna().tolist(),
    gene_sets = 'GO_Cellular_Component_2021',
    organism = "Human",
    cutoff = 0.05)
    
#Now I save the results
GO_df = GO_results.res2d.copy()
GO_df.to_csv("../results/GO_enrichment_results.csv")


#=====================================================================================================================
#Downloaded surfaceome table of proteins from https://wlab.ethz.ch/surfaceome/
surfaceome = pd.read_csv("../data/surfaceome_Wollscheid_lab.txt", 
                   sep = "\t", 
                   header = 0, 
                   skiprows = 1,
                   encoding = "latin1")

print("="*120)
print("🟢 Intersecting DGEA results with surfaceome data...")



#I have a DF with upregulated genes and their symbols and IDs
#Intersecting those with surfaceome
surfaceome["UniProt name"] = surfaceome["UniProt name"].str.slice(0, -6) 
surface_targets = upregulated[upregulated["symbol"].isin(surfaceome["UniProt name"])]

#Sorting by LFCs
surface_targets = surface_targets.sort_values("log2FoldChange", ascending = False)
surface_targets.to_csv("../results/upregulated_surface_genes.csv")

print("\n✅ All Done!")



