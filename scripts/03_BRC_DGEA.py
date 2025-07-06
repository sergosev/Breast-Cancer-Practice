# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os
import numpy as np


#set the working directory
os.chdir("/Users/smallparty/Desktop/Breast Carcinoma Practice/scripts")

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

print("="*120)
print("🟢 Making the DeseqDataSet...")


#create the meta data (a df that says which sample is which)
n_lst = ["norm"]*len(norm_data_lst)
t_lst = ["tumor"]*len(tum_data_lst)
meta = pd.DataFrame(n_lst + t_lst, columns = ["condition"], index = filtered.columns)

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
#Gene Onthology Analysis (a failure)
#I need genes that are associated with plasma membrane or cell surface
DE_res = res.results_df.copy() #get a full copy of the result data frame
#mark the over/underexpressed genes with pandas subsetting and addint a new Bool column to the results df
DE_res["significant"] = (DE_res["padj"] < 0.001) & (abs(DE_res["log2FoldChange"]) >= 2)
DE_res["padj"] = DE_res["padj"].replace(0, np.nan) #all zero p-adjusted become NAs
DE_res = DE_res.dropna(subset = ["padj"]) #get rid of NAs

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
GO_CC_res = enrichr(
    gene_list = upregulated['symbol'].dropna().tolist(),
    gene_sets = 'GO_Cellular_Component_2021',
    organism = "Human",
    cutoff = 0.05)

GO_BP_res = enrichr(
    gene_list = upregulated['symbol'].dropna().tolist(),
    gene_sets = 'GO_Biological_Process_2021',
    organism = "Human",
    cutoff = 0.05)

GO_MF_res = enrichr(
    gene_list = upregulated['symbol'].dropna().tolist(),
    gene_sets = 'GO_Molecular_Function_2021',
    organism = "Human",
    cutoff = 0.05)
    
#Now I save the results
GO_CC_df = GO_CC_res.res2d.copy()
GO_BP_df = GO_BP_res.res2d.copy()
GO_MF_df = GO_MF_res.res2d.copy()
GO_CC_df.to_csv("../results/GO_CC_results.csv")
GO_BP_df.to_csv("../results/GO_BP_results.csv")
GO_MF_df.to_csv("../results/GO_MF_results.csv")

#=====================================================================================================================
#I have a list of cell surface proteins IDs from https://wlab.ethz.ch/surfaceome/
surface_ids = pd.read_csv("../data/surfaceome_ids.txt",
                          sep = "\t",
                          names = ["gene_id"],
                          encoding = "latin1")
surface_ids["gene_id"] = surface_ids["gene_id"].str.slice(0, -6)

print("="*120)
print("🟢 Intersecting DGEA results with surfaceome data...")

#Intersecting upregulated genes with surfaceome data
candidates = upregulated[upregulated["symbol"].isin(surface_ids["gene_id"])]
candidates = candidates.sort_values("log2FoldChange", ascending = False) #sorting by LFCs
candidates.to_csv("../results/upregulated_surface_genes.csv") #exporting
#=====================================================================================================================
#Found upregulated surface genes. Now I am looking for the ones with no or low expression in normal tissue
#For that I need to normalise my raw counts, using CPM normalisation method

print("="*120)
print("🟢 Filtering out genes that are highly expressed in normal tissue...")

normal_data = normal_data.set_index("gene_id")
tumor_data = tumor_data.set_index("gene_id")

#getting library sizes
norm_lib_sizes = normal_data.sum(axis=0)
tum_lib_sizes = tumor_data.sum(axis=0)

#dividing each gene counts for each sample by sample's lib size
norm_norm_df = normal_data.div(norm_lib_sizes, axis = 1) * 1e6
norm_tum_df = tumor_data.div(tum_lib_sizes, axis = 1) * 1e6

#setting IDs in counts tables from ENSG*******.** to ENSG*******
norm_norm_df = norm_norm_df.reset_index()
norm_norm_df["gene_id"] = norm_norm_df["gene_id"].str.replace(r'\.\d+', '', regex = True)

#merging upregulated sirface genes with normalised normal gene counts table
merged = candidates.merge(norm_norm_df, left_on = "gene_id", right_on = "gene_id")

#calculatin mean gene expression in normal tissue
norm_norm_df["normal_mean_CPM"] = norm_norm_df.drop(columns=["gene_id"]).mean(axis=1)

#merging the means to the surface targets table
merged = merged.merge(norm_norm_df[["gene_id", "normal_mean_CPM"]], on="gene_id", how="left")

#filtering out final candidates
final_candidates = merged[merged["normal_mean_CPM"] < 5][["gene_id", "symbol", "log2FoldChange", "lfcSE", "padj", "normal_mean_CPM"]]
final_candidates.to_csv("../results/final_candidates.csv")


#=====================================================================================================================
#Annotating the final candidates and filtering out the ones expressed in brain or any other crucial tissue
print("="*120)
print("🟢 Annotating found candidates...")

fin_cand_symb = final_candidates["symbol"].tolist() #getting list of gene names
mg = myg.MyGeneInfo() #getting a DF from BioMart

annotations = mg.querymany(fin_cand_symb,
                           scopes = "symbol",
                           fields = ["name", "summary", "refseq"],
                           species = "human")
annotations_df = pd.json_normalize(annotations)
fin_candidates_annotation = final_candidates.merge(
    annotations_df[["query", "name", "summary"]], 
    left_on = "symbol", right_on = "query", how = "left")
fin_candidates_annotation.drop("query", axis = 1, inplace = True)

fin_candidates_annotation.to_csv("../results/annotated_fin_candidates.csv")

#=====================================================================================================================
#Visualising
print("="*120)
print("🟢 Drawing pictures...")

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import mygene as myg

#getting DE results again
DE_res = res.results_df.copy()
DE_res.reset_index(inplace = True)
DE_res["gene_id"] = DE_res["gene_id"].str.replace(r'\.\d+', '', regex = True)

#adding gene symbols
DE_res_ids = DE_res["gene_id"].tolist()
mg = myg.MyGeneInfo()
gene_symbols = mg.querymany(DE_res_ids, scopes = "ensembl.gene", fields = "symbol", species = "human")
gene_df = pd.DataFrame(gene_symbols)
gene_df = gene_df[['query', 'symbol']].dropna()
DE_res = DE_res.merge(gene_df, left_on = "gene_id", right_on = "query") #do the merge

#Marking significant LFCs and potential targets
DE_res["significant"] = (DE_res["padj"] < 0.05) & (abs(DE_res["log2FoldChange"]) >= 1)

#Filtering out 0 LFCs
DE_res_filtered = DE_res[abs(DE_res["log2FoldChange"]) > 0.01]

#Drawing VolcanoPlot_______________________________________________________________________________
plt.figure(figsize=(10, 6)) #set the plot size

DE_res["Adjusted p-value"] = DE_res["significant"].map({True: "< 0.05", False: "> 0.05"})

sns.scatterplot(
    data=DE_res, x="log2FoldChange", y=-np.log10(DE_res["padj"]),
    hue="Adjusted p-value", palette={"< 0.05": "red", "> 0.05": "grey"},
    alpha=0.5, edgecolor=None
)

plt.axhline(-np.log10(0.05), linestyle='--', color='black') #draw p-value threshold
plt.axvline(-1, linestyle='--', color='black') #draw negative LFC threshhold
plt.axvline(1, linestyle='--', color='black') #draw positive LFC threshold
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 Adjusted p-value")
plt.title("Volcano Plot")
plt.tight_layout()
plt.savefig("../results/volcano_targets.png")
plt.show()

#Drawing MA Plot___________________________________________________________________________________
plt.figure(figsize=(10,6))
sns.scatterplot(
    x = "baseMean",
    y = "log2FoldChange",
    data = DE_res,
    hue = DE_res["padj"] < 0.05,
    alpha=0.5, #transparency
)
plt.legend(title = "Adjusted p-value", labels = ["< 0.05", "> 0.05"])
plt.axhline(0, color="black", linestyle="--") #horizontal line at 0 LFC
plt.xscale("log") #make x scale logarithmic
plt.xlabel("Mean expression (baseMean)")
plt.ylabel("Log2 Fold Change")
plt.title("MA Plot")
plt.savefig("../results/ma_plot.png")
plt.show()

#Drawing LFC histogram______________________________________________________________________________
plt.figure(figsize=(8,5))
sns.histplot(DE_res["log2FoldChange"], #data
             bins=100,  #number of bars
             color="steelblue", #color of bars
             kde=True) #draw a line based on distribution

plt.axvline(0, color="black", linestyle="--") #draw a vertical at zero

plt.title("Distribution of Log2 Fold Changes")
plt.xlabel("log2 Fold Change")
plt.ylabel("Number of Genes")
plt.tight_layout()
plt.savefig("../results/log2fc_histogram.png")
plt.show()

plt.figure(figsize=(8,5))
sns.histplot(DE_res_filtered["log2FoldChange"], 
             bins=100, 
             color="steelblue", 
             kde=True)

plt.axvline(0, color="black", linestyle="--")

plt.title("Distribution of Log2 Fold Changes")
plt.xlabel("log2 Fold Change withour close to zero values")
plt.ylabel("Number of Genes")
plt.tight_layout()
plt.savefig("../results/log2fc_histogram_wo_zeros.png")
plt.show()
print("\n✅ All Done!")



