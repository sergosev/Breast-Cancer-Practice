# Breast Cancer RNA-seq Differential Expression & CAR-T Target Discovery

This project performs a full RNA-seq differential expression analysis (DEA) of breast cancer vs. normal tissues using STAR-count data from the TCGA-BRCA project. It identifies significantly upregulated genes and explores their potential as CAR-T cell therapy targets using gene ontology and surfaceome annotation.
The STAR-counts tables were downloaded from GDCA, TCGA-project. The tumor data was collected from primary solid tissue samples of Ductal and Lobular neoplasms. The normal data was collected from healthy tissues of the patients with the same cancer type.
I am using PyDESeq package for Differential Gene Expression Analysis (DGEA), Enrichr from GSEApy for Gene Onthology Enrichment Analysis. 
---

# 📁 Project Structure
|\n
|---data/					<-GDC manifests for downloading raw data and raw data gene counts (not included in the repository)\n
|---scripts/\n
|      |---01_GDC_download.bash			<-bash script for downloading raw STAR-counts from GDC (TCGA project)\n
|      |---02_data_unpacking.bash		<-bash script for preparing all the downloaded tables for differential expression analysis\n
|      |---03_BRC_DGEA.py			<-Python script for differential expression analysis and CAR-T cell therapy target search\n
|---results/					<-.csv tables and plots derived during the pipeline work\n
|---gdc-client.exe				<-gdc client used for downloading raw data for the GDC portal\n
|---README.md\n
|---LICENSE.txt\n
|---.gitignore\n

---

# Workflow
## 1: Downloading raw gene counts tables - /scripts/01_GDC_dosnload.bash
The bash script "01_GDC_download.bash" uses manifests from /data/manifests to download 30 .tsv STAR-counts tables for tumor tissue and 29 tables for normal tisue from GDC data portal. 
Tumor type: Breast Cancer, Ductal and Lobular neoplasm.
Tumor sample: primary solid tumor.
Normal sample: healthy solid tissue	

## 2: Organising the data for the analysis - /scripts/02_data_unpacking.bash
This script deletes anything non-related to the data analysis and puts all the .tsv tables in the correct folders.

## 3: Differential Gene Expression Analysis - /scripts/03_BRC_DGEA.py
Steps of this Python script:\n
	- Upload all the gene counts tables and create a gene counts matrix\n
	- Create a DESeq DataSet\n
	- Get the data frame with log2FoldChange values and their p-values\n
	- Apply LFC shrinkage\n
	- Filter out upregulated genes (p-value < 0.05 and LFC > 2.5)\n
	- Conduct a Gene Onthology Enrichment Analysis via Enrichr from GSEApy\n
	- Intersect the DGEA results data frame with the surfaceome data frame from https://wlab.ethz.ch/surfaceome/\n

---

# Requirements:
Python 3.12.17\n
Install via 'pip install -r requirements.txt'\n
Minimal dependencies:\n
'''txt
pandas
pydeseq2
matplotlib
seaborn
numpy
mygene
gseapy
bipython'''
\n
For Enrichr and mygene you will need Internet access.

---

## What's done
- Full DGEA pipeline
- Ensembl-to-symbol mapping
- GO Enrichment Analysis
- Surfaceome cross-referencing
- Potential target list generation

---

## Plans for future
- Results visualisation
- Search for genes that are not expressed in normal tissue
- Analysis of the results and literature review

---

## Git Ignore
*.tsv\n
*.log\n
.DS_Store\n
/data/*_data/*\n

all important results are uploaded as .csv files in the /results/ folder.
Raw gene_counts tables are ignored as well as .log files.

---
## Contact
Feel free to contact me via e-mail or Telegram. I would be glad to hear your feedback and ideas.\n
sergey.losev.01@gmail.com\n
TG: @small_party
