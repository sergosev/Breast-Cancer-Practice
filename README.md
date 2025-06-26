# Breast Cancer RNA-seq Differential Expression & CAR-T Target Discovery

This project performs a full RNA-seq differential expression analysis (DEA) of breast cancer vs. normal tissues using STAR-count data from the TCGA-BRCA project. It identifies significantly upregulated genes and explores their potential as CAR-T cell therapy targets using gene ontology and surfaceome annotation.

The STAR-counts tables were downloaded from GDCA, TCGA-project. The tumor data was collected from primary solid tissue samples of Ductal and Lobular neoplasms. The normal data was collected from healthy tissues of the patients with the same cancer type.

I am using PyDESeq package for Differential Gene Expression Analysis (DGEA), Enrichr from GSEApy for Gene Onthology Enrichment Analysis. 

---

# 📁 Project Structure
- data/: GDC manifests for downloading raw data and raw data gene counts (not included in the repository)
- scripts/
  - 01_GDC_download.bash: bash script for downloading raw STAR-counts from GDC (TCGA project)
  - 02_data_unpacking.bash: bash script for preparing all the downloaded tables for differential expression analysis
  - 03_BRC_DGEA.py: Python script for differential expression analysis and CAR-T cell therapy target search
  - setup.bash: a script for creating a virtual environment, installing all needed libraries from requirements.txt
- results/: .csv tables and plots derived during the pipeline work
- gdc-client.exe: gdc client used for downloading raw data for the GDC portal
- README.md
- LICENSE.txt
- .gitignore

---

# Workflow
All scripts are run in the Terminal from the main folder via ./scripts/"script_name". They worked on my MacOS with Python 3.10. I tried adapting them for Ubuntu, but still gotta check whether they will work on other computers.
The scripts are location-dependent.

## **Setup**
Run /scripts/setup.bash. It will create and activate a virtual environment and install all required libraries. Run this script only once!

After running the script run "source ./env/bin/activate". When finished working run "deactivate". Activate the venv each time you need to run the Python script.

You can skip this step if you have all needed libraries installed or go to Requirements and install the needed libraries without a virtual environment.

## **1: Downloading raw gene counts tables - /scripts/01_GDC_dosnload.bash**
The bash script "01_GDC_download.bash" uses manifests from /data/manifests to download 30 .tsv STAR-counts tables for tumor tissue and 29 tables for normal tisue from GDC data portal. 

- Tumor type: Breast Cancer, Ductal and Lobular neoplasm.
- Tumor sample: primary solid tumor.
- Normal sample: healthy solid tissue	

Make sure you have the GDC Transfer Tool installed in the main folder! Download the archive with the Tool from https://gdc.cancer.gov/access-data/gdc-data-transfer-tool. Unpack the archive to the main folder.

## **2: Organising the data for the analysis - /scripts/02_data_unpacking.bash**
This script deletes anything non-related to the data analysis and puts all the .tsv tables in the correct folders.

## **3: Differential Gene Expression Analysis - /scripts/03_BRC_DGEA.py**
Steps of this Python script:
- Upload all the gene counts tables and create a gene counts matrix
- Create a DESeq DataSet
- Get the data frame with log2FoldChange values and their p-values
- Apply LFC shrinkage
- Filter out upregulated genes (p-value < 0.05 and LFC > 2.5)
- Conduct a Gene Onthology Enrichment Analysis via Enrichr from GSEApy
- Intersect the DGEA results data frame with the surfaceome data frame from https://wlab.ethz.ch/surfaceome/


---

# Requirements:
Python 3.10

Install needed libraries via 'pip install -r requirements.txt' from the main directory or run "./scripts/setup.bash".

Minimal dependencies:
- pandas
- pydeseq2
- matplotlib
- seaborn
- numpy
- mygene
- gseapy
- biopython

For Enrichr and mygene you will need Internet access.

---

## What's done
- Full DGEA pipeline
- Ensembl-to-symbol mapping
- GO Enrichment Analysis
- Surfaceome cross-referencing
- Potential target list generation
- Tried adapting bash scripts to work on Ubuntu as well as MacOS

---

## Plans for future
- Results visualisation
- Search for genes that are not expressed in normal tissue
- Analysis of the results and literature review
- Further adaptations and optimisation of all scripts

---

## Git Ignore
- *.tsv
- *.log
- .DS_Store
- /data/*_data/*

All important results are uploaded as .csv files in the /results/ folder.

Raw gene_counts tables are ignored as well as .log files.

---
## Contact
Feel free to contact me via e-mail or Telegram. I would be glad to hear your feedback and ideas.

sergey.losev.01@gmail.com

TG: @small_party
