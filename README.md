This is my bioinformatics self-study project.
I am trying to conduct a differential gene expression analysis of Breast Cancer in order to determine whether there is a potential target for CAR-T cell therapy.

For that I am downloading 30 STAR-counts for normal tissues and 30 STAR-counts tables for tumor tissue, specifically primary tumor, solid tissues only.
The tables come from the GDC database, the TCGA-BRCA project, all cases associated with Ductal and Lobular Neoplasms.
The STAR-Count tables contain the amount of Illumina-seq reads that got mapped for each annotation in human genome.
Manifests with references to all the tables are in the "manifests" directory.
There is also a script ('GDC_download.bash) for downloading all the data in the right directories.
Bash script "file_unpacking.bash" is for preparing the data for my Python script. 

I am using Python for data analysis and target search. Specifically these libraries:
- pandas
- numpy
- PyDESeq
- matplotlib
- and others (to be  added)
These libraries need to be install in order for the pipeline to work.

The idea is to find highly expressed genes in cancerous samples, that are associated with plasma membrane or cell surface.

At this point I managed to do the DGE analysis of some random data I choose blindly at GDC.
Unfortunately that data was not good for my Gene Onthology Enrichment Analysis, thus I can't find which genes are associated with various cell components.
However I have found a table ("Surfaceome.txt") that contains all proteins, localised on the cell surface.
Untersecting my DGE results with it gave me a list of upregulated genes, probably associated with cell surface.

Now I wanna do it all on a +- OK data_set and prepare a pipeline to present my project to some friends and maybe my personal blog.
