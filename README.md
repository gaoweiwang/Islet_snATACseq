# Download raw data
Download 10x snATAC-seq raw data from NCBI Gene Expression Omnibus (GSE169453)
Download 10x multiome raw data from NCBI Gene Expression Omnibus (GSE200044)

# Download processed data to reproduce figures
## 1, download files less than 25M from github data folder:
https://github.com/gaoweiwang/Islet_snATACseq

## 2, download files from figshare:
multiome: https://figshare.com/articles/dataset/processed_multiome_zip/19497665
snATACseq: https://figshare.com/articles/dataset/processed_snATACseq_zip/19497656

## 3, download large files:
snATAC-seq data (http://169.228.232.194/~mmallick/o/processed_snATACseq.tar.gz)
multiome data (http://169.228.232.194/~mmallick/o/processed_multiome.tar.gz)
  
# Figures:
Use the cell_clustering.ipynb notebook to reproduce cell clustering results
Use the downstream_analysis.ipynb notebook to reproduce downstream analyses
  
# Scripts:
Use AI_Hypothesis_testing.ipynb to train classifier to distinguish beta cells from ND, pre-T2D and T2D donors
Use AI_subtype.ipynb to train classifier to distinguish beta cell subtypes enriched in ND and T2D donors
Use AI_HPAP.ipynb to predict subtype of beta cells from public HPAP isle snATAC-seq dataset
Use run_scWASP.py and run_imbalance.py to perform genetic analysis

# Session information
Platform: x86_64-conda_cos6-linux-gnu (64-bit) \
Running under: CentOS Linux 7 (Core) \
attached base packages: \
 grid  stats4  parallel  stats  graphics  grDevices utils  datasets  methods   base     
