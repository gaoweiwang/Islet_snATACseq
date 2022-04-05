# Download raw data
Download 10x snATAC-seq raw data from NCBI Gene Expression Omnibus (GSE169453)
Download 10x multiome raw data from NCBI Gene Expression Omnibus (GSE200044)

# Download processed snATAC-seq and multiome data to reproduce figures
## 1, download files less than 25M from github data folder: https://github.com/gaoweiwang/Islet_snATACseq

or, from figshare:
multiome: https://figshare.com/articles/dataset/processed_multiome_zip/19497665
snATACseq: https://figshare.com/articles/dataset/processed_snATACseq_zip/19497656

## 2, download large processed or intermeidate files:
snATAC-seq data (http://169.228.232.194/~mmallick/o/processed_snATACseq.tar.gz)
multiome data (http://169.228.232.194/~mmallick/o/processed_multiome.tar.gz)
  
## 1. Use the cell_clustering.ipynb notebook to cluster cells.
### Session information
anndata     0.7.4 \
scanpy      1.6.0 \
sinfo       0.3.1 \
PIL                 7.2.0 \
anndata             0.7.4 \
backcall            0.2.0 \
bottleneck          1.3.2 \
cairo               1.19.1 \
cffi                1.14.0 \
cloudpickle         1.6.0 \
colorama            0.4.3 \
cycler              0.10.0 \
cytoolz             0.11.0 \
dask                2.30.0 \
dateutil            2.8.1 \
decorator           4.4.2 \
h5py                2.10.0 \
igraph              0.8.3 \
importlib_metadata  1.7.0 \
ipykernel           5.3.4 \
ipython_genutils    0.2.0 \
ipywidgets          7.5.1 \
jedi                0.17.2 \
jinja2              2.11.2 \
joblib              0.16.0 \
kiwisolver          1.2.0 \
legacy_api_wrap     0.0.0 \
leidenalg           0.8.3 \
louvain             0.6.1 \
markupsafe          1.1.1 \
matplotlib          2.2.3 \
natsort             7.0.1 \
numba               0.50.1 \
numexpr             2.7.1 \
numpy               1.17.0 \
packaging           20.4 \
pandas              1.1.1 \
parso               0.7.0 \
patsy               0.5.1 \
pexpect             4.8.0 \
pickleshare         0.7.5 \
prompt_toolkit      3.0.5 \
psutil              5.6.3 \
ptyprocess          0.6.0 \
pycparser           2.20 \
pygments            2.6.1 \
pyparsing           2.4.7 \
pytz                2020.1 \
rpy2                3.1.0 \
scanpy              1.6.0 \
scipy               1.5.2 \
seaborn             0.10.1 \
sinfo               0.3.1 \
six                 1.15.0 \
sklearn             0.22 \
statsmodels         0.11.1 \
tables              3.6.1 \
texttable           1.6.3 \
tlz                 0.11.0 \
toolz               0.11.1 \
tornado             6.0.4 \
traitlets           4.3.3 \
umap                0.4.6 \
wcwidth             0.2.5 \
yaml                5.3.1 \
zmq                 19.0.1 \
IPython             7.17.0 \
jupyter_client      6.1.6 \
jupyter_core        4.6.3 \
jupyterlab          2.2.6 \
notebook            6.0.3 \

## 2. Use the downstream_analysis.ipynb notebook to reproduce figures in the manuscript.
### Session information
Platform: x86_64-conda_cos6-linux-gnu (64-bit) \
Running under: CentOS Linux 7 (Core) \
attached base packages: \
 grid  stats4  parallel  stats  graphics  grDevices utils  datasets  methods   base     

other attached packages: \
 ggpubr_0.2.5 \
 magrittr_1.5           
 EnhancedVolcano_1.2.0  \
 ggrepel_0.8.2      
 pheatmap_1.0.12 \
 VennDiagram_1.6.20         
 futile.logger_1.4.3 \
 DESeq2_1.22.1              
 SummarizedExperiment_1.12.0  \
 DelayedArray_0.8.0         
 matrixStats_0.56.0   \
 GenomicRanges_1.34.0       
 GenomeInfoDb_1.18.2   \
 IRanges_2.16.0             
 S4Vectors_0.20.1  \
 ggfortify_0.4.10           
 ggplot2_3.3.2   \
 Biobase_2.42.0             
 BiocGenerics_0.28.0 \
 gPCA_1.0                   
 sva_3.30.1   \
 BiocParallel_1.16.6        
 genefilter_1.64.0  \
 mgcv_1.8-33                
 nlme_3.1-149         
