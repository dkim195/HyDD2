# HyDD2
HyDD2: Analysis of Hypothalamic Development Database Version 2

**Project Description**
HyDD2 is an integrative bioinformatics project aimed at studying gene expression dynamics across multiple stages of hypothalamic development. This repository includes analysis scripts for processing, clustering, and visualizing single-cell RNA sequencing (scRNA-seq) data, focusing on the regulatory networks governing hypothalamic differentiation and maturation over key developmental stages.

This project is structured to support researchers interested in single-cell sequencing, developmental biology, and the underlying mechanisms of neurogenesis in the hypothalamus.

**#Repository Structure (/Codes)**
.:
1. scRNA/
2. scATAC/
3. Analysis/
4. SCENICPlus/
5. Xenium/
6. DeepLabCut/
repo_structure.txt

./1. scRNA:
1. Each Age/
2. Merging/
3. Mutant/
4. Useful/

./1. scRNA/1. Each Age:
E11_Preprocessing.R
E12_Preprocessing.R
E13_Preprocessing.R
E14_Preprocessing.R
E15_Preprocessing.R
E16_Preprocessing.R
E18_Preprocessing.R
P4_Preprocessing.R
P8_Preprocessing.R

./1. scRNA/2. Merging:
E11_E14_Merging.R
E15_P8_Merging.R
HyDD2_E12_Atlas.R
HyDD2_E16_Atlas.R
P8_Part1.R
P8_Part2.R
Pattern_Neurogenesis_Merge.R

./1. scRNA/3. Mutant:
Dlx1_2_Preprocessing.R
Isl1_Preprocessing.R
Lhx1_Preprocessing.R
Nkx2-2_Preprocessing.R
P8_Atlas-Dlx.R
ZI_TRN_snRNA_P21_DLX.R

./1. scRNA/4. Useful:
DEG_Calculating.R
Integration/
LoupeR.R
Mutant_SubCluster_Comparison.R

./1. scRNA/4. Useful/Integration:
integration.R
objects.R
utilities.R

./2. scATAC:
1. Each Age/
2. Merge/
3. Mutant/
4. Useful/

./2. scATAC/1. Each Age:
Signac_Processing.R

./2. scATAC/2. Merge:
Signac_Merge_Pattern.R

./2. scATAC/3. Mutant:
Signac_Merging_Dlx.R
Signac_Merging_Foxd1.R
Signac_Merging_Isl1.R

./2. scATAC/4. Useful:
PEAK_Calculating.R
Signac_Peak_Motif_FP.R
Signac_Peak_Mutant.R
Signac_RNA_TransferLabel.R

./3. Analysis:
GRN/
GWAS/
Mutant/
NPC/
Plot/
Th/
Xenium/

./3. Analysis/GRN:
GRN_Update_Region_SENIC.R
Neurogenesis_GRN.R
Pattern_GRN.R

./3. Analysis/GWAS:
Compile.R
GWAS.R

./3. Analysis/Mutant:
E12Mutant/
MutantPlot/
P21Dlx/
P8Dlx/
Plot/

./3. Analysis/Mutant/E12Mutant:
Dlx/
Isl1/
Lhx1/
Nkx2-2/

./3. Analysis/Mutant/E12Mutant/Dlx:
Fig3_scATAC_Mutant_Dlx.R
Fig3_scRNA_Mutant_Dlx.R
Script_Dlx.R

./3. Analysis/Mutant/E12Mutant/Isl1:
Fig3_scATAC_Mutant_Isl1.R
Fig3_scRNA_Mutant_Isl1.R
Script_Isl1.R

./3. Analysis/Mutant/E12Mutant/Lhx1:
Fig3_scRNA_Mutant_Lhx1.R
Script_Lhx1.R

./3. Analysis/Mutant/E12Mutant/Nkx2-2:
Fig3_scRNA_Mutant_Nkx22.R
Script_Nkx22.R

./3. Analysis/Mutant/MutantPlot:
CKO.R
MOTIF_DLX.R
MOTIF_FOXD1.R
MOTIF_ISL1.R
MOTIF_LHX2.R

./3. Analysis/Mutant/P21Dlx:
Fig4_P21Dlx_ZI.R
Script_P21Dlx.R

./3. Analysis/Mutant/P8Dlx:
Fig4_P8Dlx.R
Script_P8Dlx.R

./3. Analysis/Mutant/Plot:
Mutant_SCENIC_Update.R

./3. Analysis/NPC:
Fig_NPC_Age_V2.R
NPC_Scenic.R

./3. Analysis/Plot:
Atlas/
NEUROPEPTIDE.R
Plotting.R
umap_contour_plot.R

./3. Analysis/Plot/Atlas:
Fig1_Script_Addition.R
Fig1_Script_scATAC_Server_Addition.R
Fig1_Script_scATAC_Server_Addition2.R
Fig1_scRNA_Update.R

./3. Analysis/Th:
merge_Th.R

./3. Analysis/Xenium:
Plot.py

./4. SCENICPlus:
1. RNA/
2. ATAC/
3. SCENIC/
Plotting/

./4. SCENICPlus/1. RNA:
Examples/
OBJECT_SCANPYPROCESSING.py
OBJECT_SEURAT2SCANPY.R
SCENIC_RNA_Preparation.sh*

./4. SCENICPlus/1. RNA/Examples:
Mutant/
NEUROGENESIS_SCANPYPROCESSING.py
NEUROGENESIS_SEURAT2SCANPY.R
PATTERN/
SCENIC_NEUROGENESIS_RNA.sh*

./4. SCENICPlus/1. RNA/Examples/Mutant:
MUTANT_SCANPYPROCESSING.py
MUTANT_SEURAT2SCANPY.R
SCENIC_MUTANT_RNA.sh*

./4. SCENICPlus/1. RNA/Examples/PATTERN:
PATTERN_SCANPYPROCESSING.py
PATTERN_SEURAT2SCANPY.R
SCENIC_PATTERN_RNA.sh*

./4. SCENICPlus/2. ATAC:
2A. pycisTOPIC/
2B. pycisTARGET/

./4. SCENICPlus/2. ATAC/2A. pycisTOPIC:
MUTANT/
Mastercopy/
NEUROGENESIS_PYCISTOPIC.py
NEUROGENESIS_SIGNAC2PYCISTOPIC.R
PATTERN/
SCENIC_NEUROGENESIS_ATAC.sh*

./4. SCENICPlus/2. ATAC/2A. pycisTOPIC/MUTANT:
DLX_PYCISTOPIC.py
FOXD1_PYCISTOPIC.py
ISL1_PYCISTOPIC.py
LHX2_PYCISTOPIC.py
MUTANT_SIGNAC2PYCISTOPIC.R
SCENIC_MUTANT_ATAC.sh*

./4. SCENICPlus/2. ATAC/2A. pycisTOPIC/Mastercopy:
OBJECT_PYCISTOPIC.py
OBJECT_SIGNAC2PYCISTOPIC.R
SCENIC_ATAC_pycisTOPIC.sh*

./4. SCENICPlus/2. ATAC/2A. pycisTOPIC/PATTERN:
PATTERN_PYCISTOPIC.py
PATTERN_SIGNAC2PYCISTOPIC.R
SCENIC_PATTERN_ATAC.sh*

./4. SCENICPlus/2. ATAC/2B. pycisTARGET:
Mastercopy/
Mutant/
NEUROGENESIS_Consensus.py
PATTERN/
SCENIC_NEUROGENESIS_CISTARGET.sh*

./4. SCENICPlus/2. ATAC/2B. pycisTARGET/Mastercopy:
OBJECT_Consensus.py
SCENIC_ATAC_pycisTARGET.sh*

./4. SCENICPlus/2. ATAC/2B. pycisTARGET/Mutant:
Dlx_Consensus.py
Foxd1_Consensus.py
Isl1_Consensus.py
Lhx2_Consensus.py
SCENIC_DLX_CISTARGET.sh*
SCENIC_FOXD1_CISTARGET.sh*
SCENIC_ISL1_CISTARGET.sh*
SCENIC_LHX2_CISTARGET.sh*

./4. SCENICPlus/2. ATAC/2B. pycisTARGET/PATTERN:
PATTERN_Consensus.py
SCENIC_PATTERN_CISTARGET.sh*

./4. SCENICPlus/3. SCENIC:
MUTANT/
Mastercopy/
Neurogenesis/
PATTERN/
Plotting/
SCENIC_NEUROGENESIS_CLI.sh*
chromsizes.tsv
genome_annotation.tsv
genome_annotation_with_chr.tsv
scenic things to fix.txt

./4. SCENICPlus/3. SCENIC/MUTANT:
SCENIC_DLX_CLI.sh*
SCENIC_FOXD1_CLI_V2.sh*
SCENIC_ISL1_CLI.sh*
SCENIC_LHX2_CLI.sh*

./4. SCENICPlus/3. SCENIC/Mastercopy:
SCENIC_DLX_Ctrl_CLI.sh*
chromsizes.tsv
genome_annotation.tsv
genome_annotation_with_chr.tsv

./4. SCENICPlus/3. SCENIC/Neurogenesis:
SCENIC_NEUROGENESIS_CLI.sh*

./4. SCENICPlus/3. SCENIC/PATTERN:
SCENIC_PATTERN_CLI.sh*

./4. SCENICPlus/3. SCENIC/Plotting:
SCENIC_Dlx.py
SCENIC_Foxd1.py
SCENIC_Isl1.py
SCENIC_Lhx2.py
SCENIC_Neurogenesis.py
SCENIC_Neurogenesis2.py
SCENIC_Pattern.py
SCENIC_Pattern.sh*
SCENIC_Pattern2.py
SCENIC_Pattern2.sh*
SCENIC_PatternX.sh*
SCENIC_Patternx.py

./4. SCENICPlus/Plotting:
SCENIC_Plot.py

./5. Xenium:
SpatialData.py
SpatialData_HET.py
SpatialData_HOMO.py
Squidpy.py
Squidpy_Plotting.py

./6. DeepLabCut:
Deeplabcut_Loca.py
Modified DeepLabCut/
Plot/
deeplabcut.sh*
deeplabcut2.sh*
deeplabcut3.sh*
deeplabcut_analyzing.py
deeplabcut_analyzing2.py
deeplabcut_create_extract.py
deeplabcut_evaluating.py
deeplabcut_training.py

./6. DeepLabCut/Modified DeepLabCut:
analyze_videos.py
trackingutils.py

./6. DeepLabCut/Plot:
DLCAnalyzer_Functions_final.R
Deeplabcut_ProcessingV2.R
Deeplabcut_Processing_LargeBatch.R



**RSession here**
> sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Copenhagen
tzcode source: system (glibc)

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] BSgenome.Mmusculus.UCSC.mm10_1.4.3 BSgenome_1.71.4                    rtracklayer_1.63.3                
 [4] BiocIO_1.13.1                      Biostrings_2.71.6                  XVector_0.43.1                    
 [7] TFBSTools_1.41.1                   JASPAR2020_0.99.10                 presto_1.0.0                      
[10] EnsDb.Mmusculus.v79_2.99.0         ensembldb_2.27.1                   AnnotationFilter_1.27.0           
[13] GenomicFeatures_1.55.4             AnnotationDbi_1.65.2               Biobase_2.63.1                    
[16] GenomicRanges_1.55.4               GenomeInfoDb_1.39.14               IRanges_2.37.1                    
[19] S4Vectors_0.41.7                   BiocGenerics_0.49.1                Signac_1.13.0                     
[22] pheatmap_1.0.12                    ggrepel_0.9.5                      ggraph_2.2.1                      
[25] igraph_2.0.3                       data.table_1.15.4                  colorspace_2.1-0                  
[28] reshape2_1.4.4                     stringr_1.5.1                      harmony_1.2.0                     
[31] Rcpp_1.0.12                        future_1.33.2                      patchwork_1.2.0                   
[34] ggplot2_3.5.1                      RColorBrewer_1.1-3                 Seurat_5.0.3                      
[37] SeuratObject_5.0.1                 sp_2.1-3                           Matrix_1.7-0                      
[40] dplyr_1.1.4                        cowplot_1.1.3                     

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.4.0               later_1.3.2                 bitops_1.0-7               
  [5] R.oo_1.26.0                 tibble_3.2.1                polyclip_1.10-6             DirichletMultinomial_1.45.0
  [9] XML_3.99-0.16.1             fastDummies_1.7.3           lifecycle_1.0.4             pwalign_0.99.2             
 [13] globals_0.16.3              lattice_0.22-6              MASS_7.3-60.2               magrittr_2.0.3             
 [17] plotly_4.10.4               yaml_2.3.8                  httpuv_1.6.15               sctransform_0.4.1          
 [21] spam_2.10-0                 spatstat.sparse_3.0-3       reticulate_1.36.1           CNEr_1.39.1                
 [25] pbapply_1.7-2               DBI_1.2.2                   abind_1.4-5                 zlibbioc_1.49.3            
 [29] Rtsne_0.17                  R.utils_2.12.3              purrr_1.0.2                 RCurl_1.98-1.14            
 [33] pracma_2.4.4                tweenr_2.0.3                GenomeInfoDbData_1.2.12     irlba_2.3.5.1              
 [37] listenv_0.9.1               spatstat.utils_3.0-5        seqLogo_1.69.0              goftest_1.2-3              
 [41] RSpectra_0.16-1             annotate_1.81.2             spatstat.random_3.2-3       fitdistrplus_1.1-11        
 [45] parallelly_1.37.1           DelayedArray_0.29.9         leiden_0.4.3.1              codetools_0.2-20           
 [49] RcppRoll_0.3.0              ggforce_0.4.2               tidyselect_1.2.1            UCSC.utils_0.99.7          
 [53] farver_2.1.1                viridis_0.6.5               matrixStats_1.3.0           spatstat.explore_3.2-7     
 [57] GenomicAlignments_1.39.5    jsonlite_1.8.8              tidygraph_1.3.1             progressr_0.14.0           
 [61] ggridges_0.5.6              survival_3.6-4              tools_4.4.0                 TFMPvalue_0.0.9            
 [65] ica_1.0-3                   glue_1.7.0                  SparseArray_1.3.5           gridExtra_2.3              
 [69] MatrixGenerics_1.15.1       withr_3.0.0                 fastmap_1.1.1               fansi_1.0.6                
 [73] caTools_1.18.2              digest_0.6.35               R6_2.5.1                    mime_0.12                  
 [77] GO.db_3.19.1                scattermore_1.2             poweRlaw_0.80.0             gtools_3.9.5               
 [81] tensor_1.5                  spatstat.data_3.0-4         RSQLite_2.3.6               R.methodsS3_1.8.2          
 [85] utf8_1.2.4                  tidyr_1.3.1                 generics_0.1.3              S4Arrays_1.3.7             
 [89] graphlayouts_1.1.1          httr_1.4.7                  htmlwidgets_1.6.4           uwot_0.2.2                 
 [93] pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                  lmtest_0.9-40              
 [97] htmltools_0.5.8.1           dotCall64_1.1-1             ProtGenerics_1.35.4         scales_1.3.0               
[101] png_0.1-8                   rstudioapi_0.16.0           tzdb_0.4.0                  rjson_0.2.21               
[105] nlme_3.1-164                curl_5.2.1                  cachem_1.0.8                zoo_1.8-12                 
[109] KernSmooth_2.23-22          parallel_4.4.0              miniUI_0.1.1.1              restfulr_0.0.15            
[113] pillar_1.9.0                vctrs_0.6.5                 RANN_2.6.1                  promises_1.3.0             
[117] xtable_1.8-4                cluster_2.1.6               readr_2.1.5                 cli_3.6.2                  
[121] compiler_4.4.0              Rsamtools_2.19.4            rlang_1.1.3                 crayon_1.5.2               
[125] future.apply_1.11.2         plyr_1.8.9                  stringi_1.8.3               viridisLite_0.4.2          
[129] deldir_2.0-4                BiocParallel_1.37.1         munsell_0.5.1               lazyeval_0.2.2             
[133] spatstat.geom_3.2-9         RcppHNSW_0.6.0              hms_1.1.3                   bit64_4.0.5                
[137] KEGGREST_1.43.1             shiny_1.8.1.1               SummarizedExperiment_1.33.3 ROCR_1.0-11                
[141] memoise_2.0.1               fastmatch_1.1-4             bit_4.0.5   


**Data Access
Raw Data (scRNA-seq)**
NCBI GEO: Raw matrix files, including count matrices for each sample and metadata, are accessible at NCBI GEO under accession number **GSEXXXXXX.**

**Processed Data (.Robj Files)**
Figshare: Processed Seurat objects (e.g., HyDD2_Processed.Robj) can be accessed and downloaded from Figshare at **Figshare link**. These files include pre-processed and normalized data suitable for running the scripts in this repository.

**Contact**
For questions or collaboration inquiries, please contact Thomas Kim at tkim@dandrite.au.dk
