# HyDD2
HyDD2: Analysis of Hypothalamic Development Database Version 2

**Project Description**
HyDD2 is an integrative bioinformatics project aimed at studying gene expression dynamics across multiple stages of hypothalamic development. This repository includes analysis scripts for processing, clustering, and visualizing single-cell RNA sequencing (scRNA-seq) data, focusing on the regulatory networks governing hypothalamic differentiation and maturation over key developmental stages.

This project is structured to support researchers interested in single-cell sequencing, developmental biology, and the underlying mechanisms of neurogenesis in the hypothalamus.

**#Repository Structure**
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



**#Mention RSession here**



**Data Access
Raw Data (scRNA-seq)**
NCBI GEO: Raw matrix files, including count matrices for each sample and metadata, are accessible at NCBI GEO under accession number **GSEXXXXXX.**

**Processed Data (.Robj Files)**
Figshare: Processed Seurat objects (e.g., HyDD2_Processed.Robj) can be accessed and downloaded from Figshare at **Figshare link**. These files include pre-processed and normalized data suitable for running the scripts in this repository.

**Contact**
For questions or collaboration inquiries, please contact Thomas Kim at tkim@dandrite.au.dk
