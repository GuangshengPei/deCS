# 1. deCS introduction
&#8194;&#8194;Single-cell RNA sequencing (scRNA-seq) is rapidly accelerate our understanding of the cellular compositions of complex tissues. Yet, one major limitation for current protocols rely on manual annotations, which are subjectivity and time-consuming. The increasing numbers of scRNA-seq data sets, as well as numerous genetic studies, enable us to build a comprehensive cell type reference atlas. Here, we present deCS, for automatic cell type annotations based on a comprehensive collection of human cell type expression profiles or list of marker genes. We applied deCS to single-cell data sets from various tissues, and systematically evaluated the annotation accuracy under different conditions. Under the same conditions, deCS runs faster and have comparable even better accuracy than the competitive tools.
# 2. Usage
## 2.1 Installing deCS
### Requirements of other dependencies
deCS relies on R (>= 3.5), reshape2 (>= 1.4.4), ggplot2 (>= 3.3.2). Please follow their installation instruction.   

```
install.packages("reshape2")
install.packages("ggplot2")
```
### deCS R package can be easily installed from Github using devtools:
```
# install.packages("devtools")        
devtools::install_github("GuangshengPei/deCS")
```
### Getting Started 
Once we have the package installed, we can load the package. 
```
library(deCS) 
``` 
## 2.2 Built-in data loading
&#8194;&#8194;deCS collected several cell type reference panels, including BlueprintEncode, the Database of Immune Cell Expression (DICE), MonacoImmune, human cell landscape, human cell atlas of fetal et al. After installation of deCS package, one can load the build-in references using the following commands:  
```
# BlueprintEncode_main_t_score
data(BlueprintEncode_main)
# BlueprintEncode_fine_t_score
data(BlueprintEncode_fine)
# DICE_main_t_score
data(DICE_main)
# DICE_fine_t_score
data(DICE_fine)
# MonacoImmune_main_t_score
data(MonacoImmune_main)
# MonacoImmune_fine_t_score
data(MonacoImmune_fine)
# Human_cell_landscape (HCL_z_score)
data(Human_cell_landscape)
# Human_cell_atlas_of_fetal (HCAF_z_score)
data(Human_cell_atlas_of_fetal)
```
&#8194;&#8194;In addition, one can also load the cell type marker gene list from CellMatch database. 
```
data(CellMatch)
head(CellMatch_markers)
                       Cell_type Marker_gene
1 1-Cell Stage Cell (Blastomere)       ACCSL
2 1-Cell Stage Cell (Blastomere)      ACVR1B
3 1-Cell Stage Cell (Blastomere)    ARHGEF16
4 1-Cell Stage Cell (Blastomere)       ASF1B
5 1-Cell Stage Cell (Blastomere)     BCL2L10
6 1-Cell Stage Cell (Blastomere)       BLCAP
```
&#8194;&#8194;The cell type-marker genes list of `"CellMatch_markers"` will be loaded.
## 2.3 Input data
&#8194;&#8194;Depending on the type of query data, we implemented two test approaches: Correlation analysis and Fisher's exact test for cell type enrichment analysis.
In this tutorial, we will run deCS on [preprocessed PBMC data](https://github.com/GuangshengPei/deCS/tree/main/Example_data/1.1.PBMC/pbmc_example.rda). Or start with Cellranger output.

### Load the PBMC dataset
```
library(Seurat)
library(dplyr)
pbmc.data <- Read10X(data.dir = "Example_data/1.1.PBMC/hg19/")
```
### Initialize the Seurat object with the raw data, then conduct standard normalization.
```
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() 
```
### Standard cell clustering analysis.
```
pbmc <- FindNeighbors(pbmc, dims = 1:10) 
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10) 
```
### Identify cluster specific expressed genes.
```
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pbmc_top10_markers_list = pbmc.markers[which(pbmc.markers$gene %in% top10$gene),] 
```
### Union of marker genes and z-score calculation.
```
pbmc_cluster_average = log(AverageExpression(pbmc)[[1]] + 1, 2)
pbmc_cluster_marker_average = pbmc_cluster_average[which(rownames(pbmc_cluster_average) %in% top10$gene), ]
pbmc_cluster_marker_z_score = t(scale(t(pbmc_cluster_marker_average)))
```
### 2.3.1 deCS correlation analysis for expression profiles
&#8194;&#8194;If the query is gene expression profile, we provide function `deCS.correlation()` for cell type enrichment analysis. Simply, we calculate pearson correlation coefficient (PCC) or Spearman's rank correlation coefficient with each of the cell type in the reference dataset, then the most relevant cell type(s) will be identified. 
```
load("pbmc_example.rda")
head(pbmc_cluster_marker_z_score)
       Naive CD4 T Memory CD4 T CD14+ Mono          B      CD8 T FCGR3A+ Mono          NK         DC   Platelet
RHOC    -0.7581827   -0.6728063 -0.5539563 -0.6884791 -0.3323893    2.2580565  0.77225930 -0.4513599  0.4268578
CD2      0.4694667    1.7349728 -0.7898039 -0.7591715  1.4013552   -0.7432280  0.09901601 -0.5552728 -0.8573346
S100A9  -0.4071876   -0.4062238  2.6504337 -0.4045475 -0.4059206   -0.1321780 -0.40858303 -0.1606083 -0.3251849
S100A8  -0.3739926   -0.3716080  2.6617595 -0.3746993 -0.3819118   -0.2019757 -0.37228956 -0.2841058 -0.3011767
FCER1A  -0.3220375   -0.3432463 -0.3138008 -0.3430898 -0.3326199   -0.3414830 -0.32525747  2.6665050 -0.3449702
FCER1G  -0.7553841   -0.7614134  0.6567761 -0.7648104 -0.7323036    2.2576395  0.30423366  0.1057573 -0.3104950

# deCS.correlation(markers_expression, ref_panel)
# pdf ("deCS_Cor_annotation_A.pdf", 8, 8)
pbmc_deCS_cor_panel_A <- deCS.correlation(pbmc_cluster_marker_z_score, MonacoImmune_main_t_score)
# dev.off()
```  
&#8194;&#8194;Here, `markers_expression` is scaled gene expression matrix for human scRNA-seq or bulk RNA-seq, users can see detail preprocess steps from [Example_code page](https://github.com/GuangshengPei/deCS/tree/main/Example_code); `ref_panel` is pre-calculated cell type specificity score reference panel, see section (2.2 Built-in data loading).  
```  
pbmc_deCS_cor_panel_A 
         Query Max_correlation            Top1            Top2            Top3     Cell_labels
1  Naive CD4 T       0.8670209    CD4+ T cells    CD8+ T cells         T cells    CD4+ T cells
2 Memory CD4 T       0.8621523    CD4+ T cells    CD8+ T cells         T cells    CD4+ T cells
3   CD14+ Mono       0.7053299       Monocytes     Neutrophils Dendritic cells       Monocytes
4            B       0.8152868         B cells Dendritic cells         T cells         B cells
5        CD8 T       0.7981408    CD8+ T cells         T cells    CD4+ T cells    CD8+ T cells
6 FCGR3A+ Mono       0.7660590       Monocytes     Neutrophils Dendritic cells       Monocytes
7           NK       0.8930805        NK cells         T cells    CD8+ T cells        NK cells
8           DC       0.8151785 Dendritic cells       Monocytes         B cells Dendritic cells
9     Platelet       0.7281423     Progenitors       Basophils     Neutrophils     Progenitors

write.table(pbmc_deCS_cor_panel_A, file = "pbmc_deCS_result.txt", sep = "\t", quote = F)
```  
Users can change other reference panels and parameters, e.g. give `cor_threshold` to avoid mis-annotation (regard as "Undetermined cells").
```  
 deCS.correlation(pbmc_cluster_marker_z_score, HCL_z_score, top_n = 3, cor_threshold = 0.5, p_threshold = 0.01, cell_type_threshold = 0.5)
         Query Max_correlation                     Top1                   Top2                           Top3              Cell_labels
1  Naive CD4 T       0.6355362 C52_Proliferating T cell              C6_T cell          C77_Ureteric bud cell C52_Proliferating T cell
2 Memory CD4 T       0.6794370 C52_Proliferating T cell              C6_T cell C64_Fetal skeletal muscle cell C52_Proliferating T cell
3   CD14+ Mono       0.6579867             C13_Monocyte         C51_Macrophage                 C69_Macrophage             C13_Monocyte
4            B       0.8321910               C37_B cell C3_B cell (Plasmocyte)       C9_Enterocyte progenitor               C37_B cell
5        CD8 T       0.7537835               C24_T cell       C93_Myeloid cell                      C6_T cell               C24_T cell
6 FCGR3A+ Mono       0.4364741           C69_Macrophage           C13_Monocyte                  C2_Macrophage       Undetermined cells
7           NK       0.6104342               C24_T cell       C93_Myeloid cell         C35_Smooth muscle cell               C24_T cell
8           DC       0.7169097       C22_Dendritic cell           C13_Monocyte                  C2_Macrophage       C22_Dendritic cell
9     Platelet       0.4415510     C20_Endothelial cell   C76_Mesothelial cell               C72_Stromal cell       Undetermined cells
```  
More parameters in `deCS.correlation` function is available at `help(deCS.correlation)`.   

### 2.3.2 deCS Fisher's exact test for list of genes     
&#8194;&#8194;If the query is a list of genes (e.g. union of marker genes, traits associated genes), we provide function `deCS.fisher()`, implement with Fisher’s exact test to identify query gene set enriched in cell type-specific genes (CTgenes). We allow the user to define the cutoff values, e.g., the top 5% genes with highest t-score/z-score as CTgenes. For each query gene set and CTgenes in a given cell type, deCS will identify whether a set of “candidate marker genes” are disproportionately overlapped in a specific cell type specific genes.     
```  
load("pbmc_example.rda")
head(pbmc_top10_markers_list, 20)
                  p_val avg_logFC pct.1 pct.2     p_val_adj      cluster      gene
RPS3A     1.075148e-108 0.5375120 0.991 0.977 1.474457e-104  Naive CD4 T     RPS3A
LDHB      1.963031e-107 0.7300635 0.901 0.594 2.692101e-103  Naive CD4 T      LDHB
CCR7       1.606796e-82 0.9219135 0.436 0.110  2.203560e-78  Naive CD4 T      CCR7
CD3D       4.198081e-77 0.6598608 0.838 0.406  5.757249e-73  Naive CD4 T      CD3D
CD3E       2.316054e-54 0.5999356 0.726 0.399  3.176237e-50  Naive CD4 T      CD3E
NOSIP      3.191555e-50 0.6926087 0.628 0.358  4.376899e-46  Naive CD4 T     NOSIP
LEF1       3.324866e-49 0.7296372 0.336 0.104  4.559722e-45  Naive CD4 T      LEF1
PRKCQ-AS1  2.498572e-44 0.7119736 0.331 0.110  3.426542e-40  Naive CD4 T PRKCQ-AS1
PIK3IP1    7.614450e-43 0.6500625 0.438 0.185  1.044246e-38  Naive CD4 T   PIK3IP1
IL7R       3.894244e-35 0.5029688 0.597 0.333  5.340566e-31  Naive CD4 T      IL7R
MAL        6.485255e-33 0.6420019 0.263 0.088  8.893879e-29  Naive CD4 T       MAL
TRAT1      4.014940e-15 0.3657010 0.265 0.140  5.506089e-11  Naive CD4 T     TRAT1
IL32       1.894810e-92 0.8373872 0.948 0.464  2.598542e-88 Memory CD4 T      IL32
LTB        7.953303e-89 0.8921170 0.981 0.642  1.090716e-84 Memory CD4 T       LTB
CD3D1      1.655937e-70 0.6436286 0.919 0.431  2.270951e-66 Memory CD4 T      CD3D
IL7R1      3.688893e-68 0.8147082 0.747 0.325  5.058947e-64 Memory CD4 T      IL7R
LDHB1      2.292819e-67 0.6253110 0.950 0.613  3.144372e-63 Memory CD4 T      LDHB
CD2        2.504468e-61 0.8559204 0.652 0.244  3.434627e-57 Memory CD4 T       CD2
AQP3       1.851623e-60 0.8586034 0.422 0.110  2.539316e-56 Memory CD4 T      AQP3
CD3E1      8.015029e-55 0.6006175 0.828 0.409  1.099181e-50 Memory CD4 T      CD3E

# deCS.fisher(markers_list, ref_panel)
pbmc_deCS_FET_panel_A <- deCS.fisher(pbmc_top10_markers_list, MonacoImmune_main_t_score)
```  
Here, `markers_list` is a gene list table with at least two columns (named with "cluster" and "gene"); `ref_panel` is pre-calculated cell type specificity score reference panel. Please specifiy type = "list" when the reference is cell type-marker gene lists. 
``` 
pbmc_deCS_FET_CellMatch <- deCS.fisher(pbmc_top10_markers_list, CellMatch_markers, type = "list", p.adjust.methods = "bonferroni", p_threshold = 1e-3, cell_type_threshold = 0.05)
``` 
#### Users can also create your own cell type-marker genes list, with at least two columns with names `Cell_type` and `Marker_gene`.  
More parameters in `deCS.fisher` function is available at `help(deCS.fisher)`.      
## 3. Shiny application
&#8194;&#8194;We further incorporated a shiny application to deCS package available at https://gpei.shinyapps.io/decs_cor/ or https://gpei.shinyapps.io/decs_fisher/. Users only need to upload a well clustered Seurat object (save by saveRDS function). Please see details in demo file (Shiny_app/pbmc_downsampled.rds). Within one minute, the deCS annotation result files among different references will be displayed in an interactive mode, users can change any parameters to optimize their results.     
## 4. More examples application and evaluation   
## 4.1 single cell RNA-seq data
&#8194;&#8194;To explore the deCS utility for direct annotation of single cell expression profiles, we evaluated three different scRNA-seq datasets.    
&#8194;&#8194;(1) The test dataset of 3k Peripheral Blood Mononuclear Cells (PBMC) was collected from [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/datasets). The raw data of 3k PBMCs from a healthy donor included 2,700 single cells that were sequenced on the Illumina NextSeq 500 platform. Overall this dataset with ~2200 median UMI count per cell and ~800 median genes per cell.    
&#8194;&#8194;(2) The test dataset of bronchoalveolar immune cells was collected from healthy controls (HC, n = 3) and patients with moderate (M, n = 3) and severe (S, n = 6) COVID-19 infection. The test data included 63,103 single cells sequenced on a BGI MGISEQ-2000 or Illumina platform. The raw data and predefined label are available from [Github](https://github.com/zhangzlab/covid_balf).   
&#8194;&#8194;(3) The cardiac cells scRNA-seq data were collected from 18 human embryos, ranging from 5 weeks (5W) to 25W of gestation. The UMI count of 4,948 cells and predefined label were downloaded from NCBI GEO GSE106118. Each single cell detected 4,109 genes on average.   
&#8194;&#8194;All data preprocess scripts and datasets are available at https://github.com/GuangshengPei/deCS/tree/master/Example_code.
## 4.2 Bulk RNA-seq 
&#8194;&#8194;As an extension of our previous work, deCS can also be applied on bulk RNA-seq data. Human induced pluripotent stem cell (iPSC) have revolutionized the study of the biological mechanisms of psychiatric disorders, as they allow for the establishment of brain cellular models that account for a patient’s genetic background. In most iPSC derived cells studies, a traditional routine is to assess their differentiation quality, so that the resulting iPSCs can be transitioned towards clinical applications effectively.   
&#8194;&#8194;(1) The first bulk RNA-seq data of schizophrenia (SCZ) associated human induced pluripotent stem cell (hiPSC)-derived cell lines were generated from the population isolate of the Central Valley of Costa Rica (CVCR), combined with three unrelated individuals from [Stertz el al. study](https://www.nature.com/articles/s41386-020-00924-0). One clone from each subject was differentiated into Neuronal Precursor Cells (NPC), and neuron. In total, RNA-seq from hiPSC-NPCs (n=13) and hiPSC-neurons (n=11) were collected.   
&#8194;&#8194;(2) The second hiPSC-derived cell lines bulk RNA-seq data were collected from Sry-box 9 (SOX9+) chondroprogenitors (n=4) and GDF5+ mesenchymal cells (n=4). SOX9 is the master transcription factor for chondrocytes. However, the SOX9+ cells tended to form transient cartilage that was readily mineralized in vivo. In contract, GDF5+ cells preferentially formed permanent cartilage that expressed low-to-no hypertrophic chondrocyte markers in vitro and remained unmineralized or only became partially mineralized in vivo.   
&#8194;&#8194;We provide data preprocess scripts at https://github.com/GuangshengPei/deCS/tree/master/Example_code.
## 4.3 Traits associated genes from GWAS summary data
&#8194;&#8194;Deeper understanding of causal tissues of human complex diseases is an important step towards the etiology of disease origin, yet tissues are complex milieus consisting of numerous cell types. Tissue level association failed to elucidate cell type contributions in disease. The aim of this application was to illustrate the association between cell type and disease. To this end, deCS was applied to the preprocessed [GWAS data](https://academic.oup.com/nar/article/49/1/53/6029182) using the model for Fisher’s exact test. We provide data preprocess scripts at https://github.com/GuangshengPei/deCS/tree/master/Example_code/3.GWAS_trait_associate_genes.R.  
## Gene symbol transformation and mouse application
&#8194;&#8194;deCS works only with human gene symbol, users should transform the human ensembl ids to gene symbol at first. For mouse data, users can use the ortholog genes, or just uppercase the query gene names by `toupper()` function.   
## System Requirements
&#8194;&#8194;Hardware requirements: deCS package requires only a standard laptop with enough RAM to support the in-memory operations. deCS package is supported for Windows,  macOS and Linux. deCS can be installed on a normal computer within few mins.
## Help
&#8194;&#8194;If you have any question, comment or suggestion, please contact peiguangsheng@gmail.com.


