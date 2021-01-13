# 1. deCS introduction (v1.0)
&#8194;&#8194;Single-cell RNA sequencing (scRNA-seq) is rapidly accelerate our understanding of the cellular compositions of complex tissues. Yet, one major limitation for current protocols rely on manual annotations, which are subjectivity and time-consuming. Existing cell automatically annotation methods strongly depend on the prior knowledge. With the exponential growth number of scRNA-seq data, allows us to create a comprehensive cell type specific expression profiles. On the other hands, such profiles can also be aligned with genome-wide association studies (GWASs) to implicate cell type specificity of the traits. Here, we collected several comprehensive resource of cell type expression profiles or marker genes for hundreds of cell types in human tissues. By extending our deTS algorithm, we implemented different statistical methods for detecting cell specific genes, and enrichment test of different forms of query data. Using different benchmark data sets, we demonstrate the feasibility of deCS can be widely applied on scRNA cell type definition, and outperformed SingleR on COVID-19 infection patients. In addition, we demonstrated the widespread utility of deCS on iPSC-derived bulk RNA-seq differentiation efficiency determination, and elucidation of trait-cell type associations for human complex diseases, thus providing better insight into mechanism underlying diseases pathogenesis.
# 2. Usage
## 2.1 Installing deCS
### Requirements of other dependencies
&#8194;&#8194;deCS relies on R (>= 3.5), reshape2 (>= 1.4.4), ggplot2 (>= 3.3.2)  
&#8194;&#8194;Please follow their installation instruction.  
&#8194;&#8194;Install reshape2 using `> install.packages("reshape2")  `.  
&#8194;&#8194;Install ggplot2 using`> install.packages("ggplot2")  `.
### deCS R package can be easily installed from Github using devtools:
&#8194;&#8194;`devtools::install_github("GuangshengPei/deCS")`  
### deCS loading
&#8194;&#8194;Load the deCS package   
&#8194;&#8194;`> library(deCS)`  
 
## 2.2 Built-in data loading
&#8194;&#8194;deCS requires different cell type reference panels to conduct the cell type enrichment test.    
&#8194;&#8194;Here, `ref_panel = c("BlueprintEncode_main", "BlueprintEncode_fine", "Human_cell_landscape", "DICE_main", "DICE_fine", "MonacoImmune_main", "MonacoImmune_fine", "CellMatch")`.    
&#8194;&#8194;After installation of deCS package, one can load the reference data using the following commands:  

&#8194;&#8194;&#8194;&#8194;`> data(BlueprintEncode_main)`  
&#8194;&#8194;&#8194;&#8194;`> data(BlueprintEncode_fine)`  
&#8194;&#8194;&#8194;&#8194;`> data(Human_cell_landscape)`  
&#8194;&#8194;&#8194;&#8194;`> data(DICE_main)`  
&#8194;&#8194;&#8194;&#8194;`> data(DICE_fine)`  
&#8194;&#8194;&#8194;&#8194;`> data(MonacoImmune_main)`  
&#8194;&#8194;&#8194;&#8194;`> data(MonacoImmune_fine)`  
&#8194;&#8194;Or load the cell type marker gene list:  
&#8194;&#8194;&#8194;&#8194;`> data(CellMatch)`  
&#8194;&#8194;Then `"BlueprintEncode_main_t_score"`, `"BlueprintEncode_fine_t_score"`, `"HCL_z_score"`, `"DICE_main_t_score"`, `"DICE_fine_t_score"`, `"MonacoImmune_main_t_score"`, `"MonacoImmune_fine_t_score"`,  and `"CellMatch_markers"` will be loaded to R enviroment.  
&#8194;&#8194;Users can also upload their cell type marker gene list or database, like reference panel `CellMatch`, with at least two columns names with `c("Cell_type", "Marker_gene")`. 
## 2.3 Input data
&#8194;&#8194;Depending on the type of query data, we implemented two test approaches: Correlation analysis and Fisher's exact test for cell type enrichment analysis.
### 2.3.1 deCS.correlation() for expression profiles
&#8194;&#8194;If the query is gene expression profile, we provide function `deCS.correlation()` for cell type enrichment analysis. Simply, we calculate pearson correlation coefficient (PCC) or Spearman's rank correlation coefficient with each of the cell type in the reference dataset, then the most relevant cell type(s) will be identified.   
&#8194;&#8194;`> deCS.correlation(markers_expression, ref_panel) `     
&#8194;&#8194;Here, `markers_expression` is scaled gene expression matrix for human scRNA-seq or bulk RNA-seq.    
&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;`ref_panel` is pre-calculated cell type specificity score reference panel.     
&#8194;&#8194;More parameters in `deCS.correlation` function is available at `help(deCS.correlation)`.     
#### Gene expression profiles normalization
&#8194;&#8194;Before `deCS.correlation` analysis, make sure raw expression profiles have been normalized so that each cell in query data will be scaled appropriately with the reference data. We provided two normalization approaches: `"z-score"` and `"abundance"`, which available from our previous deTS package `tsea.expression.normalization()` function:      
&#8194;&#8194;(1) `"z-score"` normalization will calculate a z-score for the query sample for each cell type in the reference panel as below: e_i=(e_0-μ_t))/sd_t, where μ_t and sd_t were the mean and SD of cell type t.      
&#8194;&#8194;(2) `"abundance"` normalization will provide an abundance correction approach for the query sample for each cell type in the reference panel as below: e_i=(log2(e_0+1)/(log2(u_t+1)+1).     

### 2.3.2 deCS.fisher() for list of genes     
&#8194;&#8194;If the query is a list of genes (e.g. union of marker genes, traits associated genes), we provide function `deCS.fisher()`, implement with Fisher’s exact test to identify query gene set enriched in cell type specific genes (CSGs). We allow the user to define the cutoff values, e.g., the top 5% genes as CSGs. For each query gene set and CSGs in a given cell type, deCS will identify whether a set of “candidate genes” of interest are disproportionately expressed in a specific cell type by using Fisher’s exact test.     
&#8194;&#8194;`> deCS.fisher(markers_list, ref_panel) `        
&#8194;&#8194;Here, `markers_list` is a gene list table with at least two columns names with "cluster" and "gene".    
&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;&#8194;`ref_panel` is pre-calculated cell type specificity score reference panel.     
&#8194;&#8194;More parameters in `deCS.fisher` function is available at `help(deCS.fisher)`.      
## 3. Examples application   
## 3.1 single cell RNA-seq data
&#8194;&#8194;To explore the deCS utility for direct annotation of single cell expression profiles, we attempted to analyze three different scRNA-seq datasets.    
&#8194;&#8194;(1) The test dataset of 3k Peripheral Blood Mononuclear Cells (PBMC) was collected from 10X Genomics website (https://support.10xgenomics.com/single-cell-gene-expression/datasets). The raw data of 3k PBMCs from a healthy donor included 2,700 single cells that were sequenced on the Illumina NextSeq 500 platform. Overall this dataset with ~2200 median UMI count per cell and ~800 median genes per cell.    
&#8194;&#8194;(2) The test dataset of bronchoalveolar immune cells was collected from healthy controls (HC, n = 3) and patients with moderate (M, n = 3) and severe (S, n = 6) COVID-19 infection. The test data included 63,103 single cells sequenced on a BGI MGISEQ-2000 or Illumina platform. The raw data and predefined label are available from https://github.com/zhangzlab/covid_balf. Overall, these datasets with ~3300, 2638, 2400 median UMI count per cell and ~1100, 144, 1000 median genes per cell for HC, M and S groups, respectively.   
&#8194;&#8194;(3) The cardiac cells scRNA-seq data were collected from 18 human embryos, ranging from 5 weeks (5W) to 25W of gestation. The UMI count of 4,948 cells and predefined label were downloaded from NCBI GEO GSE106118. Each single cell detected 4,109 genes on average.   
&#8194;&#8194;After getting above single-cell gene expression UMI matrices, we used the Seurat package for downstream unsupervised clustering and highly variable genes identification. The Seurat `FindAllMarkers()` and `AverageExpression()` functions were used to identify cluster-specific marker genes, and gene average expression of each cluster. We provide data preprocess scripts at https://github.com/GuangshengPei/deCS/tree/master/Example_code.
## 3.2 Bulk RNA-seq 
&#8194;&#8194;As an extension of our previous work, deCS can also be applied on bulk RNA-seq data. Human induced pluripotent stem cell (iPSC) have revolutionized the study of the biological mechanisms of psychiatric disorders, as they allow for the establishment of brain cellular models that account for a patient’s genetic background. In most iPSC derived cells studies, a traditional routine is to assess their differentiation quality, e.g. compare their expression level similarity to others, so that the resulting iPSCs can be transitioned towards clinical applications effectively.   
&#8194;&#8194;(1) The first bulk RNA-seq data of schizophrenia (SCZ) associated human induced pluripotent stem cell (hiPSC)-derived cell lines were generated from the population isolate of the Central Valley of Costa Rica (CVCR), combined with three unrelated individuals [Laura Stertz et al.]. One clone from each subject was differentiated into Neuronal Precursor Cells (NPC), and neuron. In total, RNA-seq from hiPSC-NPCs (n=13) and hiPSC-neurons (n=11) were collected.   
&#8194;&#8194;(2) The second hiPSC-derived cell lines bulk RNA-seq data were collected from Sry-box 9 (SOX9+) chondroprogenitors (n=4) and GDF5+ mesenchymal cells (n=4) [Azim Pothiawala et al.]. SOX9 is the master transcription factor for chondrocytes. However, the SOX9+ cells tended to form transient cartilage that was readily mineralized in vivo. In contract, GDF5+ cells preferentially formed permanent cartilage that expressed low-to-no hypertrophic chondrocyte markers in vitro and remained unmineralized or only became partially mineralized in vivo.   
&#8194;&#8194;For above bulk RNA-seq data sets, an average of 49.8 and 25.1 million clean reads per sample were generated, and then alignment to GRCh38 with STAR. Expression of the genes were normalized as TPM (Transcripts Per Kilobase Million) for accurate quantification using RSEM software. We provide data preprocess scripts at https://github.com/GuangshengPei/deCS/tree/master/Example_code.
## 3.3 Traits associated genes from GWAS summary data
&#8194;&#8194;Deeper understanding of causal tissues of human complex diseases is an important step towards the etiology of disease origin, yet tissues are complex milieus consisting of numerous cell types. Tissue level association failed to elucidate cell type contributions in disease. The aim of this application was to illustrate the association between cell type and disease. To this end, deCS was applied to the GWAS data using the model for Fisher’s exact test. We provide data preprocess scripts at https://github.com/GuangshengPei/deCS/tree/master/Example_code/3.GWAS_trait_associate_genes.R.  
## Citation
Cell type enrichment analysis for scRNA-seq cell type annotation and elucidation of trait - cell associations for human complex traits.

