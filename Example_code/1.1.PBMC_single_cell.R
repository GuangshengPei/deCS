setwd("downloaded/deCS-master")

# Load Seurat package
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "Example_data/1.1.PBMC/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Single cell RNA-seq data preprocessing
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Identify cluster specific expressed genes
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_lo2gFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Manually annotation
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Data normalization 
cell_count = pbmc@assays[[1]]@counts
cell_expression_scale = pbmc@assays[[1]]@scale.data
cell_expression_marker_scale = cell_expression_scale[which(rownames(cell_expression_scale) %in% top10$gene), ]
cell_info = pbmc@meta.data

save(cell_count, cell_expression_scale, cell_expression_marker_scale, cell_info, file = "PBMC_single_cell.rda")

# deCS single cell annotation
# load("PBMC_single_cell.rda")
library(deCS)
data(MonacoImmune_main)
data(MonacoImmune_fine)

# deCS correlation annotation by using all genes
Sys.time()	# Time costs ~ 1 second
pbmc_deCS_anno = deCS.correlation(cell_expression_scale, MonacoImmune_main_t_score, methods = "pearson", cor_threshold = 0)
Sys.time()	

pbmc_deCS_summary = table(pbmc_deCS_anno$Cell_labels, cell_info$seurat_clusters)
colnames(pbmc_deCS_summary) <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
pbmc_deCS_summary

# deCS correlation annotation by only using cluster specific expressed genes
Sys.time()	# Time costs ~ 1 second
pbmc_deCS_anno_marker = deCS.correlation(cell_expression_marker_scale, MonacoImmune_main_t_score, methods = "pearson", cor_threshold = 0)
Sys.time()

pbmc_deCS_summary_marker = table(pbmc_deCS_anno_marker$Cell_labels, cell_info$seurat_clusters)
colnames(pbmc_deCS_summary_marker) <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
pbmc_deCS_summary_marker

# SingleR single cell annotation
# load("PBMC_single_cell.rda")
library(SingleR)
library(celldex)
ref.MonacoImmune <- celldex::MonacoImmuneData()

Sys.time()	# Time costs ~ 7 seconds
pbmc_singleR_anno <- SingleR(test = cell_count, ref = ref.MonacoImmune, assay.type.test = 1, labels = ref.MonacoImmune$label.main)
Sys.time()

pbmc_singleR_summary = table(pbmc_singleR_anno$pruned.labels, cell_info$seurat_clusters)
colnames(pbmc_singleR_summary) <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
pbmc_singleR_summary
