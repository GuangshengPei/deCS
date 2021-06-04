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
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Manually annotation
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# deCS annotation
library(deCS)
#  Method A: Correlation analysis
pbmc_cluster_average = AverageExpression(pbmc)
pbmc_cluster_marker_average = pbmc_cluster_average[[1]][which(rownames(pbmc_cluster_average[[1]]) %in% top10$gene), ]
pbmc_cluster_marker_z_score = t(scale(t(pbmc_cluster_marker_average)))

data(MonacoImmune_main)
data(MonacoImmune_fine)
pbmc_Monaco_main_cor = deCS.correlation(pbmc_cluster_marker_z_score, MonacoImmune_main_t_score, methods = "pearson")
pbmc_Monaco_fine_cor = deCS.correlation(pbmc_cluster_marker_z_score, MonacoImmune_fine_t_score, methods = "pearson")

#Method B: Fisher's exact test
pbmc_top10_markers_list = pbmc.markers[which(pbmc.markers$gene %in% top10$gene),] 
pbmc_top10_markers_list$cluster <- as.character(pbmc_top10_markers_list$cluster)

for (i in 1:length(unique(pbmc_top10_markers_list$cluster))){
	pbmc_top10_markers_list$cluster[pbmc_top10_markers_list$cluster == as.character(i-1)] <- as.character(new.cluster.ids[i])
}

data(MonacoImmune_main)
data(MonacoImmune_fine)
pbmc_Monaco_main_fisher = deCS.fisher(pbmc_top10_markers_list, MonacoImmune_main_t_score, ref_panel_markers_ratio = 0.05, p_threshold = 0.05)
pbmc_Monaco_fine_fisher = deCS.fisher(pbmc_top10_markers_list, MonacoImmune_fine_t_score, ref_panel_markers_ratio = 0.05, p_threshold = 0.05)

# Users can also upload custom cell type marker genes list as new reference
# Take CellMatch Database as example
data(CellMatch)
head(CellMatch_markers, 200)

# Intersection ratio
pbmc_Cellmatch_fisher_IR <- deCS.fisher(pbmc_top10_markers_list, CellMatch_markers, type = "list", p.adjust.methods = "bonferroni", p_threshold = 1e-3, cell_type_threshold = 0.05)

# Intersection count
pbmc_Cellmatch_fisher_IC <- deCS.fisher(pbmc_top10_markers_list, CellMatch_markers, type = "list", label = "count", p.adjust.methods = "bonferroni", p_threshold = 1e-3, , cell_type_threshold = 1)

# SingleR annotation
library(SingleR)
library(celldex)
ref1.MonacoImmune <- celldex::MonacoImmuneData()

pbmc_singleR <- SingleR(test = pbmc_cluster_marker_z_score, ref = ref1.MonacoImmune, assay.type.test = 1, labels = ref1.MonacoImmune$label.main)

