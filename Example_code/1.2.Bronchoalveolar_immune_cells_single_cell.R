setwd("downloaded/deCS-master")

# Load the preprocessed bronchoalveolar immune cells dataset
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_single_cell_marker.rda")
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_single_cell_top2k.rda")
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_cell_info.rda")

dim(BIC_single_cell_marker)
dim(BIC_single_cell_top2k)
dim(cell_info)

# deCS annotation
library(deCS)
data(BlueprintEncode_main)
data(Human_cell_landscape)

Sys.time()	# Time costs ~ 7 seconds
BIC_deCS_anno_marker = deCS.correlation(BIC_single_cell_marker, BlueprintEncode_main_t_score, methods = "pearson", top_n = 3)
Sys.time()	# Time costs ~ 11 seconds
BIC_deCS_anno_top2k = deCS.correlation(BIC_single_cell_top2k, BlueprintEncode_main_t_score, methods = "pearson", top_n = 3)
Sys.time()

BIC_deCS_summary_marker <- table(BIC_deCS_anno_marker$Cell_labels, cell_info$celltype)
BIC_deCS_summary_top2k <- table(BIC_deCS_anno_top2k$Cell_labels, cell_info$celltype)
BIC_deCS_summary_marker
BIC_deCS_summary_top2k

# SingleR annotation
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_single_cell_marker.rda")
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_single_cell_top2k.rda")
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_cell_info.rda")

library(SingleR)
library(celldex)
ref.BlueprintEncode <- celldex::BlueprintEncodeData()

Sys.time()	# Time costs ~ 13 seconds
BIC_singleR_anno <- SingleR(test = BIC_single_cell_marker, ref = ref.BlueprintEncode, assay.type.test = 1, labels = ref.BlueprintEncode$label.main)
Sys.time()	# Time costs ~ 48 seconds
BIC_singleR_anno_top2k <- SingleR(test = BIC_single_cell_top2k, ref = ref.BlueprintEncode, assay.type.test = 1, labels = ref.BlueprintEncode$label.main)
Sys.time()

BIC_singleR_summary <- table(BIC_singleR_anno$pruned.labels, cell_info$celltype)
BIC_singleR_summary_top2k <- table(BIC_singleR_anno_top2k$pruned.labels, cell_info$celltype)

