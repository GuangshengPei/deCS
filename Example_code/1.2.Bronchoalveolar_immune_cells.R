setwd("downloaded/deCS-master")

# Load the preprocessed bronchoalveolar immune cells dataset
load("Example_data/1.2.Bronchoalveolar_immune_cells/BIC_cluster.rda")
dim(BIC_marker_expression)
dim(BIC_marker_list)

# deCS annotation
library(deCS)

#  Method A: Correlation analysis
data(BlueprintEncode_main)
data(Human_cell_landscape)
data(MonacoImmune_fine)
BIC_BlueprintEncode_main_cor = deCS.correlation(BIC_marker_expression, BlueprintEncode_main_t_score, methods = "pearson", top_n = 3)
BIC_HCL_cor = deCS.correlation(BIC_marker_expression, HCL_z_score, methods = "pearson", cell_type_threshold = 0.3, top_n = 3)
BIC_MonacoImmune_fine_cor = deCS.correlation(BIC_marker_expression, MonacoImmune_fine_t_score, methods = "pearson", cell_type_threshold = 0.3, top_n = 3)

BIC_anno_correlation_summary <- cbind(BIC_BlueprintEncode_main_cor[,c("Query", "Max_correlation", "Cell_labels")], BIC_HCL_cor[,c("Max_correlation", "Cell_labels")], BIC_MonacoImmune_fine_cor[,c("Max_correlation", "Cell_labels")])

#  Method B: Fisher's exact test
data(BlueprintEncode_main)
data(Human_cell_landscape)
data(MonacoImmune_fine)
BIC_BlueprintEncode_main_fisher = deCS.fisher(BIC_marker_list, BlueprintEncode_main_t_score, top_n = 3)
BIC_HCL_fisher = deCS.fisher(BIC_marker_list, HCL_z_score, cell_type_threshold = 0.3, top_n = 3)
BIC_MonacoImmune_fine_fisher = deCS.fisher(BIC_marker_list, MonacoImmune_fine_t_score, cell_type_threshold = 0.3, top_n = 3)

BIC_anno_fisher_summary <- cbind(BIC_BlueprintEncode_main_fisher[,c("Query", "Max_intersection_ratio", "Cell_labels")], BIC_HCL_fisher[,c("Max_intersection_ratio", "Cell_labels")], BIC_MonacoImmune_fine_fisher[,c("Max_intersection_ratio", "Cell_labels")])

# SingleR annotation
library(SingleR)
library(celldex)
ref1.BlueprintEncode <- celldex::BlueprintEncodeData()

BIC_singleR <- SingleR(test = BIC_marker_expression, ref = ref1.BlueprintEncode, assay.type.test = 1, labels = ref1.BlueprintEncode$label.main)
