setwd("downloaded/deCS-master")

# Load the preprocessed cardiac cells dataset
load("Example_data/1.3.Cardiac_cells/Cardiac_cells.rda")
dim(Cardiac_cells_expression)

# deCS annotation
library(deCS)
data(BlueprintEncode_fine)
Cardiac_BlueprintEncode_fine_pcc = deCS.correlation(Cardiac_cells_expression, BlueprintEncode_fine_t_score, cor_threshold = 0.1, methods = "pearson")
Cardiac_BlueprintEncode_fine_spearman = deCS.correlation(Cardiac_cells_expression, BlueprintEncode_fine_t_score, cor_threshold = 0.1, methods = "spearman")

data(Human_cell_landscape)
Cardiac_HCL_cor_pcc = deCS.correlation(Cardiac_cells_expression, HCL_z_score, methods = "pearson", cell_type_threshold = 0.3)
Cardiac_HCL_cor_spearman = deCS.correlation(Cardiac_cells_expression, HCL_z_score, methods = "spearman", cell_type_threshold = 0.3)

# SingleR annotation
library(SingleR)
library(celldex)
ref1.BlueprintEncode <- celldex::BlueprintEncodeData()

Cardiac_singleR <- SingleR(test = Cardiac_cells_expression, ref = ref1.BlueprintEncode, assay.type.test = 1, labels = ref1.BlueprintEncode$label.main)