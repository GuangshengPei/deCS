setwd("downloaded/deCS-master")

# Load trait associate genes (TAG) list
TAG = read.delim("Example_data/3.GWAS_trait_associate_genes/TAG.txt", head = T)

# deCS annotation
library(deCS)
data(BlueprintEncode_fine)
data(Human_cell_landscape)

pdf("TAG_BlueprintEncode.pdf", 15, 6)
TAG_BlueprintEncode_fine_enrich <- deCS.fisher(TAG, BlueprintEncode_fine_t_score, label = "count", intersection_threshold = 5, p_threshold = 0.05, cell_type_threshold = 5)
dev.off()

pdf("TAG_HumanCellLandscape.pdf", 18, 13)
TAG_HCL_enrich <- deCS.fisher(TAG, HCL_z_score, label = "count", intersection_threshold = 10, p_threshold = 0.01, cell_type_threshold = 10)
dev.off()
