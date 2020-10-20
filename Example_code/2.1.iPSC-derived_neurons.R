setwd("downloaded/deCS-master")

# Load the bulk RNA-seq datasets and data preprocess
data = read.delim("Example_data/2.1.iPSC-derived_neurons/Neuron_TPM.txt", head = T)
tag_temp = names(table(data$Symbol)[table(data$Symbol) > 1])
tag = as.character(data$Symbol[!(data$Symbol %in% tag_temp)])
data = data[which(data$Symbol %in% tag),]
rownames(data) = data$Symbol
data = data[, -1]
data_clean <- data[which(apply(data, 1, max) > 1),]

# Bulk RNA-seq normalization
library(deTS)
data(correction_factor)
Neurons_normalized = tsea.expression.normalization(data_clean, correction_factor, normalization = "abundance")

# deCS annotation
library(deCS)
data(Human_cell_landscape)
Neurons_HCL_anno_cor = deCS.correlation(Neurons_normalized, HCL_z_score, ref_panel_markers_ratio = 0.01, p.adjust.methods = "none", cor_threshold = 0.1, cell_type_threshold = 0.2)


