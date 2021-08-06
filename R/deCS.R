
deCS.correlation = function(markers_expression, ref_panel, methods = "pearson", ref_panel_markers_ratio = 0.05, top_n = 3, p.adjust.methods = "BH", cor_threshold = 0.1, p_threshold = 0.05, cell_type_threshold = 0, plot_figure = TRUE){
	markers_expression = na.omit(markers_expression)		
	cell_markers = cell_markers_full = c()
			for (i in 1:ncol(ref_panel)) {
				ii <- which(as.numeric(as.vector(ref_panel[, i])) > quantile(as.numeric(as.vector(ref_panel[, i])), probs = 1 - ref_panel_markers_ratio, na.rm = TRUE))
				cell_markers = rownames(ref_panel)[ii]
				cell_markers_full = c(cell_markers_full, cell_markers)
			}
	cell_markers_full = unique(cell_markers_full)
	common_gene = intersect(rownames(markers_expression), cell_markers_full)	
	query_ordered = markers_expression[match(common_gene, rownames(markers_expression)),]
	panel_ordered = ref_panel[match(common_gene, rownames(ref_panel)),]
	#Correlation analysis	
		correlation_cor = cor(query_ordered, panel_ordered, method = methods)
		correlation_cor[correlation_cor < 0] <- 0
		query_cluster <- as.data.frame(rownames(correlation_cor))
		colnames(query_cluster) <- "Query"
		query_cluster$Max_correlation <- apply(correlation_cor, 1, max)		
		query_match_order <- matrix(NA, nrow = nrow(correlation_cor), ncol = top_n)
			for (i in 1:ncol(query_match_order)){
				query_match_order[,i] <- colnames(correlation_cor)[apply(correlation_cor, 1, order)[ncol(correlation_cor)-i+1,]]
			}
		colnames(query_match_order)	= paste0("Top", 1:ncol(query_match_order))
		Cell_labels <- query_match_order[,1]
		Cell_labels[query_cluster$Max_correlation < cor_threshold] <- "Undetermined cells"
		correlation_annotation <- cbind(query_cluster, query_match_order, Cell_labels)
	#Figure plot checking	
		if (nrow(correlation_annotation) > 100){
			plot_figure = FALSE
			print ("The plot function for maximum cell clusters should less than 100.")
		}	
		if (plot_figure == TRUE) {
			correlation_pvalue = correlation_cor
			for (i in 1:nrow(correlation_cor)){
				for (j in 1:ncol(correlation_cor)){
					correlation_pvalue[i, j] = cor.test(query_ordered[, i], panel_ordered[, j], alternative = "greater", method = methods)$p.value
				}
				correlation_pvalue[i,] = p.adjust(correlation_pvalue[i,], p.adjust.methods)
			}	
			library(reshape2)
			require(ggplot2)
				correlation_cor_annotated <- t(correlation_cor[,which(apply(correlation_cor, 2, max) >= cell_type_threshold)])
				correlation_pvalue_annotated <- t(correlation_pvalue[,which(apply(correlation_cor, 2, max) >= cell_type_threshold)])		
				anno_list_cor = melt(correlation_cor_annotated)
				anno_list_pvalue = melt(correlation_pvalue_annotated)
				anno_list_combined = cbind(anno_list_cor, anno_list_pvalue)[,c(1,2,3,6)]
				colnames(anno_list_combined) = c("Var1", "Var2", "Correlation", "p_value")
				anno_list_combined$Label = round(anno_list_combined$Correlation, 2)
				anno_list_combined$Label[anno_list_combined$p_value > p_threshold] <- ""
				anno_list_combined$Label[anno_list_combined$Correlation < cor_threshold] <- ""
				anno_list_combined$Log10_Pvalue = -log(anno_list_combined$p_value, 10)
				anno_list_combined$Log10_Pvalue[anno_list_combined$Log10_Pvalue > 5] <- 5
				p <- ggplot(anno_list_combined, aes(y = (factor((Var1))), x = factor(Var2))) + geom_tile(aes(fill = Correlation)) + scale_fill_continuous(low = "white", high = "white") + guides(fill = FALSE)
				p <- p + geom_point(aes(size = Log10_Pvalue, colour = Correlation)) + scale_size(range = c(1, 5)) + theme_bw() + xlab("Cell cluster") + ylab("Cell type in reference panel") 
				p <- p + scale_colour_gradientn(colours = c("white", "green", "yellow", "orange", "red"), limits = range(anno_list_combined$Correlation)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
				p <- p + geom_text(aes(label = Label, size = 2.5)) + scale_y_discrete(limits = rev(as.character(unique(anno_list_combined$Var1)))) 
				print(p)
		}
	return(correlation_annotation)		
}


deCS.fisher = function(markers_list, ref_panel, type = "score", ref_panel_markers_ratio = 0.05, top_n = 3, p.adjust.methods = "BH", label = "ratio", intersection_threshold = 0.1, p_threshold = 0.05, cell_type_threshold = 0, plot_figure = TRUE){
		if (length(intersect(c("cluster", "gene"), colnames(markers_list))) < 2) {
			print("Please input a gene list with at least 2 columns, names with cluster and gene!")
			stop
		}
		if (label == "ratio"){
				if (type == "score"){
					intersection_ratio = matrix(0, nrow = length(unique(markers_list$cluster)), ncol = ncol(ref_panel))
					FET_pvalue = matrix(1, nrow = length(unique(markers_list$cluster)), ncol = ncol(ref_panel))
					rownames(intersection_ratio) = rownames(FET_pvalue) = unique(markers_list$cluster)
					colnames(intersection_ratio) = colnames(FET_pvalue) = colnames(ref_panel)
						for (i in 1:length(unique(markers_list$cluster))){
							for (j in 1:ncol(ref_panel)) {
								query_markers = markers_list[which(markers_list$cluster == as.character(unique(markers_list$cluster)[i])),]$gene 
								ii <- which(as.numeric(as.vector(ref_panel[, j])) > quantile(as.numeric(as.vector(ref_panel[, j])), probs = 1 - ref_panel_markers_ratio, na.rm = TRUE))
								ref_markers = rownames(ref_panel)[ii]
								common = length(intersect(query_markers, ref_markers))
								query_unique = length(setdiff(query_markers, ref_markers))
								ref_unique = length(setdiff(ref_markers, query_markers))		
								remain = length(unique(c(unique(markers_list$gene), rownames(ref_panel)))) - common - query_unique - ref_unique
								overlapped_proportion = common/length(query_markers)
								intersection_ratio[i, j] = overlapped_proportion
								FET_pvalue[i, j] = fisher.test(matrix(c(common, query_unique, ref_unique, remain), ncol = 2))$p.value	
							}
						}
						for (i in 1:nrow(FET_pvalue)){
							FET_pvalue[i,] = p.adjust(FET_pvalue[i,], p.adjust.methods)
						}	
					query_cluster <- as.data.frame(rownames(intersection_ratio))
					colnames(query_cluster) <- "Query"
					query_cluster$Max_intersection_ratio <- apply(intersection_ratio, 1, max)		
					query_match_order <- matrix(NA, nrow = nrow(intersection_ratio), ncol = top_n)
						for (i in 1:ncol(query_match_order)){
							query_match_order[,i] <- colnames(intersection_ratio)[apply(intersection_ratio, 1, order)[ncol(intersection_ratio)-i+1,]]
						}	
					colnames(query_match_order)	= paste0("Top", 1:ncol(query_match_order))
					Cell_labels <- query_match_order[,1]
					Cell_labels[query_cluster$Max_intersection_ratio < intersection_threshold] <- "Undetermined cells"
					FET_annotation <- cbind(query_cluster, query_match_order, Cell_labels)		
				}
				if (type == "list"){
					if (length(intersect(c("cluster", "gene"), colnames(markers_list))) < 2) {
						print("Please input a reference marker gene list with at least 2 columns, names with Cell_type and Marker_gene!")
						stop
					}
					intersection_ratio = matrix(0, nrow = length(unique(markers_list$cluster)), ncol = length(unique(ref_panel$Cell_type)))
					FET_pvalue = matrix(1, nrow = length(unique(markers_list$cluster)), ncol = length(unique(ref_panel$Cell_type)))
					rownames(intersection_ratio) = rownames(FET_pvalue) = unique(markers_list$cluster)
					colnames(intersection_ratio) = colnames(FET_pvalue) = unique(ref_panel$Cell_type)	
						for (i in 1:length(table(markers_list$cluster))) {
							for (j in 1:length(unique(ref_panel$Cell_type))) {
								query_markers = markers_list[which(markers_list$cluster == as.character(unique(markers_list$cluster)[i])),]$gene
								ref_markers <- as.character(ref_panel$Marker_gene)[ref_panel$Cell_type == as.character(unique(ref_panel$Cell_type)[j])]
								common = length(intersect(query_markers, ref_markers))
								query_unique = length(setdiff(query_markers, ref_markers))
								ref_unique = length(setdiff(ref_markers, query_markers))		
								remain = length(unique(c(unique(markers_list$gene), rownames(ref_panel)))) - common - query_unique - ref_unique
								overlapped_proportion = common/length(query_markers)
								intersection_ratio[i, j] = overlapped_proportion
								FET_pvalue[i, j] = fisher.test(matrix(c(common, query_unique, ref_unique, remain), ncol = 2))$p.value	
							}
						}
						for (i in 1:nrow(FET_pvalue)){
							FET_pvalue[i,] = p.adjust(FET_pvalue[i,], p.adjust.methods)
						}	
					query_cluster <- as.data.frame(rownames(intersection_ratio))
					colnames(query_cluster) <- "Query"
					query_cluster$Max_intersection_ratio <- apply(intersection_ratio, 1, max)		
					query_match_order <- matrix(NA, nrow = nrow(intersection_ratio), ncol = top_n)
						for (i in 1:ncol(query_match_order)){
							query_match_order[,i] <- colnames(intersection_ratio)[apply(intersection_ratio, 1, order)[ncol(intersection_ratio)-i+1,]]
						}
					colnames(query_match_order)	= paste0("Top", 1:ncol(query_match_order))
					Cell_labels <- query_match_order[,1]
					Cell_labels[query_cluster$Max_intersection_ratio < intersection_threshold] <- "Undetermined cells"
					FET_annotation <- cbind(query_cluster, query_match_order, Cell_labels)	
				}
			if (nrow(FET_annotation) > 100){
				plot_figure = FALSE
				print ("The plot function for maximum cell clusters should less than 100.")
			}				
				library(reshape2)
				require(ggplot2)
					intersection_ratio_annotated <- t(intersection_ratio[,which(apply(intersection_ratio, 2, max) >= cell_type_threshold)])
					FET_pvalue_annotated <- t(FET_pvalue[,which(apply(intersection_ratio, 2, max) >= cell_type_threshold)])		
					anno_list_intersection_ratio = melt(intersection_ratio_annotated)
					anno_list_FET_pvalue = melt(FET_pvalue_annotated)
					anno_list_combined = cbind(anno_list_intersection_ratio, anno_list_FET_pvalue)[,c(1,2,3,6)]	
					colnames(anno_list_combined) = c("Var1", "Var2", "Intersection_ratio", "p_value")
					anno_list_combined$Label = round(anno_list_combined$Intersection_ratio, 2)
					anno_list_combined$Label[anno_list_combined$p_value > p_threshold] <- ""
					anno_list_combined$Label[anno_list_combined$Intersection_ratio < intersection_threshold] <- ""
					anno_list_combined$Log10_Pvalue = -log(anno_list_combined$p_value, 10)
					anno_list_combined$Log10_Pvalue[anno_list_combined$Log10_Pvalue > 5] <- 5
					p <- ggplot(anno_list_combined, aes(y = (factor((Var1))), x = factor(Var2))) + geom_tile(aes(fill = Intersection_ratio)) + scale_fill_continuous(low = "white", high = "white") + guides(fill = FALSE)
					p <- p + geom_point(aes(size = Log10_Pvalue, colour = Intersection_ratio)) + scale_size(range = c(1, 5)) + theme_bw() + xlab("Cell cluster") + ylab("Cell type in reference panel") 
					p <- p + scale_colour_gradientn(colours = c("white", "green", "yellow", "orange", "red"), limits = range(anno_list_combined$Intersection_ratio)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))	
					p <- p + geom_text(aes(label = Label, size = 2.5)) + scale_y_discrete(limits = rev(as.character(unique(anno_list_combined$Var1)))) 
					print(p)
			}
		
		if (label == "count"){
				if (type == "score"){
					intersection_count = matrix(0, nrow = length(unique(markers_list$cluster)), ncol = ncol(ref_panel))
					FET_pvalue = matrix(1, nrow = length(unique(markers_list$cluster)), ncol = ncol(ref_panel))
					rownames(intersection_count) = rownames(FET_pvalue) = unique(markers_list$cluster)
					colnames(intersection_count) = colnames(FET_pvalue) = colnames(ref_panel)
						for (i in 1:length(unique(markers_list$cluster))){
							for (j in 1:ncol(ref_panel)) {
								query_markers = markers_list[which(markers_list$cluster == as.character(unique(markers_list$cluster)[i])),]$gene 
								ii <- which(as.numeric(as.vector(ref_panel[, j])) > quantile(as.numeric(as.vector(ref_panel[, j])), probs = 1 - ref_panel_markers_ratio, na.rm = TRUE))
								ref_markers = rownames(ref_panel)[ii]
								common = length(intersect(query_markers, ref_markers))
								query_unique = length(setdiff(query_markers, ref_markers))
								ref_unique = length(setdiff(ref_markers, query_markers))		
								remain = length(unique(c(unique(markers_list$gene), rownames(ref_panel)))) - common - query_unique - ref_unique
								intersection_count[i, j] = common
								FET_pvalue[i, j] = fisher.test(matrix(c(common, query_unique, ref_unique, remain), ncol = 2))$p.value	
							}
						}
						for (i in 1:nrow(FET_pvalue)){
							FET_pvalue[i,] = p.adjust(FET_pvalue[i,], p.adjust.methods)
						}	
					query_cluster <- as.data.frame(rownames(intersection_count))
					colnames(query_cluster) <- "Query"
					query_cluster$Max_intersection_count <- apply(intersection_count, 1, max)		
					query_match_order <- matrix(NA, nrow = nrow(intersection_count), ncol = top_n)
						for (i in 1:ncol(query_match_order)){
							query_match_order[,i] <- colnames(intersection_count)[apply(intersection_count, 1, order)[ncol(intersection_count)-i+1,]]
						}
					colnames(query_match_order)	= paste0("Top", 1:ncol(query_match_order))
					Cell_labels <- query_match_order[,1]
					Cell_labels[query_cluster$Max_intersection_count < intersection_threshold] <- "Undetermined cells"
					FET_annotation <- cbind(query_cluster, query_match_order, Cell_labels)		
				}
				if (type == "list"){
					if (length(intersect(c("cluster", "gene"), colnames(markers_list))) < 2) {
						print("Please input a reference marker gene list with at least 2 columns, names with Cell_type and Marker_gene!")
						stop
					}
					FET_result = list()
					intersection_count = matrix(0, nrow = length(unique(markers_list$cluster)), ncol = length(unique(ref_panel$Cell_type)))
					FET_pvalue = matrix(1, nrow = length(unique(markers_list$cluster)), ncol = length(unique(ref_panel$Cell_type)))
					rownames(intersection_count) = rownames(FET_pvalue) = unique(markers_list$cluster)
					colnames(intersection_count) = colnames(FET_pvalue) = unique(ref_panel$Cell_type)	
						for (i in 1:length(table(markers_list$cluster))) {
							for (j in 1:length(unique(ref_panel$Cell_type))) {
								query_markers = markers_list[which(markers_list$cluster == as.character(unique(markers_list$cluster)[i])),]$gene
								ref_markers <- as.character(ref_panel$Marker_gene)[ref_panel$Cell_type == as.character(unique(ref_panel$Cell_type)[j])]
								common = length(intersect(query_markers, ref_markers))
								query_unique = length(setdiff(query_markers, ref_markers))
								ref_unique = length(setdiff(ref_markers, query_markers))		
								remain = length(unique(c(unique(markers_list$gene), rownames(ref_panel)))) - common - query_unique - ref_unique
								intersection_count[i, j] = common
								FET_pvalue[i, j] = fisher.test(matrix(c(common, query_unique, ref_unique, remain), ncol = 2))$p.value	
							}
						}
						for (i in 1:nrow(FET_pvalue)){
							FET_pvalue[i,] = p.adjust(FET_pvalue[i,], p.adjust.methods)
						}	
					query_cluster <- as.data.frame(rownames(intersection_count))
					colnames(query_cluster) <- "Query"
					query_cluster$Max_intersection_count <- apply(intersection_count, 1, max)		
					query_match_order <- matrix(NA, nrow = nrow(intersection_count), ncol = top_n)
						for (i in 1:ncol(query_match_order)){
							query_match_order[,i] <- colnames(intersection_count)[apply(intersection_count, 1, order)[ncol(intersection_count)-i+1,]]
						}
					colnames(query_match_order)	= paste0("Top", 1:ncol(query_match_order))
					Cell_labels <- query_match_order[,1]
					Cell_labels[query_cluster$Max_intersection_count < intersection_threshold] <- "Undetermined cells"
					FET_annotation <- cbind(query_cluster, query_match_order, Cell_labels)					
				}
			#Figure plot checking			
			if (nrow(FET_annotation) > 100){
				plot_figure = FALSE
				print ("The plot function for maximum cell clusters should less than 100.")
			}				
				library(reshape2)
				require(ggplot2)
					intersection_count_annotated <- t(intersection_count[,which(apply(intersection_count, 2, max) >= cell_type_threshold)])
					FET_pvalue_annotated <- t(FET_pvalue[,which(apply(intersection_count, 2, max) >= cell_type_threshold)])		
					anno_list_intersection_count = melt(intersection_count_annotated)
					anno_list_FET_pvalue = melt(FET_pvalue_annotated)
					anno_list_combined = cbind(anno_list_intersection_count, anno_list_FET_pvalue)[,c(1,2,3,6)]	
					colnames(anno_list_combined) = c("Var1", "Var2", "Intersection_count", "p_value")
					anno_list_combined$Label = round(anno_list_combined$Intersection_count, 2)
					anno_list_combined$Label[anno_list_combined$p_value > p_threshold] <- ""
					anno_list_combined$Label[anno_list_combined$Intersection_count < intersection_threshold] <- ""
					anno_list_combined$Log10_Pvalue = -log(anno_list_combined$p_value, 10)
					anno_list_combined$Log10_Pvalue[anno_list_combined$Log10_Pvalue > 5] <- 5
					p <- ggplot(anno_list_combined, aes(y = (factor((Var1))), x = factor(Var2))) + geom_tile(aes(fill = Intersection_count)) + scale_fill_continuous(low = "white", high = "white") + guides(fill = FALSE)
					p <- p + geom_point(aes(size = Log10_Pvalue, colour = Intersection_count)) + scale_size(range = c(1, 5)) + theme_bw() + xlab("Cell cluster") + ylab("Cell type in reference panel") 
					p <- p + scale_colour_gradientn(colours = c("white", "green", "yellow", "orange", "red"), limits = range(anno_list_combined$Intersection_count)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))	
					p <- p + geom_text(aes(label = Label, size = 2.5)) + scale_y_discrete(limits = rev(as.character(unique(anno_list_combined$Var1)))) 
					print(p)
				}
	return(FET_annotation)
}
