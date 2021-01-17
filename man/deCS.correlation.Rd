\name{deCS.correlation}
\alias{deCS.correlation}
\title{Cell type enrichment analysis by correlation analysis
}
\description{Cell type enrichment analysis by Pearson correlation coefficient of Spearman's rank correlation coefficient for human scRNA-seq or bulk RNA-seq expression profile.
}
\usage{
deCS.correlation(markers_expression, ref_panel, methods = "pearson", ref_panel_markers_ratio = 0.05, top_n = 3, p.adjust.methods = "BH", 
cor_threshold = 0.1, p_threshold = 0.05, cell_type_threshold = 0, plot_figure = TRUE)
}
\arguments{
  \item{markers_expression}{Gene scaled expression matrix for human scRNA-seq or bulk RNA-seq.
}
  \item{ref_panel}{Gene cell type specific t-score or z-score matrix, can be loaded by data(ref_panel), ref_panel = c("BlueprintEncode_main", "BlueprintEncode_fine", "Human_cell_landscape", "Human_cell_atlas_of_fetal", "DICE_main", "DICE_fine", "MonacoImmune_main", "MonacoImmune_fine", "CellMatch").
}
  \item{methods}{methods = c("pearson", "spearman").
}
  \item{ref_panel_markers_ratio}{The threshold to define cell type specific genes (with top t-score or z-score), the default value is 0.05.
}
  \item{top_n}{deCS result will list top n most similar cell types, the default value is 3.
}
  \item{p.adjust.methods}{p.adjust.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
}
  \item{cor_threshold}{The minimum correlation coefficient threshold for cell type annotation, the default value is 0.1. Cells failed to pass the threshold will be regarded as "Undetermined cells".
}
  \item{p_threshold}{The minimum p-value for cell type annotation, the default value is 0.05. Cells failed to pass the threshold will be regarded as "Undetermined cells".
}
  \item{cell_type_threshold}{Figure plot only, only show significant correlated cell types, the default value is 0.
}
  \item{plot_figure}{Plot figure or not, the maximum query cell number is 100, the default value is TRUE.
}
}
\details{
}
\value{
}
\references{
}
\author{Guangsheng Pei	
}
\note{
}
\seealso{https://github.com/GuangshengPei/deCS
}
\examples{
deCS.correlation
}
\keyword{ "deCS" }
\keyword{ "Correlation" }
