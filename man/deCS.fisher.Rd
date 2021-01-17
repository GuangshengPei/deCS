\name{deCS.fisher}
\alias{deCS.fisher}
\title{Cell type enrichment analysis by Fisher's exact test
}
\description{Cell type enrichment analysis by Fisher's exact test for user provide human gene list.
}
\usage{
deCS.fisher(markers_list, ref_panel, type = "score", ref_panel_markers_ratio = 0.05, top_n = 3, p.adjust.methods = "BH", label = "ratio", 
intersection_threshold = 0.1, p_threshold = 0.05, cell_type_threshold = 0, plot_figure = TRUE)
}
\arguments{
  \item{markers_list}{Gene list with at least two columns with name "cluster" (cell type) and "gene" (marker gene list).
}
  \item{ref_panel}{Gene cell type specific t-score or z-score matrix, can be loaded by data(ref_panel), ref_panel = c("BlueprintEncode_main", "BlueprintEncode_fine", "Human_cell_landscape", "Human_cell_atlas_of_fetal", "DICE_main", "DICE_fine", "MonacoImmune_main", "MonacoImmune_fine", "CellMatch").
}
  \item{type}{When reference is t-score or z-score matrix, specifify type = "score", when reference is cell type-marker gene list, specifify type = "list".
}
  \item{ref_panel_markers_ratio}{The threshold to define cell type-specific genes (with top t-score or z-score), the default value is 0.05.
}
  \item{top_n}{deCS result will list top n most similar cell types, the default value is 3.
}
  \item{p.adjust.methods}{p.adjust.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
}
  \item{label}{Calculate the gene intersection ratio between query and reference CSGs genes, specifify label = "ratio"; or calculate the gene intersection count, specifify label = "count".
}
  \item{intersection_threshold}{The minimum correlation intersection threshold for cell type annotation, the default value is 0.1. When label = "count", please provide a larger number, e.g. intersection_threshold = 10. Cells failed to pass the threshold will be regarded as "Undetermined cells".
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
\author{
Guangsheng Pei
}
\note{
}
\seealso{https://github.com/GuangshengPei/deCS
}
\examples{
deCS.fisher
}
\keyword{ "deCS" }
\keyword{ "Fisher's exact test" }
