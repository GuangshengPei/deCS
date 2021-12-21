options(shiny.maxRequestSize = 300*1024^2) 
# Load R libs ####
library(deCS)
library(dplyr)
library(Seurat)

run_deCS <- function(input, min.pct = 0.25, top_n = 10, logfc.threshold = 0.25){
        Idents(input) <- input$seurat_clusters
        input.markers <- FindAllMarkers(input, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
        topN <- input.markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC)
        input_cluster_average = log(AverageExpression(input)[[1]] + 1, 2)
        input_cluster_marker_average = input_cluster_average[which(rownames(input_cluster_average) %in% topN$gene), ]
        input_cluster_marker_z_score = t(scale(t(input_cluster_marker_average)))
        input_cluster_marker_z_score
}

plot_main <- function(input, p_threshold = 2, cell_type_threshold = 0.2, 
                      methods = "pearson", ref_panel_markers_ratio = 0.05,
                      reference_panel = 'MonacoImmune_main'){
        if (reference_panel=="MonacoImmune_main"){
                data(MonacoImmune_main)
                deCS.correlation(input, MonacoImmune_main_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="MonacoImmune_fine"){
                data(MonacoImmune_fine)
                deCS.correlation(input, MonacoImmune_fine_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="BlueprintEncode_main"){
                data(BlueprintEncode_main)
                deCS.correlation(input, BlueprintEncode_main_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="BlueprintEncode_fine"){
                data(BlueprintEncode_fine)
                deCS.correlation(input, BlueprintEncode_fine_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="Database of Immune Cell Expression (DICE)_main"){
                data(DICE_main)
                deCS.correlation(input, DICE_main_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="Database of Immune Cell Expression (DICE)_fine"){
                data(DICE_fine)
                deCS.correlation(input, DICE_fine_t_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="Human cell landscape"){
                data(Human_cell_landscape)
                deCS.correlation(input, HCL_z_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } else if (reference_panel=="Human cell atlas of fetal"){
                data(Human_cell_atlas_of_fetal)
                deCS.correlation(input, HCAF_z_score, 
                                 p_threshold = -log10(p_threshold), cell_type_threshold = cell_type_threshold,
                                 methods = methods, ref_panel_markers_ratio = ref_panel_markers_ratio)
        } 
        
}

# Define ui ####
ui <- fluidPage(
        headerPanel("deCS: A tool for systematic cell type annotations of single-cell RNA sequencing data among human tissues"),
        sidebarLayout(
                sidebarPanel(
                        fileInput("file", "Choose data File",
                                  multiple = F,
                                  accept = c("rds")),
                        selectInput("reference_panel", "Select reference:", 
                                    choices = c('BlueprintEncode_main','BlueprintEncode_fine',
                                                'Database of Immune Cell Expression (DICE)_main',
                                                'Database of Immune Cell Expression (DICE)_fine',
                                                'MonacoImmune_main','MonacoImmune_fine',
                                                'Human cell landscape', 'Human cell atlas of fetal'), 
                                    selected = "BlueprintEncode_main"),
                        selectInput("methods", "Select method:", 
                                    choices = c('pearson','spearman'), 
                                    selected = "pearson"),
                        sliderInput("top_n", "top_n",
                                    min = 5, max = 100, value = 10),
                        sliderInput("min.pct", "min.pct",
                                    min = 0, max = 1, value = 0.1, step = 0.001),
                        sliderInput("logfc.threshold", "logfc.threshold",
                                    min = 0.01, max = 1, value = 0.25, step = 0.001),
                        sliderInput("p_threshold", "-log(p_threshold)",
                                    min = 0.01, max = 5, value = 1, step = 0.1),
                        sliderInput("cell_type_threshold", "cell_type_threshold",
                                    min = 0, max = 1, value = 0.5, step = 0.001),
                        sliderInput("ref_panel_markers_ratio", "ref_panel_markers_ratio",
                                    min = 0, max = 1, value = 0.05, step = 0.001)
                ),
                mainPanel(
                tabsetPanel(
                        tabPanel("deCS output",plotOutput('plot1')))
                )
        )
)

# Define the server code
server <- function(input, output) {
        data_file <- reactive({
                req(input$file)
                readRDS(input$file$datapath)
        })
        
        process_file <- reactive({
                run_deCS(input = data_file(), min.pct = input$min.pct, top_n = input$top_n)
        })
        
        output$plot1 <- renderPlot({
                plot_main(input = process_file(), p_threshold = input$p_threshold, 
                          cell_type_threshold = input$cell_type_threshold,
                          reference_panel = input$reference_panel,
                          methods = input$methods, ref_panel_markers_ratio=input$ref_panel_markers_ratio)
        })
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)