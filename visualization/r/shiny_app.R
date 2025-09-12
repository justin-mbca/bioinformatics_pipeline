# Basic Shiny app for interactive UMAP plot (single-cell RNA-Seq)
library(shiny)
library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("../../scrna_seq/results/seurat_obj.rds")

ui <- fluidPage(
  titlePanel("Single-cell UMAP Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("feature", "Feature to plot:", choices = VariableFeatures(seurat_obj))
    ),
    mainPanel(
      plotOutput("umapPlot"),
      plotOutput("featurePlot")
    )
  )
)

server <- function(input, output) {
  output$umapPlot <- renderPlot({
    DimPlot(seurat_obj, reduction = "umap", label = TRUE) + theme_minimal()
  })
  output$featurePlot <- renderPlot({
    FeaturePlot(seurat_obj, features = input$feature) + theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
