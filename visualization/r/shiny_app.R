# Shiny app for interactive single-cell RNA-Seq visualization (UMAP, Violin, Dot, Bar, Gene Expression)
library(shiny)
library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("../../scrna_seq/results/seurat_obj.rds")
features <- VariableFeatures(seurat_obj)
idents <- as.character(unique(Idents(seurat_obj)))

ui <- fluidPage(
  titlePanel("Single-cell RNA-Seq Interactive Dashboard"),
  sidebarLayout(
    sidebarPanel(
      selectInput("plot_type", "Plot type:",
                  choices = c("UMAP" = "umap", "Violin" = "violin", "Dot" = "dot", "Bar" = "bar", "Gene Expression" = "gene")),
      selectInput("feature", "Feature (gene) to plot:", choices = features),
      conditionalPanel(
        condition = "input.plot_type == 'bar'",
        selectInput("ident", "Group (cluster/ident):", choices = idents)
      )
    ),
    mainPanel(
      plotOutput("mainPlot")
    )
  )
)

server <- function(input, output) {
  output$mainPlot <- renderPlot({
    if (input$plot_type == "umap") {
      DimPlot(seurat_obj, reduction = "umap", label = TRUE) + theme_minimal()
    } else if (input$plot_type == "violin") {
      VlnPlot(seurat_obj, features = input$feature) + theme_minimal()
    } else if (input$plot_type == "dot") {
      DotPlot(seurat_obj, features = input$feature) + theme_minimal()
    } else if (input$plot_type == "bar") {
      # Bar plot: average expression of selected gene by group
      avg_expr <- AverageExpression(seurat_obj, features = input$feature, group.by = "ident")[[1]]
      df <- data.frame(Group = rownames(avg_expr), Expression = avg_expr[, input$feature])
      ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        labs(title = paste("Average Expression of", input$feature, "by Group"))
    } else if (input$plot_type == "gene") {
      FeaturePlot(seurat_obj, features = input$feature) + theme_minimal()
    }
  })
}

shinyApp(ui = ui, server = server)
