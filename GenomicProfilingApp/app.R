
library(shiny)
library(bslib)

source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Define UI for application that draws a histogram
ui <- page_navbar(
  title = "Interactive Genomic Profiling",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Home",
            card(
              "Interactive genomic profiling tool to probe internal genomic features. Begin by navigating to the data processing page for file uploading before moving to the visualisations page to generate outputs"
            )),
  nav_panel(title = "Data Processing",
            layout_sidebar(
              sidebar = sidebar(
                open = TRUE,
                title = "File Upload",
                fileInput(inputId = "Annotation1",
                          label = "Upload annotation file (.bed/.gtf)",
                          accept = c(".bed", ".gtf")),
                fileInput(inputId = "Annotation2",
                          label = "Upload a second annotation file",
                          accept = c(".bed", ".gtf")),
                fileInput(inputId = "Sequence1",
                          label = "Upload sequence data files",
                          accept = c(".bw", ".bed"),
                          multiple = TRUE)
              ),
              card(
                full_screen = FALSE,
                card_header("Uploaded Annotation Files"),
                card_body(
                  textOutput("annofile1_name"),
                  textOutput("annofile2_name")
                )
              ),
              card(
                full_screen = FALSE,
                card_header("Uploaded Sequence Data file"),
                card_body(
                  textOutput("seqfile1_name")
                )
              )
            )
  ),
  nav_panel(title = "Visualisations",
            layout_sidebar(
              sidebar = sidebar(
                card(
                  full_screen = FALSE,
                  card_header = "Row Names",
                  checkboxInput(
                    inputId = "showrownames",
                    label = "Show row names",
                    value = FALSE
                  )
                )
              ),
              card(
                full_screen = TRUE,
                card_header("Heatmap Visualisation"),
                card_body(
                  plotOutput("enrichedHeatmapPlot")
                )
              )
            )
            )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$annofile1_name <- renderText({
    if(!is.null(input$Annotation1)) {
      paste("Annotation file 1:", input$Annotation1$name)
    } else {
      "No files uploaded"
    }
  })
  
  output$annofile2_name <- renderText({
    if(!is.null(input$Annotation2)) {
      paste("Annotation file 2:", input$Annotation2$name)
    } else {
      " "
    }
  })
  
  output$seqfile1_name <- renderText({
    if(!is.null(input$Sequence1)) {
      paste("Sequence Data File:", input$Sequence1$name)
    } else {
      "No sequence data file selected"
    }
  })
  
  show_row_names_reactive <- reactive({
    input$showrownames
  })
  
  output$enrichedHeatmapPlot <- renderPlot({
    hml <- hmList(matl = matl, wins = wins, col_fun = "red0", axis_labels = c("-20b", "TSS", "TES", "+20b"), show_row_names = show_row_names_reactive())
    hml[[1]] + hml[[2]]
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
