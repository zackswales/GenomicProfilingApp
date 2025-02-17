
library(shiny)
library(bslib)

source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Define UI for application that draws a histogram
ui <- page_navbar(
  title = "Interactive Genomic Profiling",
  bg = "#E46303",
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
                fileInput(inputId = "Region1",
                          label = "Upload a region file (.bed/.gtf)",
                          accept = c(".bed", ".gtf")),
                fileInput(inputId = "Sequence1",
                          label = "Upload sequence data files (.bw)",
                          accept = ".bw",
                          multiple = TRUE),
                textInput(inputId = "firstmatrixname",
                          label = "Matrix Name"),
                actionButton(inputId = "matrixgeneration",
                             label = "Generate Matrices"),
                helpText("Warning: matrix generation may take a few seconds")
              ),
                card(
                full_screen = FALSE,
                card_header("Uploaded Region Files"),
                card_body(
                  textOutput("annofile1_name"))
              ),
              card(
                full_screen = FALSE,
                card_header("Uploaded Sequence Data files"),
                card_body(
                  textOutput("seqfile1_name")
                )
              ),
              card(
                full_screen = FALSE,
                card_header("Generated matrices"),
                card_body(
                  textOutput("matrixnames")
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
    if(!is.null(input$Region1)) {
      paste("Annotation file 1:", input$Region1$name)
    } else {
      "No files uploaded"
    }
  })
  
  output$seqfile1_name <- renderText({
    if(!is.null(input$Sequence1)) {
      paste("Sequence Data File:", input$Sequence1$name)
    } else {
      "No sequence data file selected"
    }
  })
  
  # Matrix list generation
  
  matl <- eventReactive(input$matrixgeneration, {
    req(input$Region1, input$Sequence1, input$firstmatrixname)
    region_file <- input$Region1$datapath
    bigwig_files <- input$Sequence1$datapath
    print(paste("BED/GTF File Path:", region_file))
    print(paste("BigWig File Path(s):", paste(bigwig_files, collapse = ", ")))
    b <- import(region_file)
    features <- getFeature(b)
    fbw <- bigwig_files
    names(fbw) <- basename(fbw)
    bwf <- importBWlist(fbw, names(fbw), selection = features)
    grl <- list("features" = features)
    matl <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = 10, w = 1, strand = "no")
    if(length(matl) > 0) {
      names(matl)[1] <- input$firstmatrixname
    }
    print("Generated matrices:")
    print(names(matl))
    print(matl[1])
    return(matl)
  })
  
  
  output$matrixnames <- renderText({
    req(matl())
    mat_names <- names(matl())
    if(is.null(mat_names) || length(mat_names) == 0) {
      return("No matrices generated")
    }
    paste(mat_names, collapse = ", ")
  })
  
  show_row_names_reactive <- reactive({
    input$showrownames
  })
  
  wins = c("Upstream" = 10, "Feature" = 20, "Downstream" = 10)
  
  output$enrichedHeatmapPlot <- renderPlot({
    hml <- hmList(matl = matl(), wins = wins, col_fun = "red0", axis_labels = c("-10b", "TSS", "TES", "+10b"), show_row_names = show_row_names_reactive())
    hml[[1]]
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
