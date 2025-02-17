
library(shiny)
library(bslib)

source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Define UI for application with three panels
ui <- page_navbar(
  title = "Interactive Genomic Profiling",
  bg = "#E46303",
  inverse = TRUE,
  # Panel for home page with directions for use of application
  nav_panel(title = "Home",
            card(
              "Interactive genomic profiling tool to probe internal genomic features. Begin by navigating to the data processing page for file uploading before moving to the visualisations page to generate outputs"
            )),
  # Panel for data processing and matrix formation
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
                          label = "First Matrix Name"),
                actionButton(inputId = "matrixgeneration",
                             label = "Generate Matrices"),
                helpText("Warning: only click once as matrix generation takes a few seconds to complete")
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
  # Panel for visualisations
  nav_panel(title = "Visualisations",
            layout_sidebar(
              sidebar = sidebar(
                open = TRUE, 
                title = "Visualisation Customisation",
                radioButtons(
                  inputId = "col_fun",
                  label = "Select colour scheme:",
                  choices = list("White to red" = 1, "Blue to red" = 2, "Red scale" = 3),
                  selected = 1
                ),
                sliderInput(
                  inputId = "quantiles",
                  label = "Set quantile range:",
                  min = 0,
                  max = 0.99,
                  value = c(0,0.99)
                ),
                sliderInput(
                  inputId = "ylim",
                  label = "Set y axis limit",
                  min = 0,
                  max = 20000,
                  value = 10000
                ),
                checkboxInput(
                  inputId = "showrownames",
                  label = "Show row names",
                  value = FALSE
                ),
                downloadButton(
                  outputId = "downloadpng",
                  label = "Download as .png"),
                downloadButton(
                  outputId = "downloadpdf",
                  label = "Download as .pdf"
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

# Define server logic required to compute matrices and plot outputs
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
      paste(input$Sequence1$name)
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
  
  # Creating reactive objects for the visualisation customisation options
  
  show_row_names_reactive <- reactive({
    input$showrownames
  })
  
  col_fun_reactive <- reactive({
    switch(input$col_fun,
           "1" = "red",
           "2" = "bl2rd",
           "3" = "red0")
  })
  
  min_quantile_reactive <- reactive({
    input$quantiles[1]
  })
  
  max_quantile_reactive <- reactive({
    input$quantiles[2]
  })
  
  ylim_reactive <- reactive({
    input$ylim
  })
  
  # Creating the heatmaps
  
  wins = c("Upstream" = 10, "Feature" = 20, "Downstream" = 10)
  
  output$enrichedHeatmapPlot <- renderPlot({
    req(matl())
    hml <- hmList(matl = matl(), wins = wins, col_fun = col_fun_reactive(), axis_labels = c("-10b", "TSS", "TES", "+10b"), show_row_names = show_row_names_reactive(), min_quantile = min_quantile_reactive(), max_quantile = max_quantile_reactive(), ylim = c(0, ylim_reactive()))
    req(length(hml) > 0)
    combined_hm <- Reduce(`+`, hml)
    draw(combined_hm, merge_legend = TRUE)
  })
  
  output$downloadpng <- downloadHandler(
    filename = function() { paste("heatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1200, height = 800, res = 150)  # Adjust quality
      hml <- hmList(
        matl = matl(), 
        wins = wins, 
        col_fun = col_fun_reactive(), 
        axis_labels = c("-10b", "TSS", "TES", "+10b"), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = min_quantile_reactive(), 
        max_quantile = max_quantile_reactive(), 
        ylim = c(0, ylim_reactive())
      )
      
      req(length(hml) > 0)
      combined_hm <- Reduce(`+`, hml)
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  # Close the PNG file
    }
  )
  
  output$downloadpdf <- downloadHandler(
    filename = function() { paste("heatmap_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 12, height = 8)  # Adjust quality
      hml <- hmList(
        matl = matl(), 
        wins = wins, 
        col_fun = col_fun_reactive(), 
        axis_labels = c("-10b", "TSS", "TES", "+10b"), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = min_quantile_reactive(), 
        max_quantile = max_quantile_reactive(), 
        ylim = c(0, ylim_reactive())
      )
      
      req(length(hml) > 0)
      combined_hm <- Reduce(`+`, hml)
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  # Close the pdf file
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
