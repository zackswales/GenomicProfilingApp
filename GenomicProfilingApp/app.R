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
                width = 325,
                title = "File Upload",
                fileInput(inputId = "Region1",
                          label = "Upload a region file (.bed/.gtf)",
                          accept = c(".bed", ".gtf"),
                          multiple = FALSE),
                fileInput(inputId = "Sequence1",
                          label = "Upload sequence data files (.bw)",
                          accept = ".bw",
                          multiple = TRUE),
                radioButtons(inputId = "strand",
                             label = "Select strandedness of data:",
                             choices = list("Unstranded" = 1, "Forward" = 2, "Reverse" = 3),
                             selected = 1),
                helpText("If all files have the same strandedness - select unstranded"),
                numericInput(inputId = "flank",
                             label = "Specify flank around feature:",
                             value = 20,
                             step = 10),
                sliderInput(inputId = "windowsize",
                            label = "Select window size for signal aggregation:",
                            min = 1,
                            max = 20,
                            value = 1),
                helpText("Warning: flank value must be divisible by window size"),
                actionButton(inputId = "matrixgeneration",
                             label = "Generate Matrices"),
                helpText("Warning: only click once as matrix generation takes a few seconds to complete")
              ),
              card(
                full_screen = FALSE,
                card_header("Uploaded Region Files"),
                card_body(
                  textOutput("regionfile1_name"))
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
                numericInput(
                  inputId = "maxylim",
                  label = "Upper y-axis limit for metaplot",
                  value = 5000,
                  step = 50
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
  output$regionfile1_name <- renderText({
    if(!is.null(input$Region1)) {
      paste(input$Region1$name)
    } else {
      "No files uploaded"
    }
  })
  
  output$seqfile1_name <- renderText({
    if(!is.null(input$Sequence1)) {
      paste(input$Sequence1$name, collapse = ", ")
    } else {
      "No files uploaded"
    }
  })
  
  # Creating reactive objects for the visualisation customisation options
  
  strand_reactive <- reactive({
    switch(input$strand,
           "1" = "no",
           "2" = "for",
           "3" = "rev")
  })
  
  flank_reactive <- reactive({
    input$flank
  })
  
  windowsize_reactive <- reactive({
    input$windowsize
  })
  
  # Matrix list generation
  
  matl <- eventReactive(input$matrixgeneration, {
    req(input$Region1, input$Sequence1)
    
    region_file <- input$Region1$datapath
    bigwig_files <- input$Sequence1$datapath
    bigwig_file_names <- input$Sequence1$name
    
    print(paste("BED/GTF File Path:", region_file))
    print(paste("BigWig File Path(s):", paste(bigwig_files, collapse = ", ")))

    # Import region files and get features
    
    b <- import(region_file)
    features <- getFeature(b)
    flank <- getFeature(b, start_flank = input$flank, end_flank = input$flank)
    
    # Filtering names of bigwig files
    
    fbw <- bigwig_files[grepl("\\.f\\.bw$", bigwig_file_names)]
    rbw <- bigwig_files[grepl("\\.r\\.bw$", bigwig_file_names)]
    names(fbw) <- sub("\\.f\\.bw$", "", bigwig_file_names[grepl("\\.f\\.bw$", bigwig_file_names)])
    names(rbw) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
    
    # Import bigwig files as a list
    
    bwf <- importBWlist(fbw, names(fbw), selection = flank)
    bwr <- importBWlist(rbw, names(rbw), selection = flank)
    
    grl <- list("features" = features)
    
    matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive())
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
  
  max_ylim_reactive <- reactive({
    input$maxylim
  })
  

  # Defining the windows and axis labels
  
  wins_reactive <- reactive({
    req(flank_reactive())
    c("Upstream" = flank_reactive() / windowsize_reactive(), 
      "Feature" = (2 * flank_reactive()) / windowsize_reactive(), 
      "Downstream" = flank_reactive() / windowsize_reactive())
  })
  
  axis_labels_reactive <- reactive({
    req(flank_reactive())
    x <- flank_reactive()
    c(paste0("-", x, "b"), "TSS", "TES", paste0("+", x, "b"))
  })
  
  
  # Creating the heatmaps
  
  output$enrichedHeatmapPlot <- renderPlot({
    req(matl())
    hml <- hmList(matl = matl(), wins = wins_reactive(), col_fun = col_fun_reactive(), axis_labels = axis_labels_reactive(), show_row_names = show_row_names_reactive(), min_quantile = min_quantile_reactive(), max_quantile = max_quantile_reactive(), ylim = c(0, max_ylim_reactive()))
    req(length(hml) > 0)
    combined_hm <- Reduce(`+`, hml)
    draw(combined_hm, merge_legend = TRUE)
  })
  
  output$downloadpng <- downloadHandler(
    filename = function() { paste("heatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)  # Adjust quality
      hml <- hmList(
        matl = matl(), 
        wins = wins_reactive(), 
        col_fun = col_fun_reactive(), 
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = min_quantile_reactive(), 
        max_quantile = max_quantile_reactive(), 
        ylim = c(0, max_ylim_reactive())
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
      pdf(file, width = 13.5, height = 9)  # Adjust quality
      hml <- hmList(
        matl = matl(), 
        wins = wins_reactive(), 
        col_fun = col_fun_reactive(), 
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = min_quantile_reactive(), 
        max_quantile = max_quantile_reactive(), 
        ylim = c(0, max_ylim_reactive())
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