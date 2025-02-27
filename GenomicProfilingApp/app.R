library(shiny)
library(bslib)

source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Define UI for application with three panels
ui <- page_navbar(
  title = "Interactive Genomic Profiling",
  bg = "#E46303",
  inverse = TRUE,
  
  # Panel for home page with directions for use of application
  nav_panel(
    title = "Home",
    card(
      "Interactive genomic profiling tool to probe internal genomic features. 
       Begin by navigating to the data processing page for file uploading and 
       matrix computation before moving to the visualisations page to generate outputs."
    )
  ),
  
  # Panel for data processing and matrix formation
  nav_panel(
    title = "Data Processing",
    navset_tab(
      
      # File Upload Panel
      nav_panel(
        title = "File Upload",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "Matrix Generation",
            fileInput(
              inputId = "Region1",
              label = "Upload a region file (.bed/.gtf)",
              accept = c(".bed", ".gtf"),
              multiple = FALSE
            ),
            fileInput(
              inputId = "Sequence1",
              label = "Upload sequence data files (.bw)",
              accept = ".bw",
              multiple = TRUE
            ),
            radioButtons(
              inputId = "strand",
              label = "Select strandedness of data:",
              choices = list("Unstranded" = 1, "Forward" = 2, "Reverse" = 3),
              selected = 1
            ),
            helpText("If all files have the same strandedness - select unstranded")
          ),
          card(
            full_screen = FALSE,
            card_header("Uploaded Region Files"),
            card_body(textOutput("regionfile1_name"))
          ),
          card(
            full_screen = FALSE,
            card_header("Uploaded Sequence Data Files"),
            card_body(textOutput("seqfile1_name"))
          )
        )
      ),
      
      # Feature Specification Panel
      nav_panel(
        title = "Feature Specification",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "Feature Specification",
            radioButtons(
              inputId = "getFeature",
              label = "Specify feature of interest:",
              choices = list("Full gene" = 1, "TSS" = 2, "TES" = 3)
            ),
            conditionalPanel(
              condition = "input.getFeature == 1",
              sliderInput(
                inputId = "windowsize",
                label = "Select window size for signal aggregation:",
                min = 1,
                max = 20,
                value = 1
              ),
              helpText("Warning: flank value must be divisible by window size")
            ),
            numericInput(
              inputId = "flank",
              label = "Specify flank around feature:",
              value = 20,
              step = 10
            )
          )
        )
      ),
      
      # Further Customisation Panel
      nav_panel(
        title = "Further Customisation",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "Further Customisation",
            checkboxInput(
              inputId = "smooth",
              label = "Smoothen data",
              value = FALSE
            ),
            actionButton(
              inputId = "matrixgeneration",
              label = "Generate Matrices"
            ),
            helpText("Warning: only click once as matrix generation takes a few seconds to complete")
          ),
          card(
            full_screen = FALSE,
            card_header("Generated Matrices"),
            card_body(textOutput("matrixnames"))
          )
        )
      ),
      
      # Saved Matrices Panel
      nav_panel(
        title = "Saved Matrices",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "Saved Matrices"
          )
        )
      )
    )
  ),
  
  # Panel for visualisations
  nav_panel(
    title = "Visualisations",
    navset_tab(
      
      # Heatmap Panel
      nav_panel(
        title = "Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE, 
            title = "Heatmap Customisation",
            radioButtons(
              inputId = "heatmap_col_fun",
              label = "Select colour scheme:",
              choices = list("White to red" = 1, "Blue to red" = 2, "Red scale" = 3),
              selected = 1
            ),
            sliderInput(
              inputId = "heatmapquantiles",
              label = "Set quantile range:",
              min = 0,
              max = 0.99,
              value = c(0, 0.99)
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
            actionButton(
              inputId = "heatmapplotbutton",
              label = "Plot Output"
            ),
            helpText("Output plotting will take a few seconds"),
            helpText("Use download buttons after clicking Plot Output"),
            downloadButton(
              outputId = "heatmapdownloadpng",
              label = "Download as .png"
            ),
            downloadButton(
              outputId = "heatmapdownloadpdf",
              label = "Download as .pdf"
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Heatmap Visualisation"),
            card_body(plotOutput("enrichedHeatmapPlot"))
          )
        )
      ),
      
      # Average Profile Plot Panel
      nav_panel(
        title = "Average Profile Plot",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE, 
            title = "Average Profile Plot Customisation",
            textInput(
              inputId = "plottitle",
              label = "Enter plot title:",
              value = "Average profile plot"
            ),
            sliderInput(
              inputId = "averageprofilequantiles",
              label = "Set quantile range:",
              min = 0,
              max = 1,
              value = c(0, 1)
            ),
            sliderInput(
              inputId = "alpha",
              label = "Select alpha value:",
              min = 0,
              max = 1,
              value = 1
            ),
            actionButton(
              inputId = "averageprofileplotbutton",
              label = "Plot Output"
            ),
            helpText("Output plotting will take a few seconds"),
            helpText("Use download buttons after clicking Plot Output"),
            downloadButton(
              outputId = "averageprofiledownloadpng",
              label = "Download as .png"
            ),
            downloadButton(
              outputId = "averageprofiledownloadpdf",
              label = "Download as .pdf"
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Average Profile Plot"),
            card_body(plotOutput("averageprofileplot"))
          )
        )
      )
    )
  )
)


# Define server logic required to compute matrices and plot outputs
server <- function(input, output) {
  
  ## Data processing
  
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
  
  smooth_reactive <- reactive({
    input$smooth
  })
  
  windowsize_reactive <- reactive({
    input$windowsize
  })
  
  observe({
    if (input$flank %% input$windowsize != 0) {
      showNotification("Flank value must be divisible by window size!", type = "error")
    }
  })
  
  observe({
    if (input$matrixgeneration > 0) {
      showNotification("Generating matrices...", type = "message")
    }
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
    
    if(input$getFeature == 1){
      features <- getFeature(b)
      flank <- getFeature(b, start_flank = input$flank, end_flank = input$flank)
    }
    
    if(input$getFeature == 2){
      features <- getFeature(b, start_feature = "TSS", end_feature = "TSS")
      flank <- getFeature(b, start_feature = "TSS", end_feature = "TSS", start_flank = input$flank, end_flank = input$flank)
    }
    
    if(input$getFeature == 3){
      features <- getFeature(b, start_feature = "TES", end_feature = "TES")
      flank <- getFeature(b, start_feature = "TES", end_feature = "TES", start_flank = input$flank, end_flank = input$flank)
    }
    
    # Filtering names of bigwig files
    
    if(any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))){
      fbw <- bigwig_files[grepl("\\.f\\.bw$", bigwig_file_names)]
      rbw <- bigwig_files[grepl("\\.r\\.bw$", bigwig_file_names)]
      names(fbw) <- sub("\\.f\\.bw$", "", bigwig_file_names[grepl("\\.f\\.bw$", bigwig_file_names)])
      names(rbw) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
    }
    
    else {
      fbw <- bigwig_files[grepl("\\.bw$", bigwig_file_names)]
      names(fbw) <- sub("\\.bw$", "", bigwig_file_names[grepl("\\.bw$", bigwig_file_names)])
    }
    
    
    # Import bigwig files as a list
    
    bwf <- importBWlist(fbw, names(fbw), selection = flank)
    bwr <- importBWlist(rbw, names(rbw), selection = flank)
    
    grl <- list("features" = features)
    
    if(any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))){
      if(input$getFeature == 1){
        matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        return(matl)
      } else if(input$getFeature == 2 || input$getFeature == 3){
        matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
      return(matl)
    }
    } else {
      if(input$getFeature == 1){
        matl <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        return(matl)
      } else if(input$getFeature == 2 || input$getFeature == 3){
        matl <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
      return(matl)
      }
    }
  })
  
  
  output$matrixnames <- renderText({
    req(matl())
    mat_names <- names(matl())
    if(is.null(mat_names) || length(mat_names) == 0) {
      return("No matrices generated")
    }
    paste(mat_names, collapse = ", ")
  })
  
  ### Heatmaps
  
  # Creating reactive objects for the heatmaps customisation options
  
  show_row_names_reactive <- reactive({
    input$showrownames
  })
  
  heatmap_col_fun_reactive <- reactive({
    switch(input$heatmap_col_fun,
           "1" = "red",
           "2" = "bl2rd",
           "3" = "red0")
  })
  
  heatmap_min_quantile_reactive <- reactive({
    input$heatmapquantiles[1]
  })
  
  heatmap_max_quantile_reactive <- reactive({
    input$heatmapquantiles[2]
  })
  
  max_ylim_reactive <- reactive({
    input$maxylim
  })
  
  wins_reactive <- reactive({
    req(flank_reactive(), windowsize_reactive())
    
    if (input$getFeature == 1) {
      return(c(
        "Upstream" = flank_reactive() / windowsize_reactive(), 
        "Feature" = (2 * flank_reactive()) / windowsize_reactive(), 
        "Downstream" = flank_reactive() / windowsize_reactive()
      ))
    } else if (input$getFeature == 2 || input$getFeature == 3) {
      return(c(
        "Upstream" = flank_reactive() / windowsize_reactive(), 
        "Feature" = 1, 
        "Downstream" = flank_reactive() / windowsize_reactive()
      ))
    } else {
      return(NULL)
    }
  })
  
  
  axis_labels_reactive <- reactive({
    req(flank_reactive(), input$getFeature)
    
    if(input$getFeature == 1){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", "TES", paste0("+", x, "b"))
    }
    
    else if(input$getFeature == 2){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", paste0("+", x, "b"))
    }
    
    else if(input$getFeature == 3){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TES", paste0("+", x, "b"))
    }
  })
  
  
  # Creating the heatmaps
  
  output$enrichedHeatmapPlot <- renderPlot({
    req(hml())
    req(length(hml()) > 0)
    combined_hm <- Reduce(`+`, hml())
    draw(combined_hm, merge_legend = TRUE)
  })
  
  hml <- eventReactive(input$heatmapplotbutton, {
    req(matl())
    hml <- hmList(matl = matl(), wins = wins_reactive(), col_fun = heatmap_col_fun_reactive(), axis_labels = axis_labels_reactive(), show_row_names = show_row_names_reactive(), min_quantile = heatmap_min_quantile_reactive(), max_quantile = heatmap_max_quantile_reactive(), ylim = c(0, max_ylim_reactive()))
    return(hml)
  })
 
  ## Download options for the heatmaps
   
  output$heatmapdownloadpng <- downloadHandler(
    filename = function() { paste("heatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)
      hml <- hmList(
        matl = matl(), 
        wins = wins_reactive(), 
        col_fun = heatmap_col_fun_reactive(), 
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = heatmap_min_quantile_reactive(), 
        max_quantile = heatmap_max_quantile_reactive(), 
        ylim = c(0, max_ylim_reactive())
      )
      
      req(length(hml) > 0)
      combined_hm <- Reduce(`+`, hml)
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  
    }
  )
  
  output$heatmapdownloadpdf <- downloadHandler(
    filename = function() { paste("heatmap_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)
      hml <- hmList(
        matl = matl(), 
        wins = wins_reactive(), 
        col_fun = heatmap_col_fun_reactive(), 
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(), 
        min_quantile = heatmap_min_quantile_reactive(), 
        max_quantile = heatmap_max_quantile_reactive(), 
        ylim = c(0, max_ylim_reactive())
      )
      
      req(length(hml) > 0)
      combined_hm <- Reduce(`+`, hml)
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  
    }
  )
  
  ## Average profile plot
  
  ## Reactive customisation options for the average profile plots
  
  title_reactive <- reactive({
    input$plottitle
  })
  
  averageprofile_min_quantile_reactive <- reactive({
    input$averageprofilequantiles[1]
  })
  
  averageprofile_max_quantile_reactive <- reactive({
    input$averageprofilequantiles[2]
  })
  
  alpha_reactive <- reactive({
    req(input$alpha)
    input$alpha
  })
  
  pal = c(RColorBrewer::brewer.pal(n = 9,name = "Set1"))
  pal2 = c(wes_palette("Darjeeling2")[2],wes_palette("Zissou1")[1],wes_palette("Darjeeling1")[4],wes_palette("Darjeeling1")[3])
  colmap = c(pal[c(1:5,7:9)])
  
  breaks_reactive <- reactive({
    if(input$getFeature == 1){
      x_range <- ncol(matl()[[1]])
      break_positions <- c(0,0.25,0.75,1) * x_range
    }
    else if(input$getFeature == 2){
      x_range <- ncol(matl()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
    }
    else if(input$getFeature == 3){
      x_range <- ncol(matl()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
    }
    
    return(break_positions)
  })
  
  labels_reactive <- reactive({
    req(flank_reactive(), input$getFeature)
    
    if(input$getFeature == 1){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", "TES", paste0("+", x, "b"))
    }
    
    else if(input$getFeature == 2){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", paste0("+", x, "b"))
    }
    
    else if(input$getFeature == 3){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TES", paste0("+", x, "b"))
    }
  })

  
  
  # Creating the average profile plot
  
  output$averageprofileplot <- renderPlot({
    req(average_profile())
    average_profile()
  })
  
  
  average_profile <- eventReactive(input$averageprofileplotbutton, {
    req(matl())
    average_profile <- mplot(matl = matl(), colmap = colmap, title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
    return(average_profile)
  })
  
  ## Downloading average profile plots
  
  output$averageprofiledownloadpng <- downloadHandler(
    filename = function() { paste("averageprofileplot_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)  
      average_profile <- mplot(matl = matl(), colmap = colmap, title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
      
      print(average_profile)
      dev.off() 
    }
  )
  
  output$averageprofiledownloadpdf <- downloadHandler(
    filename = function() { paste("averageprofileplot_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)  
      average_profile <- mplot(matl = matl(), colmap = colmap, title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
      
      print(average_profile)
      dev.off()  
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)