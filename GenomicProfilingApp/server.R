server <- function(input, output, session) {
  
  ## Data processing
  
  output$regionfile1_name <- renderText({
    if(!is.null(input$Region1)) {
      paste(input$Region1$name, collapse = ", ")
    } else {
      "No files uploaded"
    }
  })
  
  output$Region1splitting <- renderUI({
    if(input$split){
      textInput(
        inputId = "Region1splitting",
        label = "Splitting identifier in first region file"
      )
    }
  })
  
  output$conditionalRegion2 <- renderUI({
    if(input$split){
      fileInput(
        inputId = "Region2",
        label = "Upload second region file for splitting (.bed/.gtf)",
        accept = c(".bed", ".gtf"),
        multiple = FALSE
      )
    }
  })
  
  output$conditionalRegion2splitting <- renderUI({
    if(input$split){
      textInput(
        inputId = "Region2splitting",
        label = "Splitting identifier in second region file"
      )
    }
  })
  
  region_files_count <- reactive({
    req(input$Region1)
    length(input$Region1$datapath)
  })
  
  output$seqfile1_name <- renderText({
    if(!is.null(input$Sequence1)) {
      paste(input$Sequence1$name, collapse = ", ")
    } else {
      "No files uploaded"
    }
  })
  
  # Creating reactive objects for the matrix customisation options
  
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
    
    if(region_files_count() == 1){
      region_file <- input$Region1$datapath
      if(!is.null(input$Region2)){
        region2_file <- input$Region2$datapath
      }
      bigwig_files <- input$Sequence1$datapath
      bigwig_file_names <- input$Sequence1$name
      
      print(paste("BED/GTF File Path:", region_file))
      print(paste("BigWig File Path(s):", paste(bigwig_files, collapse = ", ")))
      
      # Import region files and get features
      
      b <- readBed(region_file)
      print(head(b))
      if(!is.null(input$Region2)){
        b2 <- readBed(region2_file)
        combined_b <- c(b, b2)
        print(head(b2))
      }
      
      if(input$getFeature == 1){
        if(!is.null(input$Region2)){
        features <- getFeature(combined_b, start_feature = "TSS", end_feature = "TES")
        flank <- getFeature(combined_b, start_flank = input$flank, end_flank = input$flank)
        } else {
          features <- getFeature(b, start_feature = "TSS", end_feature = "TES")
          flank <- getFeature(b, start_flank = input$flank, end_flank = input$flank)
        }
      }
      
      if(input$getFeature == 2){
        if(!is.null(input$Region2)){
        features <- getFeature(combined_b, start_feature = "TSS", end_feature = "TSS")
        flank <- getFeature(combined_b, start_feature = "TSS", end_feature = "TSS", start_flank = input$flank, end_flank = input$flank)
        } else {
          features <- getFeature(b, start_feature = "TSS", end_feature = "TSS")
          flank <- getFeature(b, start_feature = "TSS", end_feature = "TSS", start_flank = input$flank, end_flank = input$flank)
        }
      }
      
      if(input$getFeature == 3){
        if(!is.null(input$Region2)){
        features <- getFeature(combined_b, start_feature = "TES", end_feature = "TES")
        flank <- getFeature(combined_b, start_feature = "TES", end_feature = "TES", start_flank = input$flank, end_flank = input$flank)
        } else {
          features <- getFeature(b, start_feature = "TES", end_feature = "TES")
          flank <- getFeature(b, start_feature = "TES", end_feature = "TES", start_flank = input$flank, end_flank = input$flank)
        }
      }
      
      if(input$getFeature == 4){
        if(!is.null(input$Region2)){
        flank <- getFeature(combined_b, start_flank = input$flank, end_flank = input$flank)
        up <- getFeature(combined_b, start_flank = input$flank, end_feature = "TSS")
        ex1 <- getFeature(combined_b, end_feature = "Exon", end_exon = 1, end_exon_boundary = "3prime")
        in1 <- getFeature(combined_b, start_feature = "Exon", start_exon = 1,start_exon_boundary = "3prime",end_feature = "Exon",end_exon = 2,end_exon_boundary = "5prime")
        body <- getFeature(combined_b, start_feature = "Exon",start_exon = 2,start_exon_boundary = "5prime")
        down <- getFeature(combined_b, start_feature = "TES", end_flank = input$flank)
        } else {
          flank <- getFeature(b, start_flank = input$flank, end_flank = input$flank)
          up <- getFeature(b, start_flank = input$flank, end_feature = "TSS")
          ex1 <- getFeature(b, end_feature = "Exon", end_exon = 1, end_exon_boundary = "3prime")
          in1 <- getFeature(b, start_feature = "Exon", start_exon = 1,start_exon_boundary = "3prime",end_feature = "Exon",end_exon = 2,end_exon_boundary = "5prime")
          body <- getFeature(b, start_feature = "Exon",start_exon = 2,start_exon_boundary = "5prime")
          down <- getFeature(b, start_feature = "TES", end_flank = input$flank)
        }
      }
      
      # Filtering names of bigwig files
      
      if(any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        fbw <- bigwig_files[grepl("\\.f\\.bw$", bigwig_file_names)]
        rbw <- bigwig_files[grepl("\\.r\\.bw$", bigwig_file_names)]
        
        # Assign names ONLY if fbw or rbw is non-empty
        if (length(fbw) > 0) {
          names(fbw) <- sub("\\.f\\.bw$", "", bigwig_file_names[grepl("\\.f\\.bw$", bigwig_file_names)])
          bwf <- importBWlist(fbw, names(fbw), selection = flank)
        }
        
        if (length(rbw) > 0 && length(fbw) > 0) {
          names(rbw) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
          bwr <- importBWlist(rbw, names(rbw), selection = flank)
        }
        
        if (length(fbw) <= 0 && length(rbw) > 0){
          fbw <- rbw
          names(fbw) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
          bwf <- importBWlist(fbw, names(fbw), selection = flank)
        }
      } else {
        fbw <- bigwig_files[grepl("\\.bw$", bigwig_file_names)]
        if (length(fbw) > 0) {
          names(fbw) <- sub("\\.bw$", "", bigwig_file_names[grepl("\\.bw$", bigwig_file_names)])
          bwf <- importBWlist(fbw, names(fbw), selection = flank)
        }
      }
      
      
      
      # Import bigwig files as a list
      
      if(any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))){
        if(input$getFeature == 1){
          grl <- list("features" = features)
          matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
          return(matl)
        }
        if(input$getFeature == 2 || input$getFeature == 3){
          grl <- list("features" = features)
          matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
          return(matl)
        }
        if (input$getFeature == 4){
          grl <- list("up" = up, "exon1" = ex1, "intron1" = in1, "body" = body, "down" = down)
          wins <- c("up" = input$up, "exon1" = input$exon1, "intron1" = input$intron1, "body" = input$body, "down" = input$down)
          matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), wins = wins, strand = strand_reactive(), smooth = smooth_reactive())
          return(matl)
        }
        
        
      } else {
        if(input$getFeature == 1){
          matl <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
          return(matl)
        }
        if(input$getFeature == 2 || input$getFeature == 3){
          matl <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
          return(matl)
        }
        if(input$getFeature == 4){
          grl <- list("up" = up, "exon1" = ex1, "intron1" = in1, "body" = body, "down" = down)
          wins <- c("up" = input$up, "exon1" = input$exon1, "intron1" = input$intron1, "body" = input$body, "down" = input$down)
          matl <- matList(bwf = bwf, grl = grl, names = names(fbw), wins = wins, strand = strand_reactive(), smooth = smooth_reactive())
        }
      }
    }})
  
  
  output$matrixnames <- renderText({
    req(matl())
    mat_names <- names(matl())
    if(is.null(mat_names) || length(mat_names) == 0) {
      return("No matrices generated")
    }
    paste(mat_names, collapse = ", ")
  })
  
  # Saving Matrices
  
  saved_matrices <- reactiveValues(list = list())
  save_dir = tempdir()
  
  observeEvent(input$savematrices, {
    req(matl())
    mat_list <- matl()
    
    for (mat_name in names(mat_list)) {
      file_path <- file.path(save_dir, paste0(mat_name, ".rds"))
      saveRDS(mat_list[[mat_name]], file = file_path)
      saved_matrices$list[[mat_name]] <- file_path
    }
    
    showNotification("Matrices saved successfully!", type = "message")
  })
  
  output$savedmatrices <- renderText({
    if (length(saved_matrices$list) == 0) {
      return("No matrices saved")
    }
    paste(names(saved_matrices$list), collapse = ", ")
  })
  
  observeEvent(input$clearmatrices, {
    for (file_path in saved_matrices$list) {
      if (file.exists(file_path)) {
        file.remove(file_path)
      }
    }
    saved_matrices$list <- list()  # Clear list
    showNotification("Saved matrices cleared", type = "warning")
    updateCheckboxGroupInput(session, "selectedmatrices", choices = character(0), selected = NULL)
  })
  
  observe({
    updateCheckboxGroupInput(session, "selectedmatrices",
                             choices = names(saved_matrices$list),
                             selected = NULL
    )
  })
  
  ### Heatmaps
  
  # Creating reactive objects for the heatmaps customisation options
  
  # Creating the conditional splitting object
  split_reactive <- reactive({
    req(input$split)
    region_file <- input$Region1$datapath
    region2_file <- input$Region2$datapath
    b <- readBed(region_file)
    b2 <- readBed(region2_file)
    Anno <- data.frame(name = rownames(matl()[[1]])) |>
      left_join(data.frame(name = b$name, b_group = input$Region1splitting), by = "name") |>
      left_join(data.frame(name = b2$name, b2_group = input$Region2splitting), by = "name") |>
      mutate(Group = coalesce(b_group, b2_group)) |>
      select(name, Group) |>
      column_to_rownames("name")
    return(Anno)
  })
  
  split_cols_reactive <- reactive({
    anno_cols <- list(group = c(`input$Region1splitting` = "#DA4167", `input$Region2splitting` = "#083D77"))
    return(anno_cols)
  })
  
  # Creating other reactive objects for heatmap customisation

  selected_matrices_reactive <- reactive({
    req(input$selectedmatrices)
    
    selected_mats <- lapply(input$selectedmatrices, function(mat_name) {
      mat_path <- saved_matrices$list[[mat_name]]
      if (file.exists(mat_path)) {
        return(readRDS(mat_path))
      } else {
        return(NULL)
      }
    })
    selected_mats <- selected_mats[!sapply(selected_mats, is.null)]
    names(selected_mats) <- input$selectedmatrices[!sapply(selected_mats, function(x) is.null(x))]
    
    return(selected_mats)
  })
  
  
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
    }
    if (input$getFeature == 2 || input$getFeature == 3) {
      return(c(
        "Upstream" = flank_reactive() / windowsize_reactive(), 
        "Feature" = 1, 
        "Downstream" = flank_reactive() / windowsize_reactive()
      ))
    }
    if (input$getFeature == 4){
      wins <- c("up" = input$up, "exon1" = input$exon1, "intron1" = input$intron1, "body" = input$body, "down" = input$down)
      return(wins)
    }
    else {
      return(NULL)
    }
  })
  
  
  axis_labels_reactive <- reactive({
    req(flank_reactive(), input$getFeature)
    
    if(input$getFeature == 1){
      x <- flank_reactive()
      return(c(paste0("-", x, "b"), "TSS", "TES", paste0("+", x, "b")))
    }
    
    if(input$getFeature == 2){
      x <- flank_reactive()
      return(c(paste0("-", x, "b"), "TSS", paste0("+", x, "b")))
    }
    
    if(input$getFeature == 3){
      x <- flank_reactive()
      return(c(paste0("-", x, "b"), "TES", paste0("+", x, "b")))
    }
    
    if(input$getFeature == 4){
      x <- flank_reactive()
      return(c(paste0("-", x, "b"), paste0("+", x, "b")))
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
    req(selected_matrices_reactive())
    
    if (isTRUE(input$split)) {  # Ensure reactive context
      hml <- hmList(
        matl = selected_matrices_reactive(),
        wins = wins_reactive(),
        split = split_reactive(), 
        split_cols = split_cols_reactive(), 
        col_fun = heatmap_col_fun_reactive(),
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(),
        min_quantile = heatmap_min_quantile_reactive(), 
        max_quantile = heatmap_max_quantile_reactive(),
        ylim = c(0, max_ylim_reactive())
      )
    } else {
      hml <- hmList(
        matl = selected_matrices_reactive(),
        wins = wins_reactive(),
        col_fun = heatmap_col_fun_reactive(),
        axis_labels = axis_labels_reactive(), 
        show_row_names = show_row_names_reactive(),
        min_quantile = heatmap_min_quantile_reactive(), 
        max_quantile = heatmap_max_quantile_reactive(),
        ylim = c(0, max_ylim_reactive())
      )
    }
    
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
    if(input$getFeature == 2){
      x_range <- ncol(matl()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
    }
    if(input$getFeature == 3){
      x_range <- ncol(matl()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
    }
    if(input$getFeature == 4){
      x_range <- ncol(matl()[[1]])
      sum = input$up + input$exon1 + input$intron1 + input$body + input$down
      break_positions <- c(0, (input$up/sum), ((input$exon1/sum) + (input$up/sum)), ((input$intron1/sum) + (input$exon1/sum) + (input$up/sum)), ((input$body/sum) + (input$intron1/sum) + (input$exon1/sum) + (input$up/sum)), 1) * x_range
    }
    
    return(break_positions)
  })
  
  labels_reactive <- reactive({
    req(flank_reactive(), input$getFeature)
    
    if(input$getFeature == 1){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", "TES", paste0("+", x, "b"))
    }
    
    if(input$getFeature == 2){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TSS", paste0("+", x, "b"))
    }
    
    if(input$getFeature == 3){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "TES", paste0("+", x, "b"))
    }
    
    if(input$getFeature == 4){
      x <- flank_reactive()
      c(paste0("-", x, "b"), "Exon 1", "Intron 1", "Body", "TES", paste0("+", x, "b"))
    }
  })
  
  
  
  # Creating the average profile plot
  
  output$averageprofileplot <- renderPlot({
    req(average_profile())
    average_profile()
  })
  
  
  average_profile <- eventReactive(input$averageprofileplotbutton, {
    req(matl())
    if(input$split){
    average_profile <- mplot(matl = selected_matrices_reactive(), split = split_reactive(),colmap = colmap, title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
    return(average_profile)
    } else {
      average_profile <- mplot(matl = selected_matrices_reactive(), colmap = colmap, title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
      return(average_profile)
    }
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