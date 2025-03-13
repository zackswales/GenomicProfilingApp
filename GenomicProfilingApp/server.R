server <- function(input, output, session) {
  
  ## Data processing
  
  output$regionfile1_name <- renderText({
    region1_names <- if (!is.null(input$Region1)) input$Region1$name else NULL
    region2_names <- if (!is.null(input$Region2)) input$Region2$name else NULL
    
    if (!is.null(region1_names) && !is.null(region2_names)) {
      paste(c(region1_names, region2_names), collapse = ", ")
    } else if (!is.null(region1_names)) {
      paste(region1_names, collapse = ", ")
    } else {
      "No files uploaded"
    }
  })
  
  ## Saving sequence data files
  
  saved_sequence <- reactiveValues(list = list())
  sequence_dir <- file.path("users", "s2274585", "sequence data files")
  
  read_sequence_files <- function() {
    if (dir.exists(sequence_dir)) {
      save_files <- list.files(sequence_dir, pattern = "\\.rds$")
      sequence_names <- gsub("\\.rds$", "", save_files)
      return(sequence_names)
    } else {
      return(character(0))
    }
  }
  
  observeEvent(input$savesequencedata, {
    sequence_name <- input$sequencenames # Get the desired file name from textInput
    sequence_dir <- file.path("users", "s2274585", "sequence data files")
    
    if (!dir.exists(sequence_dir)) {
      dir.create(sequence_dir, recursive = TRUE)
    }
    
    if (nchar(sequence_name) == 0) {
      showNotification("Please enter a sequence name.", type = "error")
      return()
    }
    
    # Copy uploaded files to the sequence data files directory
    uploaded_files <- input$Sequence1
    new_file_name <- paste0(input$sequencenames, ".", tools::file_ext(uploaded_files$name)) # Use the name from textInput
    new_file_path <- file.path(sequence_dir, new_file_name)
    file.copy(uploaded_files$datapath, new_file_path)
    
    # No RDS saving
    
    tryCatch({
      # ... any other actions you want to perform after saving ...
      showNotification(paste("Sequence", sequence_name, "saved successfully!"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving sequence", sequence_name, ":", e$message), type = "error")
    })
  })
  
  observeEvent(input$clearsequence, {
    if(dir.exists(sequence_dir)) {
      saved_files <- list.files(sequence_dir, full.names = TRUE)
      for(sequence_path in saved_files) {
        if(file.exists(sequence_path)) {
          tryCatch({
            file.remove(sequence_path)
            print(paste("File deleted:", sequence_path))
          }, error = function(e) {
            print(paste("Error deleting file:", sequence_path, e$message))
            showNotification(paste("Error deleting file:", basename(sequence_path), e$message), type = "error")
          })
        }
      }
      showNotification("Saved matrices cleared", type = "warning")
    }
  })
  
  read_saved_sequence <- function() {
    if (dir.exists(sequence_dir)) {
      saved_files <- list.files(sequence_dir, pattern = "\\.bw$")
      # Remove ".bw" extension
      sequence_names <- tools::file_path_sans_ext(saved_files)
      return(sequence_names)
    } else {
      return(character(0))
    }
  }
  
  
  saved_sequence_poll <- reactivePoll(
    intervalMillis = 500,
    session = session,
    checkFunc = function() {
      if (dir.exists(sequence_dir)) {
        saved_files <- list.files(sequence_dir, pattern = "\\.bw$", full.names = TRUE) # Look for .bw files
        if (length(saved_files) > 0) {
          # Return a vector of file modification times
          file.info(saved_files)$mtime
        } else {
          # Return a unique value when the directory is empty
          "empty"
        }
      } else {
        NULL
      }
    },
    valueFunc = read_saved_sequence
  )
  
  
  observe({
    sequence_names <- saved_sequence_poll()
    updateSelectInput(session, "sequencedatafiles",
                      choices = sequence_names,
                      selected = NULL)
  })
  
  
    
  output$pickgenome <- renderUI({
    if(input$databasefetch){
      selectInput(
        inputId = "genome",
        label = "Select genome:",
        choices = c("S.cerevisiae sacCer3" = 1),
        selected = 1
      )
    }
  })
  
  output$pickgroup <- renderUI({
    if(input$databasefetch){
      selectInput(
        inputId = "group",
        label = "Select group:",
        choices = c("genes" = 1),
        selected = 1
      )
    }
  })
  
  output$picktrack <- renderUI({
    if(input$databasefetch){
      selectInput(
        inputId = "track",
        label = "Select track:",
        choices = c("Structure" = 1, "Repeats" = 2, "Mutations" = 3),
        selected = 1
      )
    }
  })
  
  output$getannotation <- renderUI({
    if(input$databasefetch){
      actionButton(
        inputId = "fetchannotation",
        label = "Fetch annotation file"
      )
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
  
  # Fetching the data from UCSC
  
  genome_reactive <- reactive({
    value <- switch(input$genome,
                    "1" = "sacCer3")
  })
  
  track_reactive <- reactive({
    value <- switch(input$track,
                    "1" = "unipStruct",
                    "2" = "unipRepeat",
                    "3" = "unipMut")
  })
  
  observeEvent(input$fetchannotation, {
    showNotification("Fetching annotation...", type = "message")
  })
  
  observeEvent(input$fetchannotation, {
    # Add other debugging print statements here
  })
  
  gr <- eventReactive(input$fetchannotation, {
    
    fetch <- fetch_UCSC_track_data(genome_reactive(), track_reactive())
    
    GRanges(
      seqnames = fetch$chrom,
      ranges = IRanges(start = fetch$chromStart + 1, end = fetch$chromEnd),
      strand = fetch$strand,
      score = fetch$score,
      name = fetch$name,
      thickStart = as.integer(fetch$thickStart + 1), # Add +1 for 1-based coordinates
      thickEnd = fetch$thickEnd,
      itemRgb = as.integer(0),
      blockCount = fetch$blockCount,
      blockSizes = fetch$blockSizes,
      blockStarts = as.character(0)
    )
  })
  
  output$annotationname <- renderPrint({
    req(gr())
    paste(genome_reactive(), track_reactive())
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
  
  # Switching the inputs for the custom feature specification to required arguments
  
  startfeaturereactive <- reactive({
    switch(input$startfeature,
           "1" = "TSS",
           "2" = "TES",
           "3" = "Exon")
  })
  
  endfeaturereactive <- reactive({
    switch(input$endfeature,
           "1" = "TSS",
           "2" = "TES",
           "3" = "Exon")
  })
  
  startboundaryreactive <- reactive({
    switch(input$startboundary,
           "1" = "3prime",
           "2" = "5prime")
  })
  
  endboundaryreactive <- reactive({
    switch(input$endboundary,
           "1" = "3prime",
           "2" = "5prime")
  })
  
  saved_regions <- reactiveVal(list())
  
  observeEvent(input$getregion, {
    region_file <- if (!is.null(input$Region1)) input$Region1$datapath else NULL
    region2_file <- if (!is.null(input$Region2)) input$Region2$datapath else NULL
    
    if (!is.null(input$Region1)) {
      b <- readBed(region_file)
      combined_b <- if (!is.null(region2_file)) c(b, readBed(region2_file)) else b
    } else if (is.null(input$Region1)) {
      combined_b <- gr()
    }
    
    region_name <- input$regionname
    start_feature <- startfeaturereactive()
    end_feature <- endfeaturereactive()
    start_flank <- input$startflank
    end_flank <- input$endflank
    window_size <- input$winsize
    
    # Create a data frame for the new region
    new_region <- data.frame(
      "Region" = region_name,
      "Start Feature" = start_feature,
      "End Feature" = end_feature,
      "Upstream Flank" = start_flank,
      "Downstream Flank" = end_flank,
      "Window size" = window_size
    )
    
    current_regions <- saved_regions()
    saved_regions(c(current_regions, list(new_region)))
  })
  
  observeEvent(input$clearregions, {
    saved_regions(list()) # Reset saved_regions to an empty list
  })
  
  output$savedRegionsTable <- renderDT({
    if (input$getFeature == 4) {
      regions_list <- saved_regions()
      if (length(regions_list) > 0) {
        # Combine data frames from the list into a single data frame
        df <- do.call(rbind, regions_list)
        datatable(df)
      } else {
        datatable(data.frame(Message = "No regions saved yet."))
      }
    }
  })
  
  # Matrix list generation
  
  
  matl <- eventReactive(input$matrixgeneration, {
    req(input$Sequence1)
    
    withProgress(message = "Matrix Generation", value = 0, {
      
      regions <- list()
      
      region_file <- if (!is.null(input$Region1)) input$Region1$datapath else NULL
      region2_file <- if (!is.null(input$Region2)) input$Region2$datapath else NULL
      bigwig_files <- input$Sequence1$datapath
      bigwig_file_names <- input$Sequence1$name
      
      print(paste("BED/GTF File Path:", region_file))
      print(paste("BigWig File Path(s):", paste(bigwig_files, collapse = ", ")))
      
      # Import region files and get features
      
      if (!is.null(input$Region1)){
        b <- readBed(region_file)
        combined_b <- if (!is.null(region2_file)) c(b, readBed(region2_file)) else b
      } else if (is.null(input$Region1)) {
        combined_b <- gr()
      }
      
      
      incProgress(0.1, detail = "Reading region files...")
      
      if (input$getFeature == 1) {
        features <- getFeature(combined_b, start_feature = "TSS", end_feature = "TES")
        flank <- getFeature(combined_b, start_flank = input$flank, end_flank = input$flank)
      } else if (input$getFeature == 2) {
        features <- getFeature(combined_b, start_feature = "TSS", end_feature = "TSS")
        flank <- getFeature(combined_b, start_feature = "TSS", end_feature = "TSS", start_flank = input$flank, end_flank = input$flank)
      } else if (input$getFeature == 3) {
        features <- getFeature(combined_b, start_feature = "TES", end_feature = "TES")
        flank <- getFeature(combined_b, start_feature = "TES", end_feature = "TES", start_flank = input$flank, end_flank = input$flank)
      } else if (input$getFeature == 4) {
        flank <- getFeature(combined_b, start_feature = startfeaturereactive(), end_feature = endfeaturereactive(), start_flank = input$startflank, end_flank = input$endflank)
        feature_objects <- reactive({
          regions <- list()  # Initialize empty list
          
          # Upstream Flank
          if (input$customflanks == 1) {
            regions$up <- getFeature(
              combined_b, 
              start_flank = input$startflank, 
              end_feature = "TSS"
            )
          }
          
          # Exon numbers
          if (input$customexons == 1) {
            for (exon_num in seq(input$startexon, input$endexon)) {
              regions[[paste0("ex", exon_num)]] <- getFeature(
                combined_b, 
                start_feature = "Exon", 
                start_exon = exon_num, 
                start_exon_boundary = "5prime",
                end_feature = "Exon", 
                end_exon = exon_num, 
                end_exon_boundary = "3prime"
              )
            }
          }
          
          # First Intron
          if (input$customexons == 1 && input$customboundary == 1) {
            for (exon_num in seq(input$startexon, input$endexon - 1)) {
              regions[[paste0("in", exon_num)]] <- getFeature(
                combined_b, 
                start_feature = "Exon", 
                start_exon = exon_num,
                start_exon_boundary = "3prime",
                end_feature = "Exon",
                end_exon = exon_num + 1,
                end_exon_boundary = "5prime"
              )
            }
          }
          
          # Gene Body (remaining exons)
          if (input$customexons == 1 && input$customboundary == 1) {
            regions$body <- getFeature(
              combined_b, 
              start_feature = "Exon",
              start_exon = input$endexon,
              start_exon_boundary = "5prime"
            )
          }
          
          # Downstream Flank
          if (input$customflanks == 1) {
            regions$down <- getFeature(
              combined_b, 
              start_feature = "TES", 
              end_flank = input$endflank
            )
          }
          print("All regions:")
          print(regions)
          return(regions)
        })
      }
      
      incProgress(0.2, detail = "Getting features...")
      
      # Filtering names of bigwig files
      if (any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        fbw <- bigwig_files[grepl("\\.f\\.bw$", bigwig_file_names)]
        rbw <- bigwig_files[grepl("\\.r\\.bw$", bigwig_file_names)]
        
        if (length(fbw) > 0) {
          names(fbw) <- sub("\\.f\\.bw$", "", bigwig_file_names[grepl("\\.f\\.bw$", bigwig_file_names)])
          bwf <- importBWlist(fbw, names(fbw), selection = flank)
        }
        
        if (length(rbw) > 0 && length(fbw) > 0) {
          names(rbw) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
          bwr <- importBWlist(rbw, names(rbw), selection = flank)
        }
        
        if (length(fbw) <= 0 && length(rbw) > 0) {
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
      
      incProgress(0.5, detail = "Generating matrices")
      
      # Import bigwig files as a list
      if (any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          grl_reactive <- reactive({
            regions <- feature_objects()
            grl <- list()
            for (name in names(regions)) {
              grl[[name]] <- regions[[name]]
            }
            return(grl)
          })
          wins_reactive <- reactive({
            regions <- feature_objects()
            wins <- c()
            for (name in names(regions)) {
              wins[name] <- 10
            }
            return(wins)
          })
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl_reactive(), names = names(fbw), wins = wins_reactive(), strand = strand_reactive(), smooth = smooth_reactive())
        }
      } else {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          
          grl_reactive <- reactive({
            regions <- feature_objects()
            grl <- list()
            for (name in names(regions)) {
              grl[[name]] <- regions[[name]]
            }
            return(grl)
          })
          
          wins_reactive <- reactive({
            regions <- feature_objects()
            wins <- c()
            for (name in names(regions)) {
              wins[name] <- 10
            }
            return(wins)
          })
          
          matl_result <- matList(bwf = bwf, grl = grl_reactive(), names = names(fbw), wins = wins_reactive(), strand = strand_reactive(), smooth = smooth_reactive())
        }
      }
      
      incProgress(0.9, detail = "Almost finished...")
      return(matl_result)
    })
  })
  
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
  save_dir = file.path("users", "s2274585", "matrices")
  
  read_saved_matrices <- function() {
    if (dir.exists(save_dir)) {
      saved_files <- list.files(save_dir, pattern = "\\.rds$")
      matrix_names <- gsub("\\.rds$", "", saved_files)
      return(matrix_names)
    } else {
      return(character(0))
    }
  }
  
  saved_matrices_poll <- reactivePoll(
    intervalMillis = 1000, 
    session = session,
    checkFunc = function() {
      if (dir.exists(save_dir)) {
        file.info(save_dir)$mtime
      } else {
        NULL
      }
    },
    valueFunc = read_saved_matrices
  )
  
  observeEvent(input$savematrices, {
    req(matl())
    mat_list <- matl()
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
    for (mat_name in names(mat_list)) {
      file_path <- file.path(save_dir, paste0(mat_name, ".rds"))
      tryCatch({
        saveRDS(mat_list[[mat_name]], file = file_path)
        Sys.sleep(0.1)
        showNotification(paste("Matrix", mat_name, "saved successfully!"), type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving matrix", mat_name, ":", e$message), type = "error")
      })
    }
  })
  
  output$savedmatrices <- renderText({
    matrix_names <- saved_matrices_poll()
    if (length(matrix_names) == 0) {
      return("No matrices saved")
    } else {
      return(paste(matrix_names, collapse = ", "))
    }
  })
  
  observeEvent(input$clearmatrices, {
    if (dir.exists(save_dir)) {
      saved_files <- list.files(save_dir, full.names = TRUE)
      for (file_path in saved_files) {
        if (file.exists(file_path)) {
          tryCatch({
            file.remove(file_path)
            print(paste("File deleted:", file_path)) 
          }, error = function(e) {
            print(paste("Error deleting file:", file_path, e$message)) 
            showNotification(paste("Error deleting file:", basename(file_path), e$message), type = "error")
          })
        }
      }
      showNotification("Saved matrices cleared", type = "warning")
    }
    updateCheckboxGroupInput(session, "selectedmatrices", choices = character(0), selected = NULL)
  })
  
  observe({
    matrix_names <- saved_matrices_poll()
    updateCheckboxGroupInput(session, "selectedmatrices",
                             choices = matrix_names,
                             selected = NULL
    )
  })
  
  
  ### enrichedHeatmap Heatmaps
  
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
      mat_path <- file.path(save_dir, paste0(mat_name, ".rds")) # Construct the correct path
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
  
  observe ({
    if (input$heatmapplotbutton > 0) {
      showNotification("Plotting output...", type = "message")
    }
  })
  
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
  
  
  ## Download options for the enrichedHeatmap heatmaps
  
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
  
  
  ## ggplot Heatmaps
  
  # Reactive customisation options for the heatmaps
  
  ggshow_row_names_reactive <- reactive({
    input$ggshowrownames
  })
  
  ggheatmap_col_fun_reactive <- reactive({
    switch(input$ggheatmap_col_fun,
           "1" = "red",
           "2" = "bl2rd",
           "3" = "red0")
  })
  
  ggheatmap_min_quantile_reactive <- reactive({
    input$ggheatmapquantiles[1]
  })
  
  ggheatmap_max_quantile_reactive <- reactive({
    input$ggheatmapquantiles[2]
  })
  
  ggmax_ylim_reactive <- reactive({
    input$ggmaxylim
  })
  
  ggwins_reactive <- reactive({
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
  
  
  ggaxis_labels_reactive <- reactive({
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
  
  
  # ggplot Heatmaps
  
  observe ({
    if (input$ggplotheatmapplotbutton > 0) {
      showNotification("Plotting output...", type = "message")
    }
  })
  
  output$ggplotHeatmapPlot <- renderPlot({
    req(gghml())
    req(length(gghml()) > 0)
    combined_plot <- wrap_plots(gghml(), nrow = 1) 
    print(combined_plot)
  })
  
  gghml <- eventReactive(input$ggplotheatmapplotbutton, {
    req(selected_matrices_reactive())
    
    if(isTRUE(input$split)) {
      gghml <- gghmList(
        matl = selected_matrices_reactive(),
        wins = ggwins_reactive(),
        split = split_reactive(),
        split_cols = split_cols_reactive(),
        col_fun = ggheatmap_col_fun_reactive(),
        axis_labels = ggaxis_labels_reactive(),
        show_row_names = ggshow_row_names_reactive(),
        min_quantile = ggheatmap_min_quantile_reactive(),
        max_quantile = ggheatmap_max_quantile_reactive()
      )
    } else {
      gghml <- gghmList(
        matl = selected_matrices_reactive(),
        wins = ggwins_reactive(),
        col_fun = ggheatmap_col_fun_reactive(),
        axis_labels = ggaxis_labels_reactive(),
        show_row_names = ggshow_row_names_reactive(),
        min_quantile = ggheatmap_min_quantile_reactive(),
        max_quantile = ggheatmap_max_quantile_reactive()
      )
    }
    
    return(gghml)
  })
  
  output$ggplotheatmapdownloadpng <- downloadHandler(
    filename = function() { paste("ggplotheatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)
      gghml <- gghmList(
        matl = selected_matrices_reactive(),
        wins = ggwins_reactive(),
        col_fun = ggheatmap_col_fun_reactive(),
        axis_labels = ggaxis_labels_reactive(),
        show_row_names = ggshow_row_names_reactive(),
        min_quantile = ggheatmap_min_quantile_reactive(),
        max_quantile = ggheatmap_max_quantile_reactive()
      )
      req(length(gghml()) > 0)
      combined_plot <- wrap_plots(gghml(), nrow = 1) 
      
      print(combined_plot)
      dev.off()
    }
  )
  
  output$ggplotheatmapdownloadpdf <- downloadHandler(
    filename = function() { paste("ggplotheatmap_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)
      gghml <- gghmList(
        matl = selected_matrices_reactive(),
        wins = ggwins_reactive(),
        col_fun = ggheatmap_col_fun_reactive(),
        axis_labels = ggaxis_labels_reactive(),
        show_row_names = ggshow_row_names_reactive(),
        min_quantile = ggheatmap_min_quantile_reactive(),
        max_quantile = ggheatmap_max_quantile_reactive()
      )
      req(length(gghml()) > 0)
      combined_plot <- wrap_plots(gghml(), nrow = 1) 
      
      print(combined_plot)
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
  
  observe ({
    if (input$averageprofileplotbutton > 0) {
      showNotification("Plotting output...", type = "message")
    }
  })
  
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