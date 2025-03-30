server <- function(input, output, session) {
  tab_login$server(input, output, session)
  
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
  
  observe({
    req(session$user$userData)
    sequence_dir <- file.path(session$userData$user_dir, "sequence data files")
    
    if(!dir.exists(sequence_dir)) {
      dir.create(sequence_dir, recursive = TRUE)
      print(paste("Created directory:", sequence_dir))
    }
  })
  
  get_sequence_dir <- function() {
    if(!is.null(session$userData$user_dir)){
      return(file.path(session$userData$user_dir, "sequence data files"))
    } else {
      return(NULL)
    }
  }
  
  saved_sequence_poll <- reactivePoll(
    intervalMillis = 1000,
    session = session,
    checkFunc = function() {
      sequence_dir <- get_sequence_dir()
      if(!is.null(sequence_dir) && dir.exists(sequence_dir)){
        return(file.info(sequence_dir)$mtime)
      } else {
        return(Sys.time())
      }
    },
    valueFunc = function() {
      sequence_dir <- get_sequence_dir()
      if(!is.null(sequence_dir) && dir.exists(sequence_dir)) {
        sequence_files <- list.files(sequence_dir, pattern = "\\.rds$")
        return(gsub("\\.rds$", "", sequence_files))
      } else {
        return(character(0))
      }
    }
  )
  
  observeEvent(input$savesequencedata, {
    req(input$Sequence1, session$userData$user_dir)
    
    sequence_dir <- get_sequence_dir()
    
    if(!dir.exists(sequence_dir)) {
      dir.create(sequence_dir, recursive = TRUE)
      print(paste("Created directory:", sequence_dir))
    }
    
    tryCatch({
      filepath <- input$Sequence1$datapath
      filename <- input$Sequence1$name
      bigwig <- import.bw(filepath)
      RDS <- saveRDS(object = bigwig, file = file.path(sequence_dir, paste0(filename, ".rds")))
      showNotification(paste0("Sequence data file", "", filename, "", "saved successfully!"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving sequence data file", e$message), type = "error")
    })
  })
  
  observe({
    sequence_selection <- saved_sequence_poll()
    updateSelectInput(session, "sequencedatafiles",
                      choices = sequence_selection,
                      selected = NULL)
  })
  
  observeEvent(input$clearsequence, {
    req(session$userData$user_dir)
    sequence_dir <- get_sequence_dir()
    
      if(dir.exists(sequence_dir)){
        rds_files <- dir_ls(sequence_dir, glob = "*.rds")
        
        if(length(rds_files) > 0){
          file_delete(rds_files)
          showNotification("Sequence data files cleared successfully", type = "message")
        } else {
          showNotification("No sequence data files found in directory", type = "warning")
        }
      } else {
        showNotification("Directory does not exist", type = "error")
      }
  })
  
  
  
  
  
  ## Database annotation fetching
  
  output$pickgenome <- renderUI({
    if(input$databasefetch){
      selectInput(
        inputId = "genome",
        label = "Select genome:",
        choices = c("S.cerevisiae sacCer3" = "sacCer3", "Human hg38" = "hg38", "Mouse mm39" = "mm39"),
        selected = "sacCer3"
      )
    }
  })
  
  output$pickgroup <- renderUI({
    if(input$databasefetch){
      selectInput(
        inputId = "group",
        label = "Select group:",
        choices = c("Genes" = "genes"),
        selected = "genes"
      )
    }
  })
  
  output$picktrack <- renderUI({
    if (input$databasefetch) {
      if (!is.null(input$genome)) {  # Check if input$genome has a value
        genome_selected <- input$genome
        
        track_choices <- switch(genome_selected,
                                "sacCer3" = c("RefSeq All" = "ncbiRefSeq", "Ensembl Genes" = "ensGene"),
                                "hg38" = c("RefSeq All" = "ncbiRefSeq", "UCSC RefSeq" = "refGene", "Ensembl GENCODE V20 Genes" = "wgEncodeGencodeV22ViewGenes"),
                                "mm39" = c("RefSeq All" = "ncbiRefSeq", "UCSC RefSeq" = "refGene", "GENCODE VM36" = "knownGene")
        )
        
        selectInput(
          inputId = "track",
          label = "Select track:",
          choices = track_choices,
          selected = track_choices[1]
        )
      } else {
        # Handle the case where input$genome is NULL (e.g., return a default UI)
        return(NULL) # Or you can return a default select input
      }
    } else {
      return(NULL) # Handle the case where databasefetch is not clicked
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
  
  output$tsvupload <- renderUI({
    if (input$split) {
      fileInput(
        inputId = "tsvsplitting",
        label = "Upload .tsv file for splitting",
        accept = ".tsv",
        multiple = FALSE
      )
    }
  })
  
  output$splitselect <- renderUI({
    if(input$split){
      selectInput(
        inputId = "splitby",
        label = "Select column to split by:",
        choices = column_choices(),
        selected = column_choices()[1]
      )
    }
  })
  
  column_choices <- reactiveVal(NULL)
  
  observeEvent(input$tsvsplitting, {
    if(!is.null(input$tsvsplitting)) {
      annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      column_choices(colnames(annotation))
    } else {
      columns_choices(NULL)
    }
  })
  
  
  # Fetching the data from UCSC
  
  observeEvent(input$fetchannotation, {
    showNotification("Fetching annotation...", type = "message")
  })
  
  observeEvent(input$fetchannotation, {
    # Add other debugging print statements here
  })
  
  gr <- eventReactive(input$fetchannotation, {
    
    fetch <- fetch_UCSC_track_data(input$genome, input$track)
    columns_to_remove <- c("name", "chrom", "txStart", "txEnd", "strand")
    mcols <- fetch[, !(names(fetch) %in% columns_to_remove)]
    
    GRanges(
      seqnames = fetch$chrom,
      ranges = IRanges(start = fetch$txStart + 1, end = fetch$txEnd),
      strand = fetch$strand,
      name = fetch$name,
      cdsStart = fetch$cdsStart,
      cdsEnd = fetch$cdsEnd,
      exonCount = fetch$exonCount,
      exonStarts = fetch$exonStarts,
      exonEnds = fetch$exonEnds,
      score = fetch$score, 
      name2 = fetch$name2, 
      cdsStartStat = fetch$cdsStartStat,
      cdsEndStat = fetch$cdsEndStat, 
      exonFrames = fetch$exonFrames
    )
  })
  
  output$annotationname <- renderPrint({
    req(input$fetchannotation)
    observeEvent(input$fetchannotation, { 
      output$annotationname <- renderPrint({ 
        paste(input$genome, input$track)
      })
    })
    return(NULL)
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
  
  log_reactive <- reactive({
    input$log
  })
  
  windowsize_reactive <- reactive({
    input$windowsize
  })
  
  observe({
    req(input$flank)
    req(input$windowsize)
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
  
  startdirectionreactive <- reactive({
    switch(input$startdirection,
           "1" = "up",
           "2" = "down")
  })
  
  enddirectionreactive <- reactive({
    switch(input$enddirection,
           "1" = "up",
           "2" = "down")
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
  region_objects_list <- reactiveVal(list())
  wins_vector <- reactiveVal(c())
  
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
    start_exon <- if(input$startfeature == 3) input$startexon else NA
    start_exon_boundary <- if(input$startfeature == 3) input$startexonboundary else NA
    end_feature <- endfeaturereactive()
    end_exon <- if(input$endfeature == 3) input$endexon else NA
    end_exon_boundary <- if(input$endfeature == 3) input$endexonboundary else NA
    start_flank <- input$startflank
    end_flank <- input$endflank
    start_direction <- startdirectionreactive()
    end_direction <- enddirectionreactive()
    window_size <- input$winsize
    
    # Create a data frame for the new region
    new_region <- data.frame(
      "Region" = region_name,
      "Start Feature" = start_feature,
      "Start Exon" = start_exon,
      "Start Exon Boundary" = start_exon_boundary,
      "End Feature" = end_feature,
      "End Exon" = end_exon,
      "End Exon Boundary" = end_exon_boundary,
      "Start Flank" = start_flank,
      "End Flank" = end_flank,
      "Start Direction" = start_direction,
      "End Direction" = end_direction,
      "Window Size" = window_size
    )
    
    current_regions <- saved_regions()
    saved_regions(c(current_regions, list(new_region)))
    
    region_object <- getFeature(
      object = combined_b,
      start_feature = start_feature,
      start_exon = if(is.na(start_exon)) NULL else start_exon,
      start_exon_boundary = if(is.na(start_exon_boundary)) NULL else start_exon_boundary,
      end_feature = end_feature,
      end_exon = if(is.na(end_exon)) NULL else end_exon,
      end_exon_boundary = if(is.na(end_exon_boundary)) NULL else end_exon_boundary,
      start_flank = start_flank,
      end_flank = end_flank,
      start_direction = start_direction,
      end_direction = end_direction
    )
    
    
    current_region_objects <- region_objects_list() 
    current_region_objects[[region_name]] <- region_object 
    region_objects_list(current_region_objects)
    
    current_wins <- wins_vector()
    current_wins[region_name] <- window_size
    wins_vector(current_wins)
  })
  
  observeEvent(input$clearregions, {
    saved_regions(list())
    region_objects_list(list())
    wins_vector(c())
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
  
  observeEvent(input$testfiles, {
    sequence_dir <- get_sequence_dir()
    
    all_rds_files <- list.files(sequence_dir, pattern = "\\.rds$", full.names = TRUE)
    
    all_rds_files_no_ext <- all_rds_file_names_no_ext <- tools::file_path_sans_ext(basename(all_rds_files))
    
    selected_rds_files <- all_rds_files[all_rds_file_names_no_ext %in% input$sequencedatafiles]
    
    bigwig_files <- lapply(selected_rds_files, readRDS)
    names(bigwig_files) <- tools::file_path_sans_ext(basename(selected_rds_files))
    bigwig_file_names <- names(bigwig_files)
    
    sample1 <- readRDS("~/GenomicProfilingApp/GenomicProfilingApp/users/s2274585/sequence data files/sample1.f.bw.rds")
    print(sample1)
    print(bigwig_files[[1]])
    print(bigwig_files[[2]])
    print(bigwig_file_names)
    print(paste("selected rds files:", selected_rds_files))
  })
  
  
  # Matrix list generation
  
  
  matl <- eventReactive(input$matrixgeneration, {
    req(input$sequencedatafiles)
    
    withProgress(message = "Matrix Generation", value = 0, {
      
      regions <- list()
      
      sequence_dir <- get_sequence_dir()
      
      region_file <- if (!is.null(input$Region1)) input$Region1$datapath else NULL
      region2_file <- if (!is.null(input$Region2)) input$Region2$datapath else NULL
      
      all_rds_files <- list.files(sequence_dir, pattern = "\\.rds$", full.names = TRUE)
      
      all_rds_file_names_no_ext <- tools::file_path_sans_ext(basename(all_rds_files))
      
      selected_rds_files <- all_rds_files[all_rds_file_names_no_ext %in% input$sequencedatafiles]
      
      bigwig_files <- lapply(selected_rds_files, readRDS)
      names(bigwig_files) <- tools::file_path_sans_ext(basename(selected_rds_files))
      bigwig_file_names <- names(bigwig_files)
      
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
        flank <- getFeature(combined_b, start_flank = 500, end_flank = 500)
        region_objects <- region_objects_list() # Get the list of region objects
        grl <- lapply(names(region_objects), function(region_name) {
          region_objects[[region_name]]
        })
        names(grl) <- names(region_objects)
        find_min_rows_gr <- function(grl) {
          num_rows <- sapply(grl, length)
          min_rows_index <- which.min(num_rows)
          return(grl[[min_rows_index]])
        }
        
        min_rows_gr <- find_min_rows_gr(grl)
        
        grl <- imap(grl, function(gr, index) {
          names_subset_to <- min_rows_gr$name
          gr_subset <- gr[gr$name %in% names_subset_to]
          return(gr_subset)
        })
      }
      
      incProgress(0.2, detail = "Getting features...")
      
      # Filtering names of bigwig files
      if (any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        bwf <- bigwig_files[grepl("\\.f\\.bw$", bigwig_file_names)]
        bwr <- bigwig_files[grepl("\\.r\\.bw$", bigwig_file_names)]
        
        if (length(bwf) > 0) {
          names(bwf) <- sub("\\.f\\.bw$", "", bigwig_file_names[grepl("\\.f\\.bw$", bigwig_file_names)])
        }
        
        if (length(bwr) > 0 && length(bwf) > 0) {
          names(bwr) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
        }
        
        if (length(bwf) <= 0 && length(bwr) > 0) {
          bwf <- bwr
          names(bwf) <- sub("\\.r\\.bw$", "", bigwig_file_names[grepl("\\.r\\.bw$", bigwig_file_names)])
        }
      } else {
        bwf <- bigwig_files[grepl("\\.bw$", bigwig_file_names)]
        if (length(bwf) > 0) {
          names(bwf) <- sub("\\.bw$", "", bigwig_file_names[grepl("\\.bw$", bigwig_file_names)])
        }
      }
      
      incProgress(0.5, detail = "Generating matrices")
      
      # Import bigwig files as a list
      if (any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive(), log = log_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), wins = wins_vector(), strand = strand_reactive(), smooth = smooth_reactive())
        }
      } else {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(fbw), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), wins = wins_vector(), strand = strand_reactive(), smooth = smooth_reactive())
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
  
  # Saving Matrices
  # Reactive values to store saved matrices
  # Ensure user directory is set when logging in
  saved_matrices <- reactiveValues(list = list())
  
  observe({
    req(session$userData$user_dir)  # Ensure user directory exists
    save_dir <- file.path(session$userData$user_dir, "matrices")
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      print(paste("Created directory:", save_dir))  # Debugging statement
    }
  })
  
  # Function to safely get save_dir
  get_save_dir <- function() {
    if (!is.null(session$userData$user_dir)) {
      return(file.path(session$userData$user_dir, "matrices"))
    } else {
      return(NULL)  # Return NULL if user is not logged in
    }
  }
  
  # Reactive polling for saved matrices
  saved_matrices_poll <- reactivePoll(
    intervalMillis = 1000, 
    session = session,
    checkFunc = function() {
      save_dir <- get_save_dir()
      if (!is.null(save_dir) && dir.exists(save_dir)) {
        return(file.info(save_dir)$mtime)
      } else {
        return(Sys.time())  # Return current time to avoid NULL issue
      }
    },
    valueFunc = function() {
      save_dir <- get_save_dir()
      if (!is.null(save_dir) && dir.exists(save_dir)) {
        saved_files <- list.files(save_dir, pattern = "\\.rds$")
        return(gsub("\\.rds$", "", saved_files))
      } else {
        return(character(0))
      }
    }
  )
  
  # Event to handle saving matrices
  observeEvent(input$savematrices, {
    req(matl(), session$userData$user_dir)  # Ensure matrices exist & user is logged in
    
    save_dir <- get_save_dir()
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      print(paste("Created directory:", save_dir))  # Debugging statement
    }
    
    mat_list <- matl()
    
    for (mat_name in names(mat_list)) {
      file_path <- file.path(save_dir, paste0(mat_name, ".rds"))
      tryCatch({
        saveRDS(mat_list[[mat_name]], file = file_path)
        Sys.sleep(0.1)
        showNotification(paste("Matrix", mat_name, "saved successfully!"), type = "message")
        print(paste("Saved matrix:", mat_name, "to", file_path))  # Debugging statement
      }, error = function(e) {
        showNotification(paste("Error saving matrix", mat_name, ":", e$message), type = "error")
        print(paste("Error saving matrix:", mat_name, e$message))  # Debugging statement
      })
    }
  })
  
  # Display saved matrices
  output$savedmatrices <- renderText({
    matrix_names <- saved_matrices_poll()
    if (length(matrix_names) == 0) {
      return("No matrices saved")
    } else {
      return(paste(matrix_names, collapse = ", "))
    }
  })
  
  # Clear saved matrices
  observeEvent(input$clearmatrices, {
    req(session$userData$user_dir)
    save_dir <- get_save_dir()
    
    if (dir.exists(save_dir)) {
      saved_files <- list.files(save_dir, full.names = TRUE)
      for (file_path in saved_files) {
        if (file.exists(file_path)) {
          tryCatch({
            file.remove(file_path)
            print(paste("File deleted:", file_path))  # Debugging statement
          }, error = function(e) {
            print(paste("Error deleting file:", file_path, e$message)) 
            showNotification(paste("Error deleting file:", basename(file_path), e$message), type = "error")
          })
        }
      }
      showNotification("Saved matrices cleared", type = "message")
    }
    
    updateCheckboxGroupInput(session, "selectedmatrices", choices = character(0), selected = NULL)
  })
  
  # Update UI with saved matrices
  observe({
    matrix_names <- saved_matrices_poll()
    updateCheckboxGroupInput(session, "selectedmatrices",
                             choices = matrix_names,
                             selected = NULL)
  })
  
  
  
  
  
  ### enrichedHeatmap Heatmaps
  
  # Creating reactive objects for the heatmaps customisation options
  
  # Creating the conditional splitting object
  split_reactive <- reactive({
    req(input$split, input$tsvsplitting$datapath, input$splitby, matl()[[1]])
    
    annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(annotation)[1] <- "gene_name"
    
    split_col <- input$splitby
    if (!split_col %in% colnames(annotation)) {
      stop("Selected column not found in annotation file.")
    }
    
    Anno <- data.frame(name = rownames(matl()[[1]])) %>%
      left_join(data.frame(name = annotation$gene_name, split_val = factor(annotation[[split_col]])), by = "name") %>%
      column_to_rownames("name")
    
    return(Anno)
  })
  
  
  split_cols_reactive <- reactive({
    req(split_reactive())
    unique_families <- unique(split_reactive()$split_val)
    
    if (length(unique_families) == 2) {
      # Create named vector directly with backticks
      colors <- c("#DA4167", "#083D77")
      names(colors) <- unique_families[1:2]
      list(split_val = colors)
    } else if (length(unique_families) > 2) {
      # Assign colors from a palette for more than two unique values
      color_palette <- grDevices::rainbow(length(unique_families))
      names(color_palette) <- unique_families
      names(color_palette) <- paste0("`", names(color_palette), "`") # Add backticks to the names.
      list(split_val = color_palette)
    } else {
      # Handle cases with 0 or 1 unique values
      list(split_val = c("default" = "gray"))
    }
  })
  
  # Creating other reactive objects for heatmap customisation
  
  selected_matrices_reactive <- reactive({
    req(input$selectedmatrices, session$userData$user_dir)  # Ensure user directory exists
    
    # Define the save directory dynamically for the logged-in user
    save_dir <- file.path(session$userData$user_dir, "matrices")
    
    selected_mats <- lapply(input$selectedmatrices, function(mat_name) {
      mat_path <- file.path(save_dir, paste0(mat_name, ".rds"))  # Construct the correct path
      if (file.exists(mat_path)) {
        return(readRDS(mat_path))
      } else {
        return(NULL)
      }
    })
    
    selected_mats <- selected_mats[!sapply(selected_mats, is.null)]  # Remove NULL values
    names(selected_mats) <- input$selectedmatrices[!sapply(selected_mats, is.null)]
    
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
      return(wins_vector())
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
  
  observeEvent(input$heatmapplotbutton, {
    if(length(input$selectedmatrices) == 0){
      showNotification("Matrices need to be selected before plotting", type = "warning")
    } else {
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
        row_km = imput$row_km,
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
        ylim = c(0, max_ylim_reactive()),
        row_km = input$row_km
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
  
  split_average_reactive <- reactive({
    req(input$split)
    region_file <- input$Region1$datapath
    annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(annotation) <- c("gene_name", "family")
    Anno <- data.frame(name = c(rownames(matl()[[1]]))) |>
      left_join(data.frame(name = annotation$gene_name, Family = factor(annotation$family)), by = "name") |>
      column_to_rownames("name")
    return(Anno)
  })
  
  split_average_cols_reactive <- reactive({
    req(split_average_reactive())
    unique_families <- unique(split_reactive()$Family)
    if (length(unique_families) == 2) {
      # Create named vector directly with backticks
      colors <- c("#DA4167", "#083D77")
      names(colors) <- unique_families[1:2]
      if(input$colourby == "Family"){
        list(Family = colors)
      } else if(input$colourby == "Sample"){
        list(Sample = colors)
      }
    } else if (length(unique_families) > 2) {
      # Assign colors from a palette for more than two unique values
      color_palette <- grDevices::rainbow(length(unique_families))
      names(color_palette) <- unique_families
      names(color_palette) <- paste0("`",names(color_palette),"`") # Add backticks to the names.
      if(input$colourby == "Family"){
        list(Family = color_palette)
      } else if(input$colourby == "Sample"){
        list(Sample = color_palette)
      }
    } else {
      if(input$colourby == "Family"){
        list(Family = c("default" = "gray"))
      } else if(input$colourby == "Sample"){
        list(Sample = c("default" = "gray"))
      }
    }
  })
  
  facet_reactive <- reactive({
    switch(input$facetby,
           "1" = "Sample",
           "2" = "Family",
           "3" = c("Sample", "Family"))
  })
  
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
  
  unit_reactive <- reactive({
    input$unit
  })
  
  feature_reactive <- reactive({
    input$feature
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
      average_profile <- mplot(matl = selected_matrices_reactive(), split = split_reactive(),colmap = split_average_cols_reactive()[[input$colourby]], feature = feature_reactive(), unit = unit_reactive(), title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive(), colour_by = input$colourby,facet = facet_reactive(), facet_scale = "free",facet_independent = T,facet_type="grid")
      return(average_profile)
    } else {
      average_profile <- mplot(matl = selected_matrices_reactive(), colmap = colmap, feature = feature_reactive(), unit = unit_reactive(), title = title_reactive(), min_quantile = averageprofile_min_quantile_reactive(), max_quantile = averageprofile_max_quantile_reactive(), alpha = alpha_reactive(), breaks = breaks_reactive(), labels = labels_reactive())
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