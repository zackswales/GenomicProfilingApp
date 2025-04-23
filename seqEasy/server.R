library(genomation)
#library(tidyverse)
library(EnrichedHeatmap)
library(wesanderson)
library(furrr)
library(purrr)
library(circlize)
library(RColorBrewer)
library(plotrix)
library(ggh4x)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(stats)
library(viridis)
library(tibble)
library(cowplot)
library(ComplexHeatmap)
library(shiny)
library(bslib)
library(shinyjs)
library(glue)
library(ggplot2)
library(patchwork)
library(DT)
library(scales)
library(fs)
library(seqEasyFunctions)
library(readr)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1024^3) # 1 GB
  tab_login$server(input, output, session)
  
  ## Logic for downloading test files
  
  output$testfiledownload <- downloadHandler(
    filename = function() {
      paste0("seqEasy_test_files", Sys.Date(), ".zip")
    },
    content = function(file) {
      withProgress(message = 'Creating download package...', value = 0, {
        tmpdir <- tempdir()
        test_files <- list.files("data", full.names = TRUE)
        num_files <- length(test_files)
        
        for (i in 1:num_files) {
          file.copy(test_files[i], tmpdir, overwrite = TRUE)
          incProgress(1/num_files, detail = paste("Copying file", i, "of", num_files))
        }
        
        zip(zipfile = file, files = file.path(tmpdir, basename(test_files)), flags = "-j")
      })
    },
    contentType = "application/zip"
  )
  
  
  # Resetting all inputs when a new user logs in
  
  observeEvent(input$logout_button, {
    shinyjs::reset("databasefetch")
    shinyjs::reset("genome")
    shinyjs::reset("group")
    shinyjs::reset("track")
    shinyjs::reset("Region1")
    shinyjs::reset("newsequence")
    shinyjs::reset("Sequence1")
    shinyjs::reset("flipstrand")
    shinyjs::reset("getFeature")
    shinyjs::reset("flank")
    shinyjs::reset("windowsize")
    shinyjs::reset("regionname")
    shinyjs::reset("startfeature")
    shinyjs::reset("startexon")
    shinyjs::reset("startexonboundary")
    shinyjs::reset("endfeature")
    shinyjs::reset("endexon")
    shinyjs::reset("endexonboundary")
    shinyjs::reset("startflank")
    shinyjs::reset("endflank")
    shinyjs::reset("startdirection")
    shinyjs::reset("enddirection")
    shinyjs::reset("clearregions")
    shinyjs::reset("smooth")
    shinyjs::reset("split")
    shinyjs::reset("tsvsplitting")
    shinyjs::reset("splitby")
    shinyjs::reset("heatmap_col_fun")
    shinyjs::reset("heatmapquantiles")
    shinyjs::reset("maxylim")
    shinyjs::reset("showrownames")
    shinyjs::reset("kmeansclustering")
    shinyjs::reset("logenriched")
    shinyjs::reset("row_km")
    shinyjs::reset("colorpalette")
    shinyjs::reset("averageprofile")
    shinyjs::reset("autoz")
    shinyjs::reset("log2")
    shinyjs::reset("dottedlines")
    shinyjs::reset("plottitle")
    shinyjs::reset("unit")
    shinyjs::reset("feature")
    shinyjs::reset("averageprofilequantiles")
    shinyjs::reset("alpha")
  })
  
  ## File upload logic
  
  # Region file upload text output
  
  output$regionfile1_name <- renderText({
    if (input$databasefetch) {
      "No files uploaded"
    } else if (!is.null(input$Region1)) {
      paste(input$Region1$name, collapse = ", ")
    } else {
      "No files uploaded"
    }
  })
  
  # Resetting text if user logs out
  
  observeEvent(input$logout_button, {
    output$regionfile1_name <- renderText({
      "No files uploaded"
    })
  })
  
  # Creating the file input for region files
  
  output$Region1_ui <- renderUI({
    if (input$databasefetch) {
      fileInput(
        inputId = "Region1", 
        label = "Upload region files (.bed/.gtf)", 
        accept = c(".bed", ".gtf"),
        multiple = TRUE)
    } else {
      fileInput(
        inputId = "Region1", 
        label = "Upload region files (.bed/.gtf)", 
        accept = c(".bed", ".gtf"),
        multiple = TRUE)
    }
  })
  
  ## Saving sequence data files
  
  # Creating the file input for sequence files
  
  output$Sequence1_ui <- renderUI({
    fileInput(
      inputId = "Sequence1",
      label = "Upload sequence data file (.bw)",
      accept = ".bw",
      multiple = FALSE
    )
  })
  
  # Reading in and calling user-specific pre-uploaded sequence files
  
  saved_sequence <- reactiveValues(list = list())
  
  observe({
    req(session$user$userData)
    sequence_dir <- file.path(session$userData$user_dir, "sequence data files")
    
    if(!dir.exists(sequence_dir)) {
      dir.create(sequence_dir, recursive = TRUE)
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
  
  # If new sequence file uploaded then saved to their directory as an RDS file
  
  observeEvent(input$savesequencedata, {
    req(input$Sequence1, session$userData$user_dir)
    
    sequence_dir <- get_sequence_dir()
    
    if (!dir.exists(sequence_dir)) {
      dir.create(sequence_dir, recursive = TRUE)
    }
    
    tryCatch({
      filepath <- input$Sequence1$datapath
      filename <- input$Sequence1$name
      
      # Add progress bar
      withProgress(message = 'Saving sequence data...', value = 0, {
        bigwig <- import.bw(filepath) # Import the file
        incProgress(0.5, detail = "Importing bigwig file") #update progress bar
        
        RDS <- saveRDS(object = bigwig, file = file.path(sequence_dir, paste0(filename, ".rds")))
        incProgress(0.75, detail = "Saving RDS file") #update progress bar
      })
      
      showNotification(paste0("Sequence data file ", filename, " saved successfully!"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving sequence data file", e$message), type = "error")
    })
  })
  
  # Updating the select input with the users sequnence files in their directory
  
  observe({
    sequence_selection <- saved_sequence_poll()
    updateSelectInput(session, "sequencedatafiles",
                      choices = sequence_selection,
                      selected = NULL)
  })
  
  # Clearing sequence files from the user's directory
  
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
  
  # Output of the sequence file name to text
  
  output$seqfile1_name <- renderText({
    if(!is.null(input$Sequence1)) {
      paste(input$Sequence1$name, collapse = ", ")
    } else {
      "No file uploaded"
    }
  })
  
  observeEvent(input$logout_button, {
    output$seqfile1_name <- renderText({
      "No file uploaded"
    })
  })
  
  
  ## Database annotation fetching
  
  # Selecting genome
  
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
  
  # Selecting group
  
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
  
  # Selecting track
  
  output$picktrack <- renderUI({
    if (input$databasefetch) {
      if (!is.null(input$genome)) { 
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
        return(NULL) 
      }
    } else {
      return(NULL) 
    }
  })
  
  # Adding an action button where the user then fetches the annotation file they have specified
  
  output$getannotation <- renderUI({
    if(input$databasefetch){
      actionButton(
        inputId = "fetchannotation",
        label = "Fetch annotation file"
      )
    }
  })
  
  ## Code to fetch the data once fetch annotation has been clicked
  
  gr <- reactiveVal(GRanges())
  
  observeEvent(input$fetchannotation, {
    showNotification("Fetching annotation...", type = "message")
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
  
  # Resetting the annotation name output
  
  observeEvent(input$logout_button, {
    output$annotationname <- renderPrint({
      NULL
    })
  })
  
  region_files_count <- reactive({
    req(input$Region1)
    length(input$Region1$datapath)
  })
  
  
  
  
  # Creating reactive objects for the matrix customisation options
  
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
  
  observeEvent(input$startflank, {
    if (is.na(input$startflank)) {
      updateNumericInput(session, "startflank", value = 0)
    }
  })
  
  observeEvent(input$endflank, {
    if (is.na(input$endflank)) {
      updateNumericInput(session, "endflank", value = 0)
    }
  })
  
  # Creating the getFeature region list
  
  saved_regions <- reactiveVal(list())
  region_objects_list <- reactiveVal(list())
  wins_vector <- reactiveVal(c())
  
  observeEvent(input$getregion, {
    if(is.null(input$Region1) && input$databasefetch == FALSE){
      showNotification("No region data detected.", type = "error")
    }
  })
  
  # Logic for once a region has been specified and the getregion action button is clicked
  
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
    
    tryCatch({
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
      current_wins[region_name] <- as.numeric(window_size) 
      wins_vector(current_wins)
      
      current_regions <- saved_regions()
      saved_regions(c(current_regions, list(new_region)))
      
    }, error = function(e) {
      showNotification(
        "There are no genes with that feature. Please enter another feature.",
        type = "error"
      )
    })
  })
  
  # Clearing the region list 
  
  observeEvent(input$clearregions, {
    saved_regions(list())
    region_objects_list(list())
    wins_vector(c())
  })
  
  # Creating a data table to store the regions created 
  
  output$savedRegionsTable <- renderDT({
    if (input$getFeature == 4) {
      regions_list <- saved_regions()
      if (length(regions_list) > 0) {
        df <- do.call(rbind, regions_list)
        datatable(df)
      } else {
        datatable(data.frame(Message = "No regions saved yet."))
      }
    }
  })
  
  # Resetting the datatable upon logout
  
  observeEvent(input$logout_button, {
    output$savedRegionsTable <- renderDT({
      datatable(data.frame(Message = "No regions saved yet."))
    })
  })
  
  
  
  # Matrix list generation from all the inputs
  
  
  matl <- eventReactive(input$matrixgeneration, {
    req(input$sequencedatafiles)
    
    withProgress(message = "Matrix Generation", value = 0, {
      
      regions <- list()
      
      # Fetch the sequence directory and region files
      
      sequence_dir <- get_sequence_dir()
      
      region_file <- if (!is.null(input$Region1)) input$Region1$datapath else NULL
      region2_file <- if (!is.null(input$Region2)) input$Region2$datapath else NULL
      all_rds_files <- list.files(sequence_dir, pattern = "\\.rds$", full.names = TRUE)
      
      all_rds_file_names_no_ext <- tools::file_path_sans_ext(basename(all_rds_files))
      
      selected_rds_files <- all_rds_files[all_rds_file_names_no_ext %in% input$sequencedatafiles]
      
      bigwig_files <- lapply(selected_rds_files, readRDS)
      names(bigwig_files) <- tools::file_path_sans_ext(basename(selected_rds_files))
      bigwig_file_names <- names(bigwig_files)
      
      # Import region files
      
      if (!is.null(input$Region1)){
        b <- readBed(region_file)
        combined_b <- if (!is.null(region2_file)) c(b, readBed(region2_file)) else b
      } else if (is.null(input$Region1)) {
        combined_b <- gr()
      }
      
      
      incProgress(0.1, detail = "Reading region files...")
      
      # Getting features the region file
      
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
      
      # Filtering names of bigwig files and storing in a forward and reverse strand object
      
      if (all(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
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
        print("bwf:")
        print(bwf)
      } else {
        bwr <- NULL
        bwf <- bigwig_files[grepl("\\.bw$", bigwig_file_names)]
        if (length(bwf) > 0) {
          names(bwf) <- sub("\\.bw$", "", bigwig_file_names[grepl("\\.bw$", bigwig_file_names)])
        }
        print(bwf)
      }
      
      # Setting strandedness
      
      strand_reactive <- reactive({
        if(length(bwf) > 0 && length(bwr) > 0){
          if(input$flipstrand){
            return("rev")
          } else {
            return("for")
          }
        } else {
          return("no")
        }
      })
      print(strand_reactive())
      
      incProgress(0.5, detail = "Generating matrices")
      
      # Deploying the matList function and generating the matrices
      
      if (any(grepl("\\.f\\.bw$|\\.r\\.bw$", bigwig_file_names))) {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(bwf), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(bwf), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(bwf), wins = wins_vector(), strand = strand_reactive(), smooth = smooth_reactive())
        }
      } else {
        if (input$getFeature == 1) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(bwf), extend = input$flank, w = input$windowsize, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 2 || input$getFeature == 3) {
          grl <- list("features" = features)
          matl_result <- matList(bwf = bwf, grl = grl, names = names(bwf), extend = input$flank, w = 1, strand = strand_reactive(), smooth = smooth_reactive())
        } else if (input$getFeature == 4) {
          matl_result <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(bwf), wins = wins_vector(), strand = strand_reactive(), smooth = smooth_reactive())
        } 
      }
      
      incProgress(0.9, detail = "Almost finished...")
      return(matl_result)
    })
  })
  
  # Rendering the names of the matrices
  
  output$matrixnames <- renderText({
    req(matl())
    mat_names <- names(matl())
    if(is.null(mat_names) || length(mat_names) == 0) {
      return("No matrices generated")
    }
    paste(mat_names, collapse = ", ")
  })
  
  observeEvent(input$logout_button, {
    output$matrixnames <- renderText({
      "No matrices generated"
    })
  })
  
  ## Saving Matrices
  
  # Matrices saved as .rds files

  saved_matrices <- reactiveValues(list = list())
  
  observe({
    req(session$userData$user_dir) 
    save_dir <- file.path(session$userData$user_dir, "matrices")
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      print(paste("Created directory:", save_dir))
    }
  })
  
  # Function to get save_dir
  get_save_dir <- function() {
    if (!is.null(session$userData$user_dir)) {
      return(file.path(session$userData$user_dir, "matrices"))
    } else {
      return(NULL)
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
        return(Sys.time())
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
  
  # Event to handle saving matrices once action button clicked
  
  observeEvent(input$savematrices, {
    req(matl(), session$userData$user_dir)
    
    save_dir <- get_save_dir()
    
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      print(paste("Created directory:", save_dir))
    }
    
    mat_list <- matl()
    
    for (mat_name in names(mat_list)) {
      file_path <- file.path(save_dir, paste0(mat_name, ".rds"))
      tryCatch({
        saveRDS(mat_list[[mat_name]], file = file_path)
        Sys.sleep(0.1)
        showNotification(paste("Matrix", mat_name, "saved successfully!"), type = "message")
        print(paste("Saved matrix:", mat_name, "to", file_path))
      }, error = function(e) {
        showNotification(paste("Error saving matrix", mat_name, ":", e$message), type = "error")
        print(paste("Error saving matrix:", mat_name, e$message))  
      })
    }
  })
  
  # Display saved matrices as text
  
  output$savedmatrices <- renderText({
    matrix_names <- saved_matrices_poll()
    if (length(matrix_names) == 0) {
      return("No matrices saved")
    } else {
      return(paste(matrix_names, collapse = ", "))
    }
  })
  
  
  output$renamematrixselect <- renderUI({
    matrix_names <- saved_matrices_poll()
    selectInput(
      inputId = "selectmatrixtorename",
      label = "Select matrix",
      choices = matrix_names,
      selected = NULL
    )
  })
  
  observeEvent(input$rename, {
    req(input$selectmatrixtorename, req(input$newmatrixname))
    if(is.null(input$newmatrixname)){
      showNotification("Please enter a new name", type = "warning")
    } else {
      save_dir <- get_save_dir()
    
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        print(paste("Created directory:", save_dir))
      }
    
      clean_name <- gsub("[^a-zA-Z0-9_-]", "_", input$newmatrixname)
    
      new_filename <- paste0(clean_name, ".rds")
    
      old_path <- file.path(save_dir, paste0(input$selectmatrixtorename, ".rds"))
      new_path <- file.path(save_dir, new_filename)
    
      if(file.exists(new_path)) {
        showNotification("A file with this name already exists", type = "error")
      } else {
        success <- file.rename(from = old_path, to = new_path)
        
        if(success) {
          showNotification("Matrix renamed successfully", type = "message")
        } else {
          showNotification("Failed to rename matrix", type = "error")
        }
      }
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
            print(paste("File deleted:", file_path))  
          }, error = function(e) {
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
  
  ## Splitting file uploads
  
  # Creating the tsv upload for a splitting file in the visualisations stage
  
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
  
  # Select input for which column in the .tsv file they are splitting by
  
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
  
  # Resetting the splitting choices upon logout
  
  observeEvent(input$logout_button, {
    column_choices(NULL)
  })
  
  # Updating column_choices with the columns in the splitting file except the first one (gene name)
  
  observeEvent(input$tsvsplitting, {
    if (!is.null(input$tsvsplitting)) {
      annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if (ncol(annotation) > 1) {
        column_choices(colnames(annotation)[-1])
      } else {
        column_choices(NULL)
      }
    } else {
      column_choices(NULL)
    }
  })
  
  
  
  
  
  ### enrichedHeatmap Heatmaps
  
  # Creating reactive objects for the heatmaps customisation options
  
  # Creating the conditional splitting object
  split_reactive <- reactive({
    req(input$split, input$tsvsplitting$datapath, input$splitby, selected_matrices_reactive()[1])
    
    annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(annotation)[1] <- "gene_name"
    
    split_col <- input$splitby
    if (!split_col %in% colnames(annotation)) {
      stop("Selected column not found in annotation file.")
    }
    
    Anno <- data.frame(name = rownames(selected_matrices_reactive()[[1]])) %>%
      left_join(data.frame(name = annotation$gene_name, split_val = factor(annotation[[split_col]])), by = "name") %>%
      column_to_rownames("name")
    
    return(Anno)
  })
  
  
  split_cols_reactive <- reactive({
    req(split_reactive())
    unique_families <- unique(split_reactive()$split_val)
    
    if (length(unique_families) == 2) {
      colors <- c("#DA4167", "#083D77")
      names(colors) <- unique_families[1:2]
      list(split_val = colors)
    } else if (length(unique_families) > 2) {
      color_palette <- grDevices::rainbow(length(unique_families))
      names(color_palette) <- unique_families
      names(color_palette) <- paste0("`", names(color_palette), "`") 
      list(split_val = color_palette)
    } else {
      list(split_val = c("default" = "gray"))
    }
  })
  
  # Adding ability to filter heatmaps rather than just splitting
  
  filter_choices <- reactiveVal(NULL)
  
  observeEvent(c(input$tsvsplitting, input$splitby), {
    req(input$split, input$tsvsplitting, input$splitby)
    if(!is.null(input$tsvsplitting) && !is.null(input$splitby)) {
      annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      split_col <- input$splitby
      split_val <- factor(annotation[[split_col]])
      unique_choices <- unique(split_val)
      filter_choices(c("All", as.character(unique_choices)))
    } else {
      filter_choices(NULL)
    }
  })
  
  output$filterselect <- renderUI({
    if(input$split){
      selectInput(
        inputId = "filterby",
        label = "Filter by:",
        choices = filter_choices(),
        selected = "All",
        multiple = FALSE
      )
    }
  })
  
  filtered_data <- reactive({
    req(input$filterby, split_reactive(), selected_matrices_reactive())
    
    if (input$filterby == "All") {
      return(selected_matrices_reactive())
    } else {
      selected_rows <- rownames(split_reactive())[split_reactive()$split_val == input$filterby]
      return(selected_matrices_reactive()[selected_rows, , drop = FALSE]) 
    }
  })
  
  ## Creating other reactive objects for heatmap customisation
  
  # Reading in the matrices selected by the user from their personal directory
  
  selected_matrices_reactive <- reactive({
    req(input$selectedmatrices, session$userData$user_dir)  
    
    save_dir <- file.path(session$userData$user_dir, "matrices")
    
    selected_mats <- lapply(input$selectedmatrices, function(mat_name) {
      mat_path <- file.path(save_dir, paste0(mat_name, ".rds"))  
      if (file.exists(mat_path)) {
        return(readRDS(mat_path))
      } else {
        return(NULL)
      }
    })
    
    selected_mats <- selected_mats[!sapply(selected_mats, is.null)]  
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
  
  # Windows depend on the features they are looking at - whether it is the full gene or around a given point
  
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
  
  ## Creating reactive objects for the splitting, matrix selection, and split colours for when filtering
  
  # Editing matrices and splitting object and colours so that it is only what the user has selected to filter by
  
  filtered_matrices_reactive <- reactive({
    req(input$filterby, selected_matrices_reactive(), split_reactive())
    
    if (input$filterby == "All") {
      return(selected_matrices_reactive())
    } else {
      selected_rows <- rownames(split_reactive())[split_reactive()$split_val == input$filterby]
      filtered_mats <- lapply(selected_matrices_reactive(), function(mat) {
        mat[selected_rows, , drop = FALSE]
      })
      return(filtered_mats)
    }
  })
  
  filtered_split_reactive <- reactive({
    req(input$filterby, split_reactive())
    
    if (input$filterby == "All") {
      return(split_reactive())
    } else {
      filtered_split <- split_reactive()[split_reactive()$split_val == input$filterby, , drop = FALSE]
      return(filtered_split)
    }
  })
  
  filtered_split_cols_reactive <- reactive({
    req(input$filterby, split_reactive())
    
    if (input$filterby == "All") {
      return(split_cols_reactive())
    } else {
      filtered_colors <- split_cols_reactive()
      filtered_colors$split_val <- filtered_colors$split_val[names(filtered_colors$split_val) == input$filterby]
      return(filtered_colors)
    }
  })
  
  
  # Creating the heatmaps
  
  observeEvent(input$heatmapplotbutton, {
    if(input$split == TRUE && is.null(input$tsvsplitting)){
      showNotification("No splitting file detected", type = "error")
    }
  })
  
  observeEvent(input$heatmapplotbutton, {
    if(length(input$selectedmatrices) == 0){
      showNotification("Matrices need to be selected before plotting", type = "warning")
    } 
  })
  
  hml <- eventReactive(input$heatmapplotbutton, {
    req(selected_matrices_reactive())
    
    tryCatch({
      withProgress(message = 'Creating Heatmap...', value = 0.4, {
        
        if (isTRUE(input$split)) {
          hml <- hmList(
            matl = filtered_matrices_reactive(),
            wins = wins_reactive(),
            split = filtered_split_reactive(),
            split_cols = filtered_split_cols_reactive(),
            col_fun = heatmap_col_fun_reactive(),
            axis_labels = axis_labels_reactive(),
            show_row_names = show_row_names_reactive(),
            min_quantile = heatmap_min_quantile_reactive(),
            max_quantile = heatmap_max_quantile_reactive(),
            row_km = input$row_km,
            ylim = c(0, max_ylim_reactive()),
            log2 = input$logenriched
          )
          incProgress(1, detail = "Heatmap created") 
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
            row_km = input$row_km,
            log2 = input$logenriched
          )
          incProgress(1, detail = "Heatmap created") 
        }
      })
      
      return(hml)
      
    }, error = function(e) {
      showNotification("Plotting error", type = "error")
      return(NULL)
    })
  })
  
  observeEvent(input$logout_button, {
    output$enrichedHeatmapPlot <- renderPlot({
      NULL
    })
  })
  
  output$enrichedHeatmapPlot <- renderPlot({
    req(hml())
    req(length(hml()) > 0)
    
    withProgress(message = 'Drawing Heatmap...', value = 0, {
      combined_hm <- Reduce(`+`, hml())
      incProgress(0.5, detail = "Combining heatmaps")
      draw(combined_hm, merge_legend = TRUE)
      incProgress(0.7, detail = "Drawing plot")
    })
  })
  
  
  ## Download options for the enrichedHeatmap heatmaps as a PDF or PNG
  
  output$heatmapdownloadpng <- downloadHandler(
    filename = function() { paste("enrichedHeatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)
      req(length(hml()) > 0)
      combined_hm <- Reduce(`+`, hml())
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  
    }
  )
  
  output$heatmapdownloadpdf <- downloadHandler(
    filename = function() { paste("enrichedHeatmap_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)
      req(length(hml()) > 0)
      combined_hm <- Reduce(`+`, hml())
      draw(combined_hm, merge_legend = TRUE)
      dev.off()  
    }
  )
  
  ## ggplot heatmap
  
  # Creating splitting and filtering objects and rendering the UI
  
  output$ggfilterselect <- renderUI({
    if(input$split){
      selectInput(
        inputId = "ggfilterby",
        label = "Filter by:",
        choices = filter_choices(),
        selected = "All",
        multiple = FALSE
      )
    }
  })
  
  ggsplit_reactive <- reactive({
    if (!is.null(input$tsvsplitting)) {
      tryCatch({
        annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        if (ncol(annotation) > 0) {
          colnames(annotation)[1] <- "gene_name" # Ensure the first column is named 'gene_name'
          
          split_col <- input$splitby
          if (!split_col %in% colnames(annotation)) {
            showNotification(paste("Selected column '", split_col, "' not found in annotation file."), type = "error")
            return(NULL)
          }
          
          Anno <- data.frame(name = rownames(selected_matrices_reactive()[[1]])) %>%
            left_join(data.frame(name = annotation$gene_name, split_val = factor(annotation[[split_col]])), by = "name") %>%
            column_to_rownames("name")
          return(Anno)
        } else {
          showNotification("Annotation file is empty or has no columns.", type = "warning")
          return(NULL)
        }
      }, error = function(e) {
        showNotification(paste("Error reading or processing annotation file:", e$message), type = "error")
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })
  
  ggplot_filtered_matrices_reactive <- reactive({
    req(input$ggfilterby, selected_matrices_reactive(), ggsplit_reactive())
    
    if (input$ggfilterby == "All") {
      return(selected_matrices_reactive())
    } else {
      selected_rows <- rownames(ggsplit_reactive())[ggsplit_reactive()$split_val == input$ggfilterby]
      ggfiltered_mats <- lapply(selected_matrices_reactive(), function(mat) {
        mat[selected_rows, , drop = FALSE]
      })
      return(ggfiltered_mats)
    }
  })
  
  ggwins_reactive <- reactive({
    if (input$getFeature == 1) {
      flank_size <- input$flank
      ggwins <- c("Upstream" = flank_size, "Gene" = 2 * flank_size, "Downstream" = flank_size)
      return(ggwins)
    } else if (input$getFeature == 2){
      flank_size <- input$flank
      ggwins <- c("Upstream" = flank_size, "TSS" = 1, "Downstream" = flank_size)
      return(ggwins)
    } else if (input$getFeature == 3){
      flank_size <- input$flank
      ggwins <- c("Upstream" = flank_size, "TES" = 1, "Downstream" = flank_size)
      return(ggwins)
    } else if (input$getFeature == 4) {
      return(wins_vector())
    } else {
      return(NULL)
    }
  })
  
  ggbreaklabels_reactive <- reactive({
    if (input$getFeature == 1) {
      flank_size <- input$flank
      ggbreaklabels <- c(paste0("-", flank_size, "b"), "TSS", "TES", paste0("+", flank_size, "b"))
      return(ggbreaklabels)
    } else if (input$getFeature == 2){
      flank_size <- input$flank
      ggbreaklabels <- c(paste0("-", flank_size, "b"), "TSS", "", paste0("+", flank_size, "b"))
      return(ggbreaklabels)
    } else if (input$getFeature == 3) {
      flank_size <- input$flank
      ggbreaklabels <- c(paste0("-", flank_size, "b"), "TES", "", paste0("+", flank_size, "b"))
      return(ggbreaklabels)
    } else if (input$getFeature == 4) {
      wins <- ggwins_reactive()
      if (length(wins) > 0) {
        labels <- names(wins)
        return(c(labels, ""))
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  
  
  
  zscaling_reactive <- reactive({
    if(input$autoz == TRUE){
      z <- "auto"
    } else {
      z <- NULL
    }
    return(z)
  })
  
  # Creating the ggplot heatmap object and rendering the plot
  
  observeEvent(input$plotggheatmap, {
    if(length(input$selectedmatrices) == 0){
      showNotification("Matrices need to be selected before plotting", type = "warning")
    } 
  })
  
  observeEvent(input$plotggheatmap, {
    if(input$split == TRUE && is.null(input$tsvsplitting)){
      showNotification("No splitting file detected", type = "error")
    }
  })
  
  ggplot_heatmap_object <- eventReactive(input$plotggheatmap, {
    req(length(selected_matrices_reactive()) > 0)
    tryCatch({
      withProgress(message = 'Creating Heatmap...', value = 0.5, {
        if (isTRUE(input$split)) {
          plot <- plotggplotHeatmap(
            matl = ggplot_filtered_matrices_reactive(),
            wins = ggwins_reactive(),
            break_labels = ggbreaklabels_reactive(),
            color_palette = input$colorpalette,
            average_profile = input$averageprofile,
            zMin = zscaling_reactive(),
            zMax = zscaling_reactive(),
            log2 = input$log2,
            dottedlines = input$dottedlines,
            split = ggsplit_reactive()
          )
          incProgress(1, detail = "Heatmap created")
          print(plot)
          plot
        } else {
          plot <- plotggplotHeatmap(
            matl = selected_matrices_reactive(),
            wins = ggwins_reactive(),
            break_labels = ggbreaklabels_reactive(),
            color_palette = input$colorpalette,
            average_profile = input$averageprofile,
            zMin = zscaling_reactive(),
            zMax = zscaling_reactive(),
            log2 = input$log2,
            dottedlines = input$dottedlines
          )
          incProgress(1, detail = "Heatmap created")
          print(plot)
          plot
        }
      })
      return(plot)
    }, error = function(e) {
      showNotification("Plotting error", type = "error")
      return(NULL) # Return NULL to prevent further errors
    })
  })
  
  observeEvent(input$logout_button, {
    output$ggplotheatmap <- renderPlot({
      NULL
    })
  })
  
  output$ggplotheatmap <- renderPlot({
    withProgress(message = 'Drawing Heatmap...', value = 0.5, {
      ggplot_heatmap_object()
        incProgress(1, detail = "Heatmap rendered")
    })
  }, height = 1000)
  
  # Download options for the ggplot heatmap as either a PDF or PNG
  
  output$ggheatmapdownloadpng <- downloadHandler(
    filename = function() { paste("ggplot2heatmap_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)
      print(ggplot_heatmap_object())
      dev.off()
    }
  )
  
  output$ggheatmapdownloadpdf <- downloadHandler(
    filename = function() { paste("ggplot2heatmap_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)
      print(ggplot_heatmap_object())
      dev.off()
    }
  )
  
  
  ## Average profile plot
  # Customisation options for the average profile plot
  
  averagecolourbychoices <- reactiveVal(NULL)
  averagefacetbychoices <- reactiveVal(NULL)
  
  output$columnstoinclude <- renderUI({
    if(input$split){
      selectInput(
        inputId = "columninclude",
        label = "Select splitting columns to include in plot",
        multiple = TRUE,
        choices = column_choices(),
        selected = NULL
      )
    }
  })
  
  # Update averagecolourbychoices when input$columninclude changes
  observeEvent(input$columninclude, {
    if (!is.null(input$columninclude)) {
      averagecolourbychoices(c("Sample", input$columninclude))
      averagefacetbychoices(c("Sample", input$columninclude))
    } else {
      averagecolourbychoices(c("Sample"))
      averagefacetbychoices(c("Sample"))
    }
  })
  
  observeEvent(input$tsvsplitting, {
    if(!is.null(input$tsvsplitting)){
      annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if(ncol(annotation) > 1) {
        averagecolourbychoices(c("Sample", input$columninclude))
        averagefacetbychoices(c("Sample", input$columninclude))
      } else {
        averagefacetbychoices(NULL)
        averagecolourbychoices(NULL)
      }
    } else {
      averagefacetbychoices(NULL)
      averagecolourbychoices(NULL)
    }
  })
  
  output$colouringby <- renderUI({
    if(input$split){
      selectInput(
        inputId = "colourby",
        label = "Colour by:",
        multiple = FALSE,
        choices = averagecolourbychoices(),
        selected = averagecolourbychoices()[1]
      )
    } else {
      NULL
    }
  })
  
  output$facetingby <- renderUI({
    if(input$split){
      selectInput(
        inputId = "facetby",
        label = "Facet by:",
        multiple = TRUE,
        choices = averagefacetbychoices(),
        selected = NULL
      )
    }
  })
  
  average_split_reactive <- reactive({
    req(input$split, input$tsvsplitting$datapath, input$columninclude)
    annotation <- read.table(input$tsvsplitting$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(annotation)[1] <- "gene_name"
    
    selected_cols <- input$columninclude
    if(!all(selected_cols %in% colnames(annotation))) {
      showNotification("One or more selected columns not found in annotation file", type = "error")
      stop("One or more selected columns not found in annotation file")
    }
    
    Anno <- annotation |>
      select(gene_name, all_of(selected_cols))
    names(Anno)[1] = "name"
    
    Anno <- Anno |> column_to_rownames("name")
    
    return(Anno)
    
  })
  
  
  average_split_cols_reactive <- reactive({
    req(average_split_reactive(), input$colourby)  # Ensure both inputs are available
    
    # Get the specified column name for color by
    colourby_col <- input$colourby
    
    # Check if the selected column exists in the reactive dataframe
    if (!colourby_col %in% c("Sample",colnames(average_split_reactive()))) {
      stop("Selected column not found in the reactive data.")
    }
    
    if(colourby_col == "Sample"){
      color_palette <- grDevices::rainbow(length(names(selected_matrices_reactive())))
      names(color_palette) <- names(selected_matrices_reactive())
      return(color_palette)
    } else {
    
    # Get the unique families from the specified column
    unique_families <- unique(average_split_reactive()[[colourby_col]])
    
    # If unique families are found, generate a color palette
    if (length(unique_families) > 0) {
      color_palette <- grDevices::rainbow(length(unique_families))
      names(color_palette) <- unique_families
      return(color_palette)
    } else {
      return(NULL)
    }
    }
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
  
  avg_breaks_reactive <- reactive({
    if(input$getFeature == 1){
      x_range <- ncol(selected_matrices_reactive()[[1]])
      break_positions <- c(0, 0.25, 0.75, 1) * x_range
      return(break_positions)
    }
    if(input$getFeature == 2){
      x_range <- ncol(selected_matrices_reactive()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
      return(break_positions)
    }
    if(input$getFeature == 3){
      x_range <- ncol(selected_matrices_reactive()[[1]])
      break_positions <- c(0, 0.5, 1) * x_range
      return(break_positions)
    }
    if (input$getFeature == 4) {
      wins <- wins_vector()
      x_range <- ncol(selected_matrices_reactive()[[1]])
      cumulative_sums <- cumsum(wins)
      break_proportions <- c(0, cumulative_sums / sum(wins))
      break_positions <- break_proportions * x_range
      return(break_positions)
    }
  })
  
  
  avg_breaklabels_reactive <- reactive({
    if (input$getFeature == 1) {
      flank_size <- input$flank
      avgbreaklabels <- c(paste0("-", flank_size, "b"), "TSS", "TES", paste0("+", flank_size, "b"))
      return(avgbreaklabels)
    } else if (input$getFeature == 2){
      flank_size <- input$flank
      avgbreaklabels <- c(paste0("-", flank_size, "b"), "TSS", paste0("+", flank_size, "b"))
      return(avgbreaklabels)
    } else if (input$getFeature == 3) {
      flank_size <- input$flank
      avgbreaklabels <- c(paste0("-", flank_size, "b"), "TES", paste0("+", flank_size, "b"))
      return(avgbreaklabels)
    } else if (input$getFeature == 4) {
      wins <- ggwins_reactive()
      flank_size <- input$flank
      if (length(wins) > 0) {
        cumulative_sums <- cumsum(wins)
        total_sum <- sum(wins)
        break_proportions <- c(0, cumulative_sums / total_sum)
        x_range <- ncol(selected_matrices_reactive()[[1]])
        break_positions <- break_proportions * x_range
        labels <- c(names(wins), "")
        return(labels)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  
  unit_reactive <- reactive({
    input$unit
  })
  
  feature_reactive <- reactive({
    input$feature
  })
  
  # Defining a colour palette for the average profile plot
  
  pal = c(RColorBrewer::brewer.pal(n = 9,name = "Set1"))
  pal2 = c(wes_palette("Darjeeling2")[2],wes_palette("Zissou1")[1],wes_palette("Darjeeling1")[4],wes_palette("Darjeeling1")[3])
  colmap = c(pal[c(1:5,7:9)])
  
  
  
  # Creating the average profile plot
  
  observeEvent(input$averageprofileplotbutton, {
    if(length(input$selectedmatrices) == 0){
      showNotification("Matrices need to be selected before plotting", type = "warning")
    } 
  })
  
  observeEvent(input$logout_button, {
    output$averageprofileplot <- renderPlot({
      NULL
    })
  })
  
  output$averageprofileplot <- renderPlot({
    req(average_profile())
    average_profile()
  })
  
  observeEvent(input$averageprofileplotbutton, {
    if(input$split == TRUE && is.null(input$tsvsplitting)){
      showNotification("No splitting file detected", type = "error")
    }
  })
  
  average_profile <- eventReactive(input$averageprofileplotbutton, {
    req(selected_matrices_reactive())
    
    tryCatch({
      withProgress("Creating average profile...", value = 0.5, {
        if(input$split){
          if(is.null(input$columninclude)){
            showNotification("Select columns from splitting file to include", type = "warning")
          } else if(input$facetchoice == T && is.null(input$facetby)){
            showNotification("No faceting selected", type = "warning")
          } else {
        average_profile <- mplot(
          matl = selected_matrices_reactive(),
          split = average_split_reactive(),
          colmap = average_split_cols_reactive(),
          feature = feature_reactive(),
          unit = unit_reactive(),
          title = title_reactive(),
          min_quantile = averageprofile_min_quantile_reactive(),
          max_quantile = averageprofile_max_quantile_reactive(),
          colour_by = input$colourby,
          facet = input$facetby,
          alpha = alpha_reactive(),
          breaks = avg_breaks_reactive(),
          labels = avg_breaklabels_reactive(),
          facet_scale = "free_y"
        )
          }
        } else {
          average_profile <- mplot(
            matl = selected_matrices_reactive(),
            colmap = colmap,
            feature = feature_reactive(),
            unit = unit_reactive(),
            title = title_reactive(),
            min_quantile = averageprofile_min_quantile_reactive(),
            max_quantile = averageprofile_max_quantile_reactive(),
            alpha = alpha_reactive(),
            breaks = avg_breaks_reactive(),
            labels = avg_breaklabels_reactive()
          )
        }
        incProgress(1, "Average profile created")
      })
      return(average_profile)
    }, error = function(e) {
      showNotification("Plotting error", type = "error")
      return(NULL)
    })
  })
  
  ## Downloading average profile plots
  
  output$averageprofiledownloadpng <- downloadHandler(
    filename = function() { paste("averageprofileplot_", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 1350, height = 900, res = 150)  
      print(average_profile())
      dev.off() 
    }
  )
  
  output$averageprofiledownloadpdf <- downloadHandler(
    filename = function() { paste("averageprofileplot_", Sys.Date(), ".pdf", sep="") },
    content = function(file) {
      pdf(file, width = 13.5, height = 9)  
      print(average_profile())
      dev.off()  
    }
  )
  
}