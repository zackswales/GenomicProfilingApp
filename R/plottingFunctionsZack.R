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
#library(UCSC.utils)

## NOTES
## default to no strand?
## Add keep option? for normalizedMatrix
## Different mean options
## output DF or matrix or ggplot
## rewrite getFeature from rootGenomics

#colmap = RColorBrewer::brewer.pal(n = 4,name = "Set1")
colmap = c(wes_palette("Darjeeling2")[2],wes_palette("Zissou1")[1],wes_palette("Darjeeling1")[4],wes_palette("Darjeeling1")[3])

# future options
options(future.globals.maxSize= 30000000000)
plan(multisession, workers = 20)

importBWlist <- function(bwf,names,selection=NULL){
  if(is.null(selection)){
    bwl <- bwf |> future_map(~import(.x, format = "BigWig"))
  }
  else{
    bwl <- bwf |> future_map(~import(.x, format = "BigWig", selection = BigWigSelection(selection)))
  }
  names(bwl) <- names
  bwl
}

## Matrix
matList <- function(bwf, bwr, names, grl, wins = list("Gene" = 10), mode = "coverage", output = "norm.matrix", 
                    strand = "rev", smooth = FALSE, extend = 0, w, include_target = TRUE, target_ratio = 0.5, 
                    k = 10, keep) {
  
  fbw <- bwf
  #fbw <- bwf |> future_map(~import(.x, format = "BigWig"))
  #names(fbw) <- names
  
  if(strand %in% c("for", "rev")) {
    rbw <- bwr
    #rbw <- bwr |> future_map(~import(.x, format = "BigWig"))
    #names(rbw) <- names
    
    fgrl = grl |> map(~subset(.x, strand(.x) == "+"))
    rgrl = grl |> map(~subset(.x, strand(.x) == "-"))
  }
  
  if(length(wins) > 1) { ## Length of features to combine
    matl <- names |> future_map(function(s) {
      
      if(strand == "rev") {
        fmat <- fgrl |> imap(~normalizeToMatrix(rbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y], 
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |> 
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        rmat <- rgrl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y], 
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |> 
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        
        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "for") {
        fmat <- fgrl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y], 
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |> 
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        rmat <- rgrl |> imap(~normalizeToMatrix(rbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y], 
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |> 
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        
        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "no") {
        mat <- grl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y], 
                                              mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |> 
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)   
        rownames(mat) <- grl[[1]]$name
      }
      
      mat
      
    }, .progress = TRUE)
  }
  else {  ## Single feature - regular normalizeToMatrix
    matl <- names |> future_map(function(s) {
      
      if(strand == "rev") {
        fmat <- normalizeToMatrix(rbw[[s]], fgrl[[1]], extend = extend, value_column = "score", k = wins[1], 
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target, 
                                  target_ratio = target_ratio)
        rmat <- normalizeToMatrix(fbw[[s]], rgrl[[1]], extend = extend, value_column = "score", k = wins[1], 
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target, 
                                  target_ratio = target_ratio)
        
        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "for") {
        fmat <- normalizeToMatrix(fbw[[s]], fgrl[[1]], extend = extend, value_column = "score", k = wins[.y], 
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target, 
                                  target_ratio = target_ratio)
        rmat <- normalizeToMatrix(rbw[[s]], rgrl[[1]], extend = extend, value_column = "score", k = wins[1], 
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target, 
                                  target_ratio = target_ratio)
        
        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "no") {
        mat <- normalizeToMatrix(fbw[[s]], grl[[1]], extend = extend, value_column = "score", k = wins[.y], 
                                 mean_mode = mode, smooth = smooth, w = w, include_target = include_target, 
                                 target_ratio = target_ratio)
        rownames(mat) <- grl[[1]]$name
      }
      
      mat
      
    }, .progress = TRUE)
  }
  
  names(matl) <- names
  
  if(output == "norm.matrix") {
    matl
  }
  else if(output == "matrix") {
    matl |> map(~as.data.frame(.x) |> as.matrix())
  }
}


# encrichedHeatmap list function
hmList <- function(matl, wins, split = NULL, split_cols, max_quantile = 0.99, min_quantile = 0, col_fun = "red", show_row_names = TRUE, win_labels = NULL, ylim = NULL, summarise_by = "mean", axis_labels = "", row_km = 1, log2 = FALSE) {
  
  ## ColourMap
  reds <- RColorBrewer::brewer.pal(n = 9, name = "Reds")
  
  if (log2) {
    matl <- lapply(matl, function(x) log2(x + 1))
  }
  
  common_min <- quantile(unlist(matl), min_quantile)
  common_max <- quantile(unlist(matl), max_quantile)
  if (col_fun == "red") {
    col_fun <- circlize::colorRamp2(c(common_min, common_max), c("white", "red"))
  } else if (col_fun == "bl2rd") {
    col_fun <- circlize::colorRamp2(c(common_min, 0, common_max), c("blue", "white", "red"))
  } else if (col_fun == "red0") {
    col_fun <- circlize::colorRamp2(c(0, 1, common_max), c("white", reds[3], reds[7]))
  }
  
  if (is.null(ylim)) {
    ymin <- min(0, common_min)
    ymax <- common_max
    offset <- quantile(1:(ymax - ymin), 0.01)
    ylim <- c(ymin - offset, ymax + offset)
  }
  
  features <- vector()
  fcols <- RColorBrewer::brewer.pal(9, "Set1")[1:length(wins)]
  names(fcols) <- names(wins)
  for (i in 1:length(wins)) {
    features <- append(features, rep(names(wins)[i], wins[i]))
  }
  
  if (is.null(win_labels)) {
    win_labels <- names(wins)
  }
  
  if (!is.null(split)) {
    
    rowAnno <- rowAnnotation(df = split, col = split_cols, show_legend = FALSE)
    
    hml <- matl |> imap(~EnrichedHeatmap(
      .x,
      name = .y,
      column_title = .y,
      col = col_fun,
      show_row_names = show_row_names,
      row_names_gp = gpar(fontsize = 5),
      axis_name = axis_labels,
      row_split = split,
      row_km = row_km,
      column_title_gp = gpar(fontsize = 10),
      cluster_rows = FALSE,
      top_annotation = c(
        HeatmapAnnotation(
          Features = features,
          show_legend = TRUE,
          border = TRUE,
          col = list(Features = fcols),
          show_annotation_name = FALSE,
          annotation_legend_param = list(at = names(wins), labels = win_labels)
        ),
        HeatmapAnnotation(
          enriched = anno_enriched(
            value = summarise_by,
            gp = gpar(fontsize = 4, lwd = 2, col = split_cols[[1]]),
            axis_param = list(side = "left", facing = "inside"),
            ylim = ylim
          )
        ),
        gap = unit(2, "mm")
      )
    ))
    hml$rowAnno <- rowAnno
    hml
  } else {
    hml <- matl |> imap(~EnrichedHeatmap(
      .x,
      name = .y,
      column_title = .y,
      col = col_fun,
      show_row_names = show_row_names,
      row_names_gp = gpar(fontsize = 5),
      axis_name = axis_labels,
      column_title_gp = gpar(fontsize = 10),
      row_km = row_km,
      cluster_rows = FALSE,
      top_annotation = c(
        HeatmapAnnotation(
          Features = features,
          show_legend = TRUE,
          border = TRUE,
          col = list(Features = fcols),
          show_annotation_name = FALSE,
          annotation_legend_param = list(at = names(wins), labels = win_labels)
        ),
        HeatmapAnnotation(
          enriched = anno_enriched(
            gp = gpar(fontsize = 4, lwd = 2),
            axis_param = list(side = "left", facing = "inside"),
            ylim = ylim
          )
        ),
        gap = unit(2, "mm")
      )
    ))
    hml
  }
}

## Metaplot
mplot<-function(matl,matlc=NULL,feature = "Gene", unit = "Coverage (BPM)", title = "Gene metaplot", breaks = c(0,20,60,80), labels = c("-200","TSS","TTS","+200"), colmap, split=NULL, facet=NULL, angle = 0, strip_fill = "white",facet_scale="fixed", max_quantile = 1, min_quantile = 0, pseudo = 1,alpha=1, linewidth=0.5, error = F, alpha_error=0.5,facet_nrow=2,facet_type="wrap",facet_independent=F,summarise_by="mean",colour_by="Sample"){
  
  ml <- matl
  if(!is.null(matlc)){
    ml <- map2(matl,matlc,~log2((.x + pseudo) / (.y + pseudo)))  
  }
  
  df <- ml |> imap(~as.data.frame(.x) |>
                     mutate(Sample = .y) |>
                     rownames_to_column("name")) |>
    bind_rows()
  
  if(!is.null(split)){
    if(error == T & summarise_by == "mean"){
      df.error <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~std.error(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
    if(summarise_by == "mean"){
      df <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~mean(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
    if(summarise_by == "sum"){
      df <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~sum(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
  }
  else{
    if(error == T & summarise_by == "mean"){
      df.error <- df |> group_by(Sample) |>
        summarise(across(-name,~std.error(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T)))
    }
    df <- df |> group_by(Sample) |>
      summarise(across(-name,~mean(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T)))
  }
  
  names <- c("Sample")
  if(!is.null(split)){
    names <- c(names,names(split))
  }
  
  names(df) <- c(names,paste0("w",1:(ncol(df) - length(names))))
  df <- df |> pivot_longer(-all_of(names),names_to = "Index",values_to = "Coverage") |>
    mutate(Index = Index |> str_remove("w") |> as.numeric())
  
  if(error == T & summarise_by == "mean"){
    names(df.error) <- c(names,paste0("w",1:(ncol(df.error) - length(names))))
    df.error <- df.error |> pivot_longer(-all_of(names),names_to = "Index",values_to = "Error") |>
      mutate(Index = Index |> str_remove("w") |> as.numeric())
    df <- df |> left_join(df.error,by=join_by(!!!names,Index))
  }
  
  
  p <- df |> ggplot(aes(Index,Coverage,group=!!sym(colour_by),colour=!!sym(colour_by)))
  
  if(error == T & summarise_by == "mean"){
    p <- p +geom_ribbon(aes(ymin=Coverage-Error,ymax=Coverage+Error,fill=!!sym(colour_by)),alpha=alpha_error)
  }
  
  p <- p +  #geom_vline(xintercept = vlines, colour = "darkgrey") +
    geom_line(alpha=alpha,linewidth=linewidth) +
    theme_bw() +
    labs(x = feature,y = unit, title = title) +
    scale_colour_manual(values = colmap) +
    scale_fill_manual(values = colmap) +
    scale_x_continuous(breaks=breaks,labels=labels,minor_breaks = NULL) +
    theme(axis.text.x = element_text(angle = angle, hjust=1),
          strip.background = element_rect(fill=strip_fill))
  
  if(!is.null(facet)){
    if(facet_type=="wrap"){
      p <- p + facet_wrap(facet,scales = facet_scale,nrow = facet_nrow)
    }
    else{
      p <- p + facet_grid2(facet,scales = facet_scale,independent=facet_independent)
    }
  }
  
  p
}


# Updated getFeature function

## getFeature from rootGenomics but fixed

getFeature <- function(object, start_feature = "TSS", start_flank = 0, start_exon = NULL, 
                       start_exon_boundary = NULL, start_direction = "up", end_feature = "TES", 
                       end_flank = 0, end_direction = "down", end_exon = NULL, end_exon_boundary = NULL, 
                       return = "all") 
{
  if (!inherits(object, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  allowedFeatures = c("TSS", "TES", "Exon")
  if (!start_feature %in% allowedFeatures) {
    stop("start_feature must be either 'TSS', 'TES' or 'Exon'.")
  }
  if (!end_feature %in% allowedFeatures) {
    stop("end_feature must be either 'TSS', 'TES' or 'Exon'.")
  }
  allowed_directions = c("up", "down")
  if (!start_direction %in% allowed_directions) {
    stop("start_direction should be set to either 'up' or 'down'.")
  }
  if (!end_direction %in% allowed_directions) {
    stop("stop_direction should be set to either 'up' or 'down'")
  }
  starts = list()
  object = object
  plusStrandRanges = GRanges()
  minusStrandRanges = GRanges()
  if (start_feature == "Exon") {
    if (!is.null(start_exon_boundary)) {
      exonStartNumber = start_exon
      splitIndex = object$blockCount >= exonStartNumber
      if (TRUE %in% splitIndex) {
        object = object[splitIndex]
        if (end_feature != "Exon") {
          plusStrandRanges <- subset(object, strand(object) == 
                                       "+")
          minusStrandRanges <- subset(object, strand(object) == 
                                        "-")
        }
        else {
          if (!is.null(end_exon)) {
            exonEndNumber = end_exon
            splitIndex = object$blockCount >= exonEndNumber
            if (TRUE %in% splitIndex) {
              object = object[splitIndex]
              plusStrandRanges <- subset(object, strand(object) == 
                                           "+")
              minusStrandRanges <- subset(object, strand(object) == 
                                            "-")
            }
            else {
              stop("There are no genes with that many exons. Please enter another exon end.")
            }
          }
          else {
            stop("If end_feature is set to Exon then the end_exon must not be null.")
          }
        }
      }
      else {
        stop("There are no genes with that many exons. Please enter another exon start.")
      }
    }
    else {
      stop("If the start feature is set to Exon, start_exon_boundary must be set to either '3prime' or '5prime'.")
    }
  }
  else if (start_feature == "TSS") {
    if (is.null(start_exon) & is.null(start_exon_boundary)) {
      if (end_feature != "Exon") {
        plusStrandRanges <- subset(object, strand(object) == 
                                     "+")
        minusStrandRanges <- subset(object, strand(object) == 
                                      "-")
      }
      else {
        if (!is.null(end_exon)) {
          exonEndNumber = end_exon
          splitIndex = object$blockCount >= exonEndNumber
          if (TRUE %in% splitIndex) {
            object = object[splitIndex]
            plusStrandRanges <- subset(object, strand(object) == 
                                         "+")
            minusStrandRanges <- subset(object, strand(object) == 
                                          "-")
          }
          else {
            stop("There are no genes with that many exons. Please enter another exon end.")
          }
        }
        else {
          stop("If end_feature is set to Exon then end_exon must not be null")
        }
      }
    }
    else {
      stop("If start feature is TSS then no exon start information must be entered.")
    }
  }
  else {
    if (is.null(start_exon) & is.null(start_exon_boundary)) {
      if (end_feature != "Exon") {
        plusStrandRanges <- subset(object, strand(object) == 
                                     "+")
        minusStrandRanges <- subset(object, strand(object) == 
                                      "-")
      }
      else {
        if (!is.null(end_exon)) {
          exonEndNumber = end_exon
          splitIndex = object$blockCount >= exonEndNumber
          if (TRUE %in% splitIndex) {
            object = object[splitIndex]
            plusStrandRanges <- subset(object, strand(object) == 
                                         "+")
            minusStrandRanges <- subset(object, strand(object) == 
                                          "-")
          }
          else {
            "There are no genes with that many exons. Please enter another exon end."
          }
        }
        else {
          stop("If end_feature is set to Exon then end_exon must not be null")
        }
      }
    }
    else {
      stop("If start feature is TES then no exon start information must be entered.")
    }
  }
  minusStrandTempRanges = minusStrandRanges |> as.data.frame()
  plusStrandTempRanges = plusStrandRanges |> as.data.frame()
  if (start_feature == "TSS") {
    if (start_direction == "up") {
      plusStrandTempRanges$start = start(plusStrandRanges) - 
        start_flank
      minusStrandTempRanges$end = end(minusStrandRanges) + 
        start_flank
    }
    else {
      plusStrandTempRanges$start = start(plusStrandRanges) + 
        start_flank
      
      minusStrandTempRanges$end = end(minusStrandRanges) - 
        start_flank
    }
  }
  if (start_feature == "TES") {
    if (start_direction == "up") {
      plusStrandTempRanges$start = end(plusStrandRanges) - 
        start_flank
      minusStrandTempRanges$end = start(minusStrandRanges) + 
        start_flank
    }
    else {
      plusStrandTempRanges$start = end(plusStrandRanges) + 
        start_flank
      minusStrandTempRanges$end = start(minusStrandRanges) - 
        start_flank
    }
  }
  if (start_feature == "Exon") {
    if (start_exon_boundary == "5prime") {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][start_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus, 
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, 
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus, 
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon] + 
          block_sizes_minus[[i]][start_exon]
      }
    }
    else {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes, 
                                         ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][start_exon] + 
          block_sizes_plus[[i]][start_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus, 
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, 
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus, 
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon]
      }
    }
    if (start_direction == "up") {
      plusStrandTempRanges$start = exon_starts_plus - 
        start_flank
      minusStrandTempRanges$end = exon_starts_minus + 
        start_flank
    }
    else {
      plusStrandTempRanges$start = exon_starts_plus + 
        start_flank
      minusStrandTempRanges$end = exon_starts_minus - 
        start_flank
    }
  }
  if (end_feature == "TSS") {
    if (end_direction == "up") {
      plusStrandTempRanges$end = start(plusStrandRanges) - 
        end_flank
      minusStrandTempRanges$start = end(minusStrandRanges) + 
        end_flank
    }
    else {
      plusStrandTempRanges$end = start(plusStrandRanges) + 
        end_flank
      minusStrandTempRanges$start = end(minusStrandRanges) - 
        end_flank
    }
  }
  if (end_feature == "TES") {
    if (end_direction == "up") {
      plusStrandTempRanges$end = end(plusStrandRanges) - 
        end_flank
      minusStrandTempRanges$start = start(minusStrandRanges) + 
        end_flank
    }
    else {
      plusStrandTempRanges$end = end(plusStrandRanges) + 
        end_flank
      minusStrandTempRanges$start = start(minusStrandRanges) - 
        end_flank
    }
  }
  else if (end_feature == "Exon") {
    if (end_exon_boundary == "5prime") {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][end_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus, 
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, 
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus, 
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon] + 
          block_sizes_minus[[i]][end_exon]
      }
    }
    else {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes, 
                                         ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][end_exon] + 
          block_sizes_plus[[i]][end_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts, 
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus, 
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes, 
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus, 
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon]
      }
    }
    if (end_direction == "up") {
      plusStrandTempRanges$end = exon_starts_plus - end_flank
      minusStrandTempRanges$start = exon_starts_minus + 
        end_flank
    }
    else {
      plusStrandTempRanges$end = exon_starts_plus + end_flank
      minusStrandTempRanges$start = exon_starts_minus - 
        end_flank
    }
  }
  plusStrandTempRanges <- GRanges(plusStrandTempRanges)
  minusStrandTempRanges <- GRanges(minusStrandTempRanges)
  if (return == "all") {
    allRanges = c(minusStrandTempRanges, plusStrandTempRanges)
    return(allRanges)
  }
  else if (return == "plus") {
    return(plusStrandTempRanges)
  }
  else if (return == "minus") {
    return(minusStrandTempRanges)
  }
  else {
    print("Please set a return value to 'all', 'plus' or 'minus'")
  }
}
