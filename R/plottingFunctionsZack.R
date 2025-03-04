library(genomation)
library(tidyverse)
library(EnrichedHeatmap)
library(rootGenomics)
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
library(ggplot2)
library(patchwork)
library(scales)

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
matList <- function(bwf,bwr,names,grl,wins=list("Gene"=10),mode="coverage",output="norm.matrix",strand="rev",smooth=FALSE,
                    extend=0,w,include_target=T,target_ratio=0.5,k=10,keep,log=F){
  
  fbw <- bwf
  #fbw <- bwf |> future_map(~import(.x, format = "BigWig"))
  #names(fbw) <- names
  if(log==T){
    fbw <- fbw |> map(function(x){
      x$score <- log(x$score)
      x
    })
    names(fbw) <- names(bwf)
  }
  
  if(strand %in% c("for","rev")){
    rbw <- bwr
    #rbw <- bwr |> future_map(~import(.x, format = "BigWig"))
    #names(rbw) <- names
    if(log==T){
      rbw <- rbw |> map(function(x){
        x$score <- log(x$score)
        x
      })
      names(rbw) <- names(bwr)
    }
    
    fgrl = grl |> map(~subset(.x,strand(.x) == "+"))
    rgrl = grl |> map(~subset(.x,strand(.x) == "-"))
  }
  
  if(length(wins) > 1){ ## Length of features to combine
    matl <- names |> future_map(function(s){
      
      if(strand == "rev"){
        fmat <- fgrl |> imap(~normalizeToMatrix(rbw[[s]],.x,extend=0,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth)) |>
          bind_cols() |> as.data.frame() |> as.matrix() |> as.normalizedMatrix(k_target = sum(wins),extend=0)
        rmat<- rgrl |> imap(~normalizeToMatrix(fbw[[s]],.x,extend=0,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth)) |>
          bind_cols() |> as.data.frame() |> as.matrix() |> as.normalizedMatrix(k_target = sum(wins),extend=0)
        
        mat<-rbind(fmat,rmat)
        rownames(mat) <- c(fgrl[[1]]$name,rgrl[[1]]$name)
      }
      else if(strand == "for"){
        fmat <- fgrl |> imap(~normalizeToMatrix(fbw[[s]],.x,extend=0,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth)) |>
          bind_cols() |> as.data.frame() |> as.matrix() |> as.normalizedMatrix(k_target = sum(wins),extend=0)
        rmat<- rgrl |> imap(~normalizeToMatrix(rbw[[s]],.x,extend=0,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth)) |>
          bind_cols() |> as.data.frame() |> as.matrix() |> as.normalizedMatrix(k_target = sum(wins),extend=0)
        
        mat<-rbind(fmat,rmat)
        rownames(mat) <- c(fgrl[[1]]$name,rgrl[[1]]$name)
      }
      else if(strand == "no"){
        mat <- grl |> imap(~normalizeToMatrix(fbw[[s]],.x,extend=0,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth)) |>
          bind_cols() |> as.data.frame() |> as.matrix() |> as.normalizedMatrix(k_target = sum(wins),extend=0)   
        rownames(mat) <- grl[[1]]$name
      }
      
      mat
      
    },.progress=T)
  }
  else{  ## Single feature - regular normalizeToMatrix
    matl <- names |> future_map(function(s){
      
      if(strand == "rev"){
        fmat <- normalizeToMatrix(rbw[[s]],fgrl[[1]],extend=extend,value_column = "score",k = wins[1],mean_mode = mode,smooth = smooth,
                                  w = w,include_target = include_target,target_ratio = target_ratio)
        rmat<- normalizeToMatrix(fbw[[s]],rgrl[[1]],extend=extend,value_column = "score",k = wins[1],mean_mode = mode,smooth = smooth,
                                 w = w,include_target = include_target,target_ratio = target_ratio)
        
        mat<-rbind(fmat,rmat)
        rownames(mat) <- c(fgrl[[1]]$name,rgrl[[1]]$name)
      }
      else if(strand == "for"){
        fmat <- normalizeToMatrix(fbw[[s]],fgrl[[1]],extend=extend,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth,
                                  w = w,include_target = include_target,target_ratio = target_ratio)
        rmat <- normalizeToMatrix(rbw[[s]],rgrl[[1]],extend=extend,value_column = "score",k = wins[1],mean_mode = mode,smooth = smooth,
                                  w = w,include_target = include_target,target_ratio = target_ratio)
        
        mat<-rbind(fmat,rmat)
        rownames(mat) <- c(fgrl[[1]]$name,rgrl[[1]]$name)
      }
      else if(strand == "no"){
        mat <- normalizeToMatrix(fbw[[s]],grl[[1]],extend=extend,value_column = "score",k = wins[.y],mean_mode = mode,smooth = smooth,
                                 w = w,include_target = include_target,target_ratio = target_ratio)
        # Issue arising in k = wins[.y] changes the for and no strand from .y to 1
        rownames(mat) <- grl[[1]]$name
      }
      
      mat
      
    },.progress=T)
  }
  
  names(matl) <- names
  
  if(output == "norm.matrix"){
    matl
  }
  else if(output == "matrix"){
    matl |> map(~as.data.frame(.x) |> as.matrix())
  }
}    

# encrichedHeatmap list function

hmList <- function(matl,wins,split=NULL,split_cols,max_quantile=0.99,min_quantile=0,col_fun="red",show_row_names=T,win_labels=NULL,ylim=NULL,summarise_by="mean",axis_labels=""){
  
  ## ColourMap
  reds = RColorBrewer::brewer.pal(n=9,name="Reds")
  
  common_min = quantile(unlist(matl),min_quantile)
  common_max = quantile(unlist(matl),max_quantile)
  if(col_fun == "red"){
    col_fun = circlize::colorRamp2(c(common_min, common_max), c("white", "red"))
  }
  else if(col_fun == "bl2rd"){
    col_fun = circlize::colorRamp2(c(common_min,0,common_max), c("blue","white", "red"))
  }
  else if(col_fun == "red0"){
    col_fun = circlize::colorRamp2(c(0,1,common_max), c("white",reds[3],reds[7]))
  }
  
  if(is.null(ylim)){
    ymin = min(0,common_min)
    #ymin = matl |> map(~colMeans(.x) |> min()) |> unlist() |> min()
    ymax = common_max
    offset = 1:(ymax-ymin) |> quantile(0.01)
    ylim = c(ymin-offset,ymax+offset)
  }
  
  features = vector()
  fcols = brewer.pal(9,"Set1")[1:length(wins)]
  names(fcols) = names(wins)
  for(i in 1:length(wins)){
    features <- append(features,rep(names(wins)[i],wins[i]))
  }
  
  if(is.null(win_labels)){
    win_labels = names(wins)
  }
  
  if(!is.null(split)){
    
    rowAnno = rowAnnotation(df=split,col=split_cols,show_legend = F)
    
    hml <- matl |> imap(~EnrichedHeatmap(.x, name = .y, column_title = .y , col = col_fun,
                                         show_row_names = show_row_names,
                                         row_names_gp=gpar(fontsize=5),
                                         axis_name = axis_labels,
                                         row_split = split,
                                         column_title_gp = gpar(fontsize = 10),
                                         #show_heatmap_legend = FALSE,
                                         cluster_rows = F,
                                         top_annotation = c(HeatmapAnnotation(Features = features,show_legend = T,border = T,col = list(Features=fcols),show_annotation_name = F,
                                                                              annotation_legend_param = list(at = names(wins),labels=win_labels)),
                                                            HeatmapAnnotation(enriched = anno_enriched(value = summarise_by,gp = gpar(fontsize=4,lwd=2,col=split_cols[[1]]),
                                                                                                       axis_param = list(side = "left",facing = "inside"),
                                                                                                       ylim = ylim)), 
                                                            gap = unit(2,"mm"))
    ))
    hml$rowAnno = rowAnno
    hml
  }
  else{
    hml <- matl |> imap(~EnrichedHeatmap(.x, name = .y, column_title = .y , col = col_fun,
                                         show_row_names = show_row_names,
                                         row_names_gp=gpar(fontsize=5),
                                         axis_name = axis_labels,
                                         column_title_gp = gpar(fontsize = 10),
                                         #show_heatmap_legend = FALSE,
                                         cluster_rows = F,
                                         top_annotation = c(HeatmapAnnotation(Features = features,show_legend = T,border = T,col = list(Features=fcols),show_annotation_name = F,
                                                                              annotation_legend_param = list(at = names(wins),labels=win_labels)),
                                                            HeatmapAnnotation(enriched = anno_enriched(gp = gpar(fontsize=4,lwd=2),
                                                                                                       axis_param = list(side = "left",facing = "inside"),
                                                                                                       ylim = ylim)), 
                                                            gap = unit(2,"mm"))
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
  
  
  p <- df |> ggplot(aes(Index,Coverage,group=Sample,colour=!!sym(colour_by)))
  
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


# Creating a function for converting heatmap to a ggplot object

gghml <- function(matl, wins, split = NULL, split_cols = NULL, max_quantile = 0.99, min_quantile = 0, col_fun = "red", show_row_names = TRUE, win_labels = NULL, ylim = NULL, summarise_by = "mean", axis_labels = "") {
  
  reds <- brewer.pal(n = 9, name = "Reds")
  
  common_min <- quantile(unlist(matl), min_quantile)
  common_max <- quantile(unlist(matl), max_quantile)
  
  if (col_fun == "red") {
    color_scale <- colorRampPalette(c("white", "red"))(100)
  } else if (col_fun == "bl2rd") {
    color_scale <- colorRampPalette(c("blue", "white", "red"))(100)
  } else if (col_fun == "red0") {
    color_scale <- colorRampPalette(c("white", reds[3], reds[7]))(100)
  }
  
  if (is.null(ylim)) {
    ymin <- min(0, common_min)
    ymax <- common_max
    offset <- quantile(1:(ymax - ymin), 0.01)
    ylim <- c(ymin - offset, ymax + offset)
  }
  
  features <- unlist(lapply(names(wins), function(name) rep(name, wins[name])))
  features <- rev(features)
  
  fcols <- brewer.pal(9, "Set1")[1:length(wins)]
  fcols <- rev(fcols)
  names(fcols) <- names(wins)
  
  if (is.null(win_labels)) {
    win_labels <- names(wins)
    win_labels <- rev(win_labels)
  }
  if(!is.null(win_labels)) {
    win_labels <- rev(win_labels)
  }
  
  plot_list <- lapply(names(matl), function(name) {
    mat <- matl[[name]]
    
    mat <- mat[, rev(seq_len(ncol(mat)))]
    
    df <- as.data.frame(mat)
    df$row_name <- rownames(mat)
    df_long <- df %>%
      pivot_longer(cols = -row_name, names_to = "column", values_to = "value")
    
    df_long$column <- factor(df_long$column, levels = rev(unique(df_long$column)))
    
    if (!is.null(split)) {
      df_long <- df_long %>%
        left_join(data.frame(row_name = rownames(split), split = split[, 1]), by = "row_name")
    }
    
    summary_values <- rowMeans(mat)
    
    summary_df <- data.frame(row_name = rownames(mat), summary_value = summary_values)
    
    df_long <- df_long %>%
      left_join(summary_df, by = "row_name")
    
    df_long$scaled_value <- pmin(pmax(df_long$value, common_min), common_max)
    df_long$color_index <- round(rescale(df_long$scaled_value, to = c(1, 100)))
    df_long$fill_color <- color_scale[df_long$color_index]
    
    p <- ggplot(df_long, aes(x = column, y = row_name, fill = fill_color)) +
      geom_tile() +
      scale_fill_identity() +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = if (show_row_names) element_text(size = 5) else element_blank(),
        axis.ticks.y = if (show_row_names) element_line() else element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)
      ) +
      ggtitle(name)
    
    feature_df <- data.frame(column = colnames(mat), feature = features)
    feature_df$column <- factor(feature_df$column, levels = rev(unique(feature_df$column)))
    
    feature_annotation <- ggplot(feature_df, aes(x = column, y = 1, fill = feature)) +
      geom_tile() +
      scale_fill_manual(values = fcols, breaks = names(fcols), labels = win_labels) +
      theme_void() +
      theme(legend.position = "right", legend.justification = "top", legend.title = element_blank())
    
    x_axis_labels <- c("-50b", "TSS", "TES", "+50b")
    x_axis_df <- data.frame(column = colnames(mat), x_label = rep("", ncol(mat)))
    
    tss_pos <- ncol(mat) - cumsum(wins)[1] + 1
    tes_pos <- ncol(mat) - (cumsum(wins)[1] + wins[2]) + 1
    
    x_axis_df$x_label[ncol(mat)] <- x_axis_labels[1]
    x_axis_df$x_label[tss_pos] <- x_axis_labels[2]
    x_axis_df$x_label[tes_pos] <- x_axis_labels[3]
    x_axis_df$x_label[1] <- x_axis_labels[4]
    
    x_axis_df$column <- factor(x_axis_df$column, levels = rev(unique(x_axis_df$column)))
    
    x_axis_annotation <- ggplot(x_axis_df, aes(x = column, y = 1, label = x_label)) +
      geom_text(size = 3) +
      theme_void()
    
    combined_plot <- feature_annotation + p + x_axis_annotation + plot_layout(ncol = 1, heights = c(0.1, 1, 0.1))
    
    if (!is.null(split)) {
      split_annotation <- ggplot(df_long %>% distinct(row_name, split), aes(x = 1, y = row_name, fill = split)) +
        geom_tile() +
        scale_fill_manual(values = split_cols[[1]]) +
        theme_void() +
        theme(legend.position = "none")
      combined_plot <- combined_plot + split_annotation + plot_layout(widths = c(1, 0.05))
    }
    
    combined_plot
  })
  
  return(plot_list)
}