library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(RColorBrewer)

plotggplotHeatmap <- function(matl, color_palette = "red_white", zMin = NULL, zMax = NULL,
                              xlab = NULL, ylab = NULL, fill_label = "Signal", title = NULL,
                              wins, break_labels, average_profile = FALSE, k_clusters = NULL, log2 = FALSE, dottedlines = FALSE,
                              split = NULL) {
  if (!is.list(matl) || !all(sapply(matl, is.matrix))) {
    stop("Input 'matl' must be a list of matrices.")
  }
  if (!is.numeric(wins) || is.null(names(wins))) {
    stop("Argument 'wins' must be a named numeric vector.")
  }
  if (!is.character(break_labels) || !is.vector(break_labels)) {
    stop("Argument 'break_labels' must be a character vector.")
  }
  
  # Apply log2 transform if log2 is TRUE
  if (log2) {
    matl <- lapply(matl, function(mat) {
      mat <- mat + 1
      mat <- log2(mat)
      return(mat)
    })
  }
  
  # Compute global zMin and zMax if set to "auto"
  all_values <- unlist(matl, use.names = FALSE)
  if (is.character(zMin) && zMin == "auto") {
    zMin <- quantile(all_values, 0.01, na.rm = TRUE)
  }
  if (is.character(zMax) && zMax == "auto") {
    zMax <- quantile(all_values, 0.99, na.rm = TRUE)
  }
  
  heatmap_list <- lapply(seq_along(matl), function(i) {
    temp_df <- as.data.frame(matl[[i]])
    temp_df$gene_id <- rownames(matl[[i]]) # Assuming rownames are gene IDs
    if (is.null(temp_df$gene_id)) {
      temp_df$gene_id <- paste0("Gene", seq_len(nrow(temp_df)))
    }
    temp_df_long <- temp_df %>%
      pivot_longer(cols = -gene_id, names_to = "Bin", values_to = "Value")
    temp_df_long$Bin <- factor(temp_df_long$Bin, levels = unique(temp_df_long$Bin))
    
    matrix_name <- names(matl)[i]
    if (is.null(matrix_name)) {
      matrix_name <- paste0("Matrix", i)
    }
    
    bin_names <- levels(temp_df_long$Bin)
    breaks <- cumsum(wins)
    x_breaks <- c(1, breaks)
    x_labels <- break_labels
    
    if (length(x_labels) != length(x_breaks)) {
      stop("The length of 'break_labels' must be equal to the number of breaks defined by 'wins' plus one (for the start).")
    }
    
    base_heatmap <- ggplot(temp_df_long, aes(x = Bin, y = gene_id, fill = Value)) +
      geom_raster() +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.ticks.x = element_line(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.spacing = unit(0, "pt"),
        plot.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5) # Add border
      ) +
      scale_x_discrete(breaks = bin_names[x_breaks], labels = x_labels) +
      labs(x = xlab, y = ylab, fill = matrix_name) #Matrix name as fill label.
    
    if(dottedlines) {
      base_heatmap <- base_heatmap + geom_vline(xintercept = x_breaks, linetype = "dotted") #Add vertical lines if dottedlines is TRUE
    }
    
    # Apply color palette
    if (color_palette == "red_white") {
      base_heatmap <- base_heatmap + scale_fill_gradient(low = "white", high = "red", na.value = "white", limits = c(zMin, zMax))
    } else if (is.character(color_palette)) {
      if (color_palette %in% c("viridis", "magma", "inferno", "cividis", "Blues", "Greens", "Reds")) {
        base_heatmap <- base_heatmap + scale_fill_viridis_c(option = color_palette, na.value = "white", limits = c(zMin, zMax))
      } else if (color_palette %in% c("RdBu", "RdYlBu", "PuOr", "BrBG")) {
        base_heatmap <- base_heatmap + scale_fill_distiller(palette = color_palette, direction = 1, na.value = "white", limits = c(zMin, zMax))
      } else {
        warning("Invalid named color palette provided. Using 'viridis' as default.")
        base_heatmap <- base_heatmap + scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
      }
    } else if (is.vector(color_palette) && all(sapply(color_palette, is.character))) {
      tryCatch({
        base_heatmap <- base_heatmap + scale_fill_gradientn(colors = color_palette, na.value = "white", limits = c(zMin, zMax))
      }, error = function(e) {
        warning("Invalid custom color palette provided. Using 'viridis' as default.")
        base_heatmap <- base_heatmap + scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
      })
    } else {
      warning("Invalid color palette input. Using 'viridis' as default.")
      base_heatmap <- base_heatmap + scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
    }
    
    final_heatmap <- base_heatmap
    
    # K-means clustering
    if (!is.null(k_clusters)) {
      gene_profiles <- as.matrix(temp_df[, !(names(temp_df) %in% c("gene_id"))])
      rownames(gene_profiles) <- temp_df$gene_id
      kmeans_result <- kmeans(gene_profiles, centers = k_clusters)
      cluster_df <- data.frame(gene_id = names(kmeans_result$cluster), Cluster = factor(kmeans_result$cluster))
      
      merged_df_long <- temp_df_long %>%
        left_join(cluster_df, by = "gene_id")
      
      final_heatmap <- ggplot(merged_df_long, aes(x = Bin, y = gene_id, fill = Value)) +
        geom_raster() +
        facet_grid(Cluster ~ ., scales = "free_y", space = "free_y", switch = "y") + # Add facet_grid
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
          axis.ticks.x = element_line(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          panel.spacing = unit(0, "pt"),
          plot.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add border
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12) # Move facet labels to the left
        ) +
        scale_x_discrete(breaks = bin_names[x_breaks], labels = x_labels) +
        labs(x = xlab, y = ylab, fill = matrix_name) + #Matrix name as fill label.
        # Apply color palette
        if (color_palette == "red_white") {
          scale_fill_gradient(low = "white", high = "red", na.value = "white", limits = c(zMin, zMax))
        } else if (is.character(color_palette)) {
          if (color_palette %in% c("viridis", "magma", "inferno", "cividis", "Blues", "Greens", "Reds")) {
            scale_fill_viridis_c(option = color_palette, na.value = "white", limits = c(zMin, zMax))
          } else if (color_palette %in% c("RdBu", "RdYlBu", "PuOr", "BrBG")) {
            scale_fill_distiller(palette = color_palette, direction = 1, na.value = "white", limits = c(zMin, zMax))
          } else {
            warning("Invalid named color palette provided. Using 'viridis' as default.")
            scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
          }
        } else if (is.vector(color_palette) && all(sapply(color_palette, is.character))) {
          tryCatch({
            scale_fill_gradientn(colors = color_palette, na.value = "white", limits = c(zMin, zMax))
          }, error = function(e) {
            warning("Invalid custom color palette provided. Using 'viridis' as default.")
            scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
          })
        } else {
          warning("Invalid color palette input. Using 'viridis' as default.")
          scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
        }
      
      if(dottedlines) {
        final_heatmap <- final_heatmap + geom_vline(xintercept = x_breaks, linetype = "dotted")
      }
    }
    
    # Split functionality
    if (!is.null(split)) {
      if (!is.data.frame(split) || ncol(split) != 1) {
        stop("The 'split' argument must be a data frame with one column and rownames as gene_ids.")
      }
      if (is.null(rownames(split))) {
        stop("The 'split' data frame must have rownames representing gene_ids.")
      }
      colnames(split) <- "split_factor"
      split$gene_id <- rownames(split)
      
      merged_df_long <- temp_df_long %>%
        left_join(split, by = "gene_id")
      
      if (!is.factor(merged_df_long$split_factor)) {
        merged_df_long$split_factor <- factor(merged_df_long$split_factor)
      }
      
      final_heatmap <- ggplot(merged_df_long %>% arrange(split_factor), aes(x = Bin, y = gene_id, fill = Value)) +
        geom_raster() +
        facet_grid(split_factor ~ ., scales = "free_y", space = "free_y", switch = "y") +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
          axis.ticks.x = element_line(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          panel.spacing = unit(0, "pt"),
          plot.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          strip.text.y.left = element_text(angle = 0, hjust = 1)
        ) +
        scale_x_discrete(breaks = bin_names[x_breaks], labels = x_labels) +
        labs(x = xlab, y = ylab, fill = matrix_name) +
        # Apply color palette
        if (color_palette == "red_white") {
          scale_fill_gradient(low = "white", high = "red", na.value = "white", limits = c(zMin, zMax))
        } else if (is.character(color_palette)) {
          if (color_palette %in% c("viridis", "magma", "inferno", "cividis", "Blues", "Greens", "Reds")) {
            scale_fill_viridis_c(option = color_palette, na.value = "white", limits = c(zMin, zMax))
          } else if (color_palette %in% c("RdBu", "RdYlBu", "PuOr", "BrBG")) {
            scale_fill_distiller(palette = color_palette, direction = 1, na.value = "white", limits = c(zMin, zMax))
          } else {
            warning("Invalid named color palette provided. Using 'viridis' as default.")
            scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
          }
        } else if (is.vector(color_palette) && all(sapply(color_palette, is.character))) {
          tryCatch({
            scale_fill_gradientn(colors = color_palette, na.value = "white", limits = c(zMin, zMax))
          }, error = function(e) {
            warning("Invalid custom color palette provided. Using 'viridis' as default.")
            scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
          })
        } else {
          warning("Invalid color palette input. Using 'viridis' as default.")
          scale_fill_viridis_c(na.value = "white", limits = c(zMin, zMax))
        }
      
      if(dottedlines) {
        final_heatmap <- final_heatmap + geom_vline(xintercept = x_breaks, linetype = "dotted")
      }
    }
    
    # If average_profile is TRUE, stack profile plot above heatmap
    if (average_profile) {
      p_average <- plotBasicAverageProfile(matl[[i]], title = matrix_name, bin_names = bin_names, x_breaks = x_breaks, x_labels = x_labels) +
        theme(plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.ticks.x = element_line(),
              plot.margin = margin(0, 0, 0, 0, "pt"),
              panel.spacing = unit(0, "pt"),
              plot.background = element_blank(),
              axis.text.y = element_text(size = 10))
      
      combined_plot <- p_average / final_heatmap +
        plot_layout(heights = c(0.5, 5)) &
        theme(plot.margin = margin(0, 0, 0, 0, "pt"),
              panel.spacing = unit(0, "pt"))
      return(combined_plot)
    } else {
      final_heatmap <- final_heatmap + ggtitle(matrix_name) +
        theme(plot.title = element_text(size = 20))
      return(final_heatmap)
    }
  })
  
  # Ensure only one colorbar by extracting it from the first plot and removing from others
  if (length(heatmap_list) > 1) {
    
    combined_plot <- wrap_plots(plotlist = heatmap_list, ncol = length(heatmap_list)) + plot_layout(guides = "collect")
    return(combined_plot)
    
  } else {
    return(heatmap_list[[1]])
  }
}

plotBasicAverageProfile <- function(matrix, title = "Average Profile", bin_names = NULL, x_breaks = NULL, x_labels = NULL) {
  df <- as.data.frame(matrix)
  df$Gene <- seq_len(nrow(df))
  df_long <- df %>%
    pivot_longer(-Gene, names_to = "Bin", values_to = "Value")
  df_long$Bin <- factor(df_long$Bin, levels = unique(df_long$Bin))
  
  average_profile <- df_long %>%
    group_by(Bin) %>%
    summarise(Average_Value = mean(Value, na.rm = TRUE))
  
  p <- ggplot(average_profile, aes(x = Bin, y = Average_Value, group = 1)) +
    geom_line() +
    ggtitle(title) +
    labs(x = "", y = "") +
    theme_classic()
  
  if (!is.null(bin_names) && !is.null(x_breaks) && !is.null(x_labels)) {
    p <- p + scale_x_discrete(breaks = bin_names[x_breaks], labels = x_labels) +
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.ticks.x = element_line())
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  return(p)
}