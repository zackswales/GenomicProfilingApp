library(genomation)
library(tidyverse)
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

ui <- page_navbar(
  title = "seqEasy",
  bg = "#BCDDDF",
  inverse = TRUE,
  useShinyjs(),
  navbar_options = navbar_options(cbg = "#BCDDDF", theme = "auto"),
  
  tab_login$ui,
  
  # Source to the home panel module
  
  mod_home_ui("home"),
  
  # Source to the data processing panel module
  
  mod_data_processing_ui("data_processing"),
  
  # Source to the visualisations panel module
  
  mod_visualisations_ui("visualisations")
  
)

?page_navbar
