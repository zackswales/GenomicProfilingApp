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

mod_visualisations_ui <- function(id){
  ns <- NS(id)
  nav_panel(
    title = "Visualisations",
    navset_tab(
      # Matrix selection tab
      nav_panel(
        title = "Matrix selection + Splitting Upload",
        card(
          card_header("Matrix selection"),
          card_body(
            checkboxGroupInput(
              inputId = "selectedmatrices",
              label = "Select all that apply:",
              choices = character(0)
            )
          )
        ),
        card(
          card_header("Plot splitting"),
          card_body(
            checkboxInput(
              inputId = "split",
              label = "Heatmap splitting",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("tsvupload")
            ),
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("splitselect")
            ),
          )
        )
      ),
      
      # Heatmap Panel
      nav_panel(
        title = "EnrichedHeatmap Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE, 
            title = "Heatmap Customisation",
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("filterselect")
            ),
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
              label = "Upper y-axis limit for average profile plot",
              value = 60,
              step = 50
            ),
            checkboxInput(
              inputId = "showrownames",
              label = "Show row names",
              value = FALSE
            ),
            checkboxInput(
              inputId = "kmeansclustering",
              label = "k-Means clustering",
              value = FALSE
            ),
            checkboxInput(
              inputId = "logenriched",
              label = "Apply log2 transformation",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.kmeansclustering == 1",
              numericInput(
                inputId = "row_km",
                label = "Select number of clusters:",
                min = 1,
                value = 1
              )
            ),
            actionButton(
              inputId = "heatmapplotbutton",
              label = "Plot Output"
            ),
            helpText("Output plotting/downloading will take a few seconds"),
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
      
      # ggplot heatmap
      nav_panel(
        title = "ggplot2 Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            title = "Heatmap Customisation",
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("ggfilterselect")
            ),
            selectInput(
              inputId = "colorpalette",
              label = "Select colour palette",
              choices = c("Viridis" = "viridis", 
                          "Red-Blue" = "RdBu", 
                          "Magma" = "magma", 
                          "Red-White" = "red_white",
                          "Inferno" = "inferno",
                          "Cividis" = "cividis",
                          "Red-Yellow-Blue" = "RdYlBu",
                          "Purple-Orange" = "PuOr",
                          "Brown-Blue/Green" = "BrBG"),
              selected = "red_white"
            ),
            checkboxInput(
              inputId = "averageprofile",
              label = "Add average profile plot",
              value = FALSE
            ),
            checkboxInput(
              inputId = "autoz",
              label = "Auto z-axis scaling",
              value = FALSE
            ),
            checkboxInput(
              inputId = "log2",
              label = "Apply log2 to data",
              value = FALSE
            ),
            checkboxInput(
              inputId = "dottedlines",
              label = "Add break lines to heatmap",
              value = FALSE
            ),
            actionButton(
              inputId = "plotggheatmap",
              label = "Plot Output"
            ),
            helpText("Output plotting will take a few seconds"),
            helpText("Use download buttons after clicking Plot Output"),
            downloadButton(
              outputId = "ggheatmapdownloadpng",
              label = "Download as .png"
            ),
            downloadButton(
              outputId = "ggheatmapdownloadpdf",
              label = "Download as .pdf"
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Heatmap Visualisation"),
            card_body(plotOutput("ggplotheatmap"),
                      height = 1000)
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
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("columnstoinclude"),
              uiOutput("colouringby"),
              checkboxInput(
                inputId = "facetchoice",
                label = "Facet plot",
                value = FALSE
              ),
              conditionalPanel(
                condition = "input.facetchoice == 1",
                uiOutput("facetingby")
              )
            ),
            textInput(
              inputId = "plottitle",
              label = "Enter plot title:",
              value = "Average profile plot"
            ),
            textInput(
              inputId = "unit",
              label = "Y-axis label:",
              value = "Coverage (BPM)"
            ),
            textInput(
              inputId = "feature",
              label = "X-axis label:",
              value = "Gene"
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
}