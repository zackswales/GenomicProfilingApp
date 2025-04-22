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

mod_data_processing_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel(
    title = "Data Processing",
    navset_tab(
      
      # File Upload Panel
      nav_panel(
        title = "1. Input Selection",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 400,
            title = "File Upload",
            checkboxInput(
              inputId = "databasefetch",
              label = "Fetch database annotation",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.databasefetch == 0",
              uiOutput("Region1_ui")
            ),
            conditionalPanel(
              condition = "input.databasefetch == 1",
              uiOutput("pickgenome"),
              uiOutput("pickgroup"),
              uiOutput("picktrack"),
              uiOutput("getannotation")
            ),
            selectInput(
              inputId = "sequencedatafiles",
              label = "Select existing sequence data files:",
              choices = character(0),
              multiple = TRUE
            ),
            actionButton(
              inputId = "clearsequence",
              label = "Delete saved sequence files"
            ),
            checkboxInput(
              inputId = "newsequence",
              label = "Upload new sequence data",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.newsequence == 1",
              uiOutput("Sequence1_ui"),
              actionButton(
                inputId = "savesequencedata",
                label = "Save sequence data file"
              )
            ),
            checkboxInput(
              inputId = "flipstrand",
              label = "Flip strand association",
              value = FALSE
            )
          ),
          card(
            full_screen = FALSE,
            card_header("Region Files"),
            card_body(textOutput("regionfile1_name"))
          ),
          conditionalPanel(
            condition = "input.databasefetch == 1",
            card(
              full_screen = FALSE,
              card_header("Fetched annotation files"),
              card_body(textOutput("annotationname"))
            )
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
        title = "2. Feature Specification",
        div(
          style = "display: flex; justify-content: center; align-items: center; min-height: 50vh;",
          fluidRow(
            column(
              width = 2,
              card(
                height = "300px",
                card_header("Feature"),
                selectInput(
                  inputId = "getFeature",
                  label = "Specify feature of interest:",
                  choices = list("Full gene" = 1, "TSS" = 2, "TES" = 3, "Custom" = 4),
                  selected = 1
                )
              )
            ),
            column(
              width = 6,
              conditionalPanel(
                condition = "input.getFeature != 4",
                card(
                  card_header("Flank region"),
                  numericInput(
                    inputId = "flank",
                    label = "Specify flank around feature:",
                    value = 20,
                    step = 10,
                    min = 0
                  )
                )
              ),
              conditionalPanel(
                condition = "input.getFeature == 4",
                DTOutput("savedRegionsTable"),
                card(
                  actionButton(
                    inputId = "clearregions",
                    label = "Clear saved regions"
                  )
                )
              )
            ),
            column(
              width = 4,
              conditionalPanel(
                condition = "input.getFeature == 1",
                card(
                  card_header("Window size"),
                  sliderInput(
                    inputId = "windowsize",
                    label = "Select window size for signal aggregation:",
                    min = 1,
                    max = 20,
                    value = 1
                  )
                )
              ),
              conditionalPanel(
                condition = "input.getFeature == 4",
                card(
                  card_header("Custom specification"),
                  textInput(
                    inputId = "regionname",
                    label = "Region name:"
                  ),
                  selectInput(
                    inputId = "startfeature",
                    label = "Start feature:",
                    choices = c("TSS" = 1, "TES" = 2, "Exon" = 3),
                    selected = 1
                  ),
                  conditionalPanel(
                    condition = "input.startfeature == 3",
                    numericInput(
                      inputId = "startexon",
                      label = "Start exon:",
                      min = 1,
                      max = 100,
                      value = 1
                    ),
                    selectInput(
                      inputId = "startexonboundary",
                      label = "Start exon boundary:",
                      choices = c("5 prime" = "5prime", "3 prime" = "3prime"),
                      selected = "5prime"
                    )
                  ),
                  selectInput(
                    inputId = "endfeature",
                    label = "End feature:",
                    choices = c("TSS" = 1, "TES" = 2, "Exon" = 3),
                    selected = 2
                  ),
                  conditionalPanel(
                    condition = "input.endfeature == 3",
                    numericInput(
                      inputId = "endexon",
                      label = "End exon",
                      min = 1,
                      max = 100,
                      value = 1
                    ),
                    selectInput(
                      inputId = "endexonboundary",
                      label = "End exon boundary:",
                      choices = c("5 prime" = "5prime", "3 prime" = "3prime"),
                      selected = "3prime"
                    )
                  ),
                  numericInput(
                    inputId = "startflank",
                    label = "Start flank:",
                    min = 0,
                    max = 500,
                    value = 0,
                    step = 10
                  ),
                  numericInput(
                    inputId = "endflank",
                    label = "End flank:",
                    min = 0,
                    max = 500,
                    value = 0,
                    step = 10
                  ),
                  selectInput(
                    inputId = "startdirection",
                    label = "Start direction:",
                    choices = c("Up" = 1, "Down" =2),
                    selected = 1
                  ),
                  selectInput(
                    inputId = "enddirection",
                    label = "End direction:",
                    choices = c("Up" = 1, "Down" = 2),
                    selected = 2
                  ),
                  numericInput(
                    inputId = "winsize",
                    label = "Window size for region:",
                    min = 0,
                    value = 20,
                    step = 10
                  ),
                  actionButton(
                    inputId = "getregion",
                    label = "Get Region"
                  )
                )
              )
            )
          )
        )
      ),
      
      # Further Customisation Panel
      nav_panel(
        title = "3. Further Customisation + Matrix Generation",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "Further Customisation",
            checkboxInput(
              inputId = "smooth",
              label = "Smooth data",
              value = FALSE
            ),
            actionButton(
              inputId = "matrixgeneration",
              label = "Generate Matrices"
            ),
            verbatimTextOutput("matrixResult"),
            actionButton(
              inputId = "savematrices",
              label = "Save Matrices"
            )
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
            title = "Saved Matrices",
            actionButton(
              inputId = "clearmatrices",
              label = "Clear Saved Matrices"
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Saved Matrices"),
            card_body(textOutput("savedmatrices"))
          )
        )
      )
    )
  )
}