mod_visualisations_ui <- function(id){
  ns <- NS(id)
  nav_panel(
    title = "Visualisations",
    navset_tab(
      # Matrix selection tab
      nav_panel(
        title = "Matrix selection",
        card(
          card_header("Matrix selection"),
          card_body(
            checkboxGroupInput(
              inputId = "selectedmatrices",
              label = "Select all that apply:",
              choices = character(0)
            )
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
              label = "Upper y-axis limit for metaplot",
              value = 5000,
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
            helpText("Output plotting will take a few seconds"),
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
      nav_panel(
        title = "ggplot2 Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            title = "Heatmap Customisation",
            selectInput(
              inputId = "colorpalette",
              label = "Select colour palette",
              choices = c("viridis" = "viridis", "Red to blue" = "RdBu", "magma" = "magma", "Red to white" = "red_white"),
              selected = "viridis"
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
            )
          ),
          card(
            full_screen = TRUE,
            card_hearder = "Heatmap visualisation",
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
            conditionalPanel(
              condition = "input.split == 1",
              selectInput(
                inputId = "colourby",
                label = "Colour by:",
                choices = c("Group" = "Family", "Sample" = "Sample"),
                selected = "Family"
              ),
              selectInput(
                inputId = "facetby",
                label = "Facet by:",
                choices = c("Sample" = 1, "Group" = 2, "Sample + Group" = 3),
                selected = 2
              )
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