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
      
      # ggplot Heatmap Panel
      nav_panel(
        title = "ggplot Heatmap",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            title = "Heatmap Customisation",
            radioButtons(
              inputId = "ggheatmap_col_fun",
              label = "Select colour scheme:",
              choices = list("White to red" = 1, "Blue to red" = 2, "Red scale" = 3),
              selected = 1
            ),
            sliderInput(
              inputId = "ggheatmapquantiles",
              label = "Set quantile range:",
              min = 0,
              max = 0.99,
              value = c(0,0.99)
            ),
            checkboxInput(
              inputId = "ggshowrownames",
              label = "Show row names",
              value = FALSE
            ),
            actionButton(
              inputId = "ggplotheatmapplotbutton",
              label = "Plot Output"
            ),
            helpText("Output plotting will take a few seconds"),
            helpText("Use download buttons after clicking Plot Output"),
            downloadButton(
              outputId = "ggplotheatmapdownloadpng",
              label = "Download as .png"
            ),
            downloadButton(
              outputId = "ggplotheatmapdownloadpdf",
              label = "Download as .pdf"
            )
          ),
          card(
            full_screen = TRUE,
            card_header("Heatmap Visualisation"),
            card_body(plotOutput("ggplotHeatmapPlot"))
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