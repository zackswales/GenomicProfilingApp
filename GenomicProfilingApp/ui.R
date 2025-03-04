ui <- page_navbar(
  title = "Interactive Genomic Profiling",
  bg = "#E46303",
  inverse = TRUE,
  
  # Panel for home page with directions for use of application
  nav_panel(
    title = "Home",
    card(
      full_screen = FALSE,
      card_header("Introduction"),
      card_body("Interactive genomic profiling tool to probe internal genomic features. 
       Begin by navigating to the data processing page for file uploading and 
       matrix computation before moving to the visualisations page to generate outputs.")
    ),
    card(
      full_screen = FALSE,
      card_header("File Upload Requirements"),
      card_body("If intending on locating exons/introns ensure region file is either .bed12 format or .gtf/.gff format")
    ),
    card(
      full_screen = FALSE,
      card_header("Heatmap Splitting"),
      card_body("If intending to split each heatmap - ensure there are two region files (e.g. a file for RP and nRP if intending to split by RP)")
    )
  ),
  
  # Panel for data processing and matrix formation
  nav_panel(
    title = "Data Processing",
    navset_tab(
      
      # File Upload Panel
      nav_panel(
        title = "1. File Upload",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            width = 325,
            title = "File Upload",
            checkboxInput(
              inputId = "split",
              label = "Heatmap splitting",
              value = FALSE
            ),
            fileInput(
              inputId = "Region1",
              label = "Upload region files (.bed/.gtf)",
              accept = c(".bed", ".gtf"),
              multiple = TRUE
            ),
            uiOutput("Region1splitting"),
            uiOutput("conditionalRegion2"),
            uiOutput("conditionalRegion2splitting"),
            fileInput(
              inputId = "Sequence1",
              label = "Upload sequence data files (.bw)",
              accept = ".bw",
              multiple = TRUE
            ),
            radioButtons(
              inputId = "strand",
              label = "Select strandedness of data:",
              choices = list("Unstranded" = 1, "Forward" = 2, "Reverse" = 3),
              selected = 1
            ),
            helpText("If all files have the same strandedness - select unstranded")
          ),
          card(
            full_screen = FALSE,
            card_header("Uploaded Region Files"),
            card_body(textOutput("regionfile1_name"))
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
          style = "display: flex; justify-content: center; align-items: center; min-height: 50vh;", # Adjust min-height as needed
          layout_columns(
            columns = 3,
            width = "300px",
            title = "Feature Specification",
            card(
              height = "300px",
              card_header("Feature"),
              selectInput(
                inputId = "getFeature",
                label = "Specify feature of interest:",
                choices = list("Full gene" = 1, "TSS" = 2, "TES" = 3, "Exon 1 + Intron 1 + Gene body" = 4),
                selected = 1
              )
            ),
            card(
              card_header("Flank region"),
              numericInput(
                inputId = "flank",
                label = "Specify flank around feature:",
                value = 20,
                step = 10
              )
            ),
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
                card_header("Window size for features"),
                numericInput(
                  inputId = "up",
                  label = "Upstream flank:",
                  value = 20,
                  step = 10
                ),
                numericInput(
                  inputId = "exon1",
                  label = "Exon 1:",
                  value = 50,
                  step = 10
                ),
                numericInput(
                  inputId = "intron1",
                  label = "Intron 1:",
                  value = 50,
                  step = 10
                ),
                numericInput(
                  inputId = "body",
                  label = "Gene body:",
                  value = 50,
                  step = 10
                ),
                numericInput(
                  inputId = "down",
                  label = "Downstream flank",
                  value = 20,
                  step = 10
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
              label = "Smoothen data",
              value = FALSE
            ),
            actionButton(
              inputId = "matrixgeneration",
              label = "Generate Matrices"
            ),
            helpText("Warning: only click once as matrix generation takes a few seconds to complete"),
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
  ),
  
  # Panel for visualisations
  nav_panel(
    title = "Visualisations",
    navset_tab(
      # Matrix selection tab
      nav_panel(
        title = "Matrix selection",
        layout_sidebar(
          sidebar = sidebar(
            open = TRUE,
            title = "Select matrices for plotting:",
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
        title = "Heatmap",
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
)