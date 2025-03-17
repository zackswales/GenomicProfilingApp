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
      card_body("Upload .tsv annotation file for heatmaps splitting containing two columns - feature names + category to split by")
    )
  ),
  
  # Panel for data processing and matrix formation
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
              inputId = "split",
              label = "Heatmap splitting",
              value = FALSE
            ),
            checkboxInput(
              inputId = "databasefetch",
              label = "Fetch database annotation",
              value = FALSE
            ),
            fileInput(
              inputId = "Region1",
              label = "Upload region files (.bed/.gtf)",
              accept = c(".bed", ".gtf"),
              multiple = TRUE
            ),
            conditionalPanel(
              condition = "input.databasefetch == 1",
              uiOutput("pickgenome"),
              uiOutput("pickgroup"),
              uiOutput("picktrack"),
              uiOutput("getannotation")
            ),
            conditionalPanel(
              condition = "input.split == 1",
              uiOutput("tsvupload")
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
              fileInput(
                inputId = "Sequence1",
                label = "Upload sequence data files (.bw)",
                accept = ".bw",
                multiple = TRUE
              ),
              textInput(
                inputId = "sequencenames",
                label = "File name for uploaded sequenece data files"
              ),
              actionButton(
                inputId = "savesequencedata",
                label = "Save sequence data files"
              )
            ),
            radioButtons(
              inputId = "strand",
              label = "Select strandedness of sequence data:",
              choices = list("Unstranded" = 1, "Forward" = 2, "Reverse" = 3),
              selected = 1
            ),
            helpText("If all files have the same strandedness - select unstranded")
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
                height = "230px",
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
                    step = 10
                  )
                )
              ),
              conditionalPanel(
                condition = "input.getFeature == 4",
                DTOutput("savedRegionsTable")
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
                  actionButton(
                    inputId = "getregion",
                    label = "Get Region"
                  ),
                  textInput(
                    inputId = "regionname",
                    label = "Region name:"
                  ),
                  selectInput(
                    inputId = "startfeature",
                    label = "Start feature:",
                    choices = c("TSS" = 1, "TES" = 2),
                    selected = 1
                  ),
                  selectInput(
                    inputId = "endfeature",
                    label = "End feature:",
                    choices = c("TSS" = 1, "TES" = 2),
                    selected = 2
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
                    inputId = "clearregions",
                    label = "Clear saved regions"
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
              label = "Smoothen data",
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
  ),
  
  # Panel for visualisations
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
)