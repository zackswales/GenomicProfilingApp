mod_home_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    title = "Home",
    navset_tab(
      nav_panel(
        title = "Introduction",
        card(
          card_header("Welcome"),
          card_body("Welcome to seqEasy - an interactive genomic profiling app for the visualisation of a variety of types of sequencing data including ChIP-seq, RNA-seq, and DNA methylation sequencing.")
        ),
        card(
          card_body("To begin, please navigate to the step-by-step tutorial tab to learn how to utilise the app fully, and then the requirements tab for specification of the formats and file types that the app accepts")
        ),
        card(
          card_header("Example Outputs"),
          card_body(
            tags$div(
              style = "display: flex; justify-content: space-around;", # Arrange images horizontally
              tags$img(src = "enrichedHeatmap.png", height = "550px", width = "500px"), # Replace "image1.png" with your image file name
              tags$img(src = "ggplotHeatmap.png", height = "550px", width = "500px"), # Replace "image2.png" with your image file name
              tags$img(src = "averageprofile.png", height = "550px", width = "500px")  # Replace "image3.png" with your image file name
            )
          )
        )
      ),
      nav_panel(
        title = "Step-by-Step Tutorial",
        card(
          card_header("Select tutorial:"),
          card_body(
            checkboxInput(
              inputId = "dataprocessingtutorial",
              label = "Data processing",
              value = FALSE
            )
          )
        ),
        conditionalPanel(
          condition = "input.dataprocessingtutorial == 1",
          card(
            card_header("Input Selection"),
            card_body("The input selection tab is where you can upload your personal region files and fetch annotation files directly from the UCSC database. This is also where sequence data files are selected. If previously saved they can then be retrieved, as well as the ability to upload new sequence data and save them to your account. Customisation is also possible to flip strand association and assign reverse strand features to your forward strand sequence data files")
          ),
          card(
            card_header("Feature Specification"),
            card_body("This tab allows customisation of the features in your region file. There are standard options for the full gene, just around the TSS/TES, or the ability to create custom region blocks. When selecting Full Gene, there is the option to adjust the window size allowing you to control how many windows you divide the gene into.")
          ),
          card(
            card_header("Further Customisation + Matrix Generation"),
            card_body("Upon completion of customisation and upload, this tab is where you generate matrices that you will use to plot. Once matrices are generated they should be saved to your account using the Save Matrices button so that they can then be retrieved and used for plot generation. Before matrix generation there is also the option to smoothen the data.")
          ),
          card(
            card_header("Saved Matrices"),
            card_body("Once matrices have been generated and saved, this tab shows your saved matrices with the option to then clear these.")
          )
        ),
        card(
          card_body(
            checkboxInput(
              inputId = "visualisationstutorial",
              label = "Visualisations",
              value = FALSE
            )
          )
        ),
        conditionalPanel(
          condition = "input.visualisationstutorial == 1",
          card(
            card_header("Matrix selection + Splitting upload"),
            card_body("The matrix selection tab is where you select the matrices saved in your account that you would like to plot. There is also the option for uploading of a .tsv splitting file, for example a file which divides the heatmap into ribosomal and non-ribosomal proteins. Further details of the splitting file can be found on the requirements tab of the home page. ")
          ),
          card(
            card_header("enrichedHeatmap Heatmap"),
            card_body("This panel plots heatmaps using the enrichedHeatmap package. There are a variety of customisation options: selecting colour scheme for the plot, adjusting the quantiles of the data being plotted, altering the y-axis limit of the above average profile plot, and displaying row names for the genes. There are also a number of transformations which can be applied to the data including k-means clustering and a log2 transformation.")
          ),
          card(
            card_header("ggplot2 Heatmap"),
            card_body("The ggplot2 package is another package that can be used to generate heatmaps with additional customisation options. These include a broader range of colour palettes, optional average profile plots above the heatmap, automatic z-axis scaling, and the ability to add break lines to the plot.")
          ),
          card(
            card_header("Average Profile Plot"),
            card_body("This panel can be used when an isolated average profile plot is wanted rather than a heatmap. It provides customisation options for plot title/axis labels, as well as quantile control and alpha (line transparency) control.")
          ),
          card(
            card_header("Downloading Outputs"),
            card_body("For each of the three plots, in the sidebar, there are options to download the output as either a .png file or a .pdf file once the plot has been generated")
          )
        )
      ),
      nav_panel(
        title = "Requirements",
        card(
          card_header("Sequence data files"),
          card_body("When uploading stranded sequence data files - ensure the filename is in the format [name].[strand].bw - with a forward stranded file replacing [strand] with f, and a reverse strand file replacing [strand] with r. Ensure when uploading two sequence files from the same sample that the names are identical (e.g. sample1.f.bw and sample1.r.bw). If using unstranded sequence data files then ensure the filename is in the format [name].bw")
        ),
        card(
          card_header("Region files + Database fetching"),
          card_body("Multiple region files can be uploaded at the same time and used for matrix generation. However, when fetching a database annotation, they cannot be used in conjuction with uploaded region files. Ensure to only either use uploaded region files or fetched database annotations.")
        ),
        card(
          card_header("Plotting inputs"),
          card_body("When retrieving previously generated matrices, the feature specification information is not stored with the matrix. Therefore when replotting the matrices ensure to restore the Feature Specification panel inputs as to what they were when the matrices were generated. Inputs in the Input Selection tab and Further Customisation + Matrix generation tab do not need to be restored.")
        ),
        card(
          card_header("Splitting file"),
          card_body("The splitting file must be in a .tsv format with first column containing the identical gene names as the uploaded region file. The subsequent columns then contain the data which the user wants to split the heatmap by. When selecting a column to split by - ensure to not select the column which contains the gene names.")
        ),
        card(
          card_header("Example of a .tsv used for splitting"),
          card_body(
            tags$div(
              style = "display: flex; justify-content: space-around;", # Arrange images horizontally
              tags$img(src = "exampletsv.png", height = "125px", width = "125px"), # Replace "image1.png" with your image file name
            )
          )
        )
      )
    )
  )
}
