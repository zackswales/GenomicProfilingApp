mod_home_ui <- function(id) {
  ns <- NS(id)

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
)
}