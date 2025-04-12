ui <- page_navbar(
  title = "seqEasy",
  bg = "#BCDDDF",
  inverse = TRUE,
  useShinyjs(),
  
  tab_login$ui,
  
  # Source to the home panel module
  
  mod_home_ui("home"),
  
  # Source to the data processing panel module
  
  mod_data_processing_ui("data_processing"),
  
  # Source to the visualisations panel module
  
  mod_visualisations_ui("visualisations")
  
)