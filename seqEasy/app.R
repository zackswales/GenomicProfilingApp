## Sourcing the modules for the ui.R script
source("seqEasy/Modules/mod_home.R")
source("seqEasy/Modules/mod_data_processing_ui.R")
source("seqEasy/Modules/mod_visualisations.R")

source("seqEasy/login.R")

# Source to the ui.R script
source("seqEasy/ui.R")

# Source to the server.R script
source("seqEasy/server.R")

# Run the application
shinyApp(ui = ui, server = server)


