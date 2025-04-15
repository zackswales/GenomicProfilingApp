# Source to plotting functions script
source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")
# Source to the ggplotheatmap function
source("~/GenomicProfilingApp/R/ggplotheatmap.R")

## Sourcing the modules for the ui.R script
source("~/GenomicProfilingApp/seqEasy/Modules/mod_home.R")
source("~/GenomicProfilingApp/seqEasy/Modules/mod_data_processing_ui.R")
source("~/GenomicProfilingApp/seqEasy/Modules/mod_visualisations.R")

source("~/GenomicProfilingApp/seqEasy/login.R")

# Source to the ui.R script
source("~/GenomicProfilingApp/seqEasy/ui.R")

# Source to the server.R script
source("~/GenomicProfilingApp/seqEasy/server.R")

# Run the application
shinyApp(ui = ui, server = server)
