# Source to the plottingFunctionsZack.R script
source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Source to the ui.R script
source("~/GenomicProfilingApp/GenomicProfilingApp/ui.R")

# Source to the server.R script
source("~/GenomicProfilingApp/GenomicProfilingApp/server.R")

# Run the application
shinyApp(ui = ui, server = server)
