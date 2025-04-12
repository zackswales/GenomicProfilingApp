# Source to plotting functions script
source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")
# Source to the ggplotheatmap function
source("~/GenomicProfilingApp/R/ggplotheatmap.R")

## Sourcing the modules for the ui.R script
source("~/GenomicProfilingApp/seqEasy/Modules/mod_home.R")
source("~/GenomicProfilingApp/seqEasy/Modules/mod_data_processing_ui.R")
source("~/GenomicProfilingApp/seqEasy/Modules/mod_visualisations.R")

### Login code

user.access <- read.csv('~/Research Project/register.csv')

# Define the base directory where user folders will be stored
base_dir <- '~/GenomicProfilingApp/seqEasy/users'
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# Initialize login system
tab_login <- list()

# UI for the login page
tab_login$ui <- tabPanel(
  title = 'Login',
  # Add logo to the login page
  tags$div(style = "text-align: center;", 
           tags$img(src = 'logo.png', height = "150px", alt = "App Logo"),
           tags$h3("Welcome to seqEasy")  # Changed from "Login to continue" to "Welcome to seqEasy"
  )
)

login_dialog <- modalDialog(
  title = 'Login to continue',
  footer = actionButton('tab_login.login','Login'),
  textInput('tab_login.username','Username'),
  passwordInput('tab_login.password','Password'),
  tags$div(class = 'warn-text', textOutput('tab_login.login_msg'))
)

tab_login$server <- function(input, output, session) {
  
  # Show login dialog box when initiated
  showModal(login_dialog)
  
  observeEvent(input$tab_login.login, {
    username <- input$tab_login.username
    password <- input$tab_login.password
    row.names(user.access) <- user.access$Name
    
    # Validate login credentials
    if (username %in% user.access[, 'Name']) {
      if (password == user.access[username, 'Password']) {
        # Successfully log in
        removeModal() # Remove login dialog
        
        # Define user-specific directory
        user_dir <- file.path(base_dir, username)
        if (!dir.exists(user_dir)) {
          dir.create(user_dir, recursive = TRUE)
        }
        
        # Store user directory in session for later access
        session$userData$user_dir <- user_dir
        
      } else {
        # Password incorrect, show the warning message
        output$tab_login.login_msg <- renderText('Incorrect Password')
        shinyjs::show('tab_login.login_msg')
        shinyjs::delay(1000, shinyjs::hide('tab_login.login_msg'))
      }
    } else {
      # Username not found, show the warning message
      output$tab_login.login_msg <- renderText('Username Not Found')
      shinyjs::show('tab_login.login_msg')
      shinyjs::delay(1000, shinyjs::hide('tab_login.login_msg'))
    }
  })
}

# Source to the ui.R script
source("~/GenomicProfilingApp/seqEasy/ui.R")

# Source to the server.R script
source("~/GenomicProfilingApp/seqEasy/server.R")

# Run the application
shinyApp(ui = ui, server = server)
