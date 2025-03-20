# Source to plotting functions script
source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

## Sourcing the modules for the ui.R script
source("~/GenomicProfilingApp/GenomicProfilingApp/Modules/mod_home.R")
source("~/GenomicProfilingApp/GenomicProfilingApp/Modules/mod_data_processing_ui.R")
source("~/GenomicProfilingApp/GenomicProfilingApp/Modules/mod_visualisations.R")

### Login code

user.access <- read.csv('~/Research Project/register.csv')

# Define the base directory where user folders will be stored
base_dir <- '~/GenomicProfilingApp/GenomicProfilingApp/users'
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}

# Initialize login system
tab_login <- list()

tab_login$ui <- tabPanel(
  title = 'Login',
  shinyjs::hidden(tags$div(
    id = 'tab_login.welcome_div',
    class = 'login-text', 
    textOutput('tab_login.welcome_text', container = tags$h2))
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
        output$tab_login.welcome_text <- renderText(glue::glue('Welcome, {username}'))
        shinyjs::show('tab_login.welcome_div') # Show welcome message
        
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
source("~/GenomicProfilingApp/GenomicProfilingApp/ui.R")

# Source to the server.R script
source("~/GenomicProfilingApp/GenomicProfilingApp/server.R")

# Run the application
shinyApp(ui = ui, server = server)
