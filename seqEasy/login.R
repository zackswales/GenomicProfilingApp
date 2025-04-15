library(shiny)
library(shinyjs)

# Define the path to your CSV file
credentials_file <- '~/GenomicProfilingApp/register.csv'

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
           tags$h3("Welcome to seqEasy"),
           tags$p(textOutput("logged_in_user")), # Display logged-in user
           actionButton("logout_button", "Log out") # Add logout button
  )
)

login_dialog <- modalDialog(
  title = 'Login to continue',
  footer = tagList(
    actionButton('tab_login.login', 'Login'),
    actionButton('tab_login.register', 'Register'),
    actionButton('tab_login.delete_account', 'Delete Account') # Delete account button
  ),
  textInput('tab_login.username', 'Username'),
  passwordInput('tab_login.password', 'Password'),
  tags$div(class = 'warn-text', textOutput('tab_login.login_msg'))
)

register_dialog <- modalDialog(
  title = 'Register',
  footer = actionButton('tab_login.create_account', 'Create Account'),
  textInput('tab_login.new_username', 'New Username'),
  passwordInput('tab_login.new_password', 'New Password'),
  tags$div(class = 'warn-text', textOutput('tab_login.register_msg'))
)

delete_dialog <- modalDialog(
  title = 'Delete Account',
  footer = actionButton('tab_login.confirm_delete', 'Delete'),
  textInput('tab_login.delete_username', 'Username'),
  passwordInput('tab_login.delete_password', 'Password'),
  tags$div(class = 'warn-text', textOutput('tab_login.delete_msg'))
)

tab_login$server <- function(input, output, session) {
  
  # Show login dialog box when initiated
  showModal(login_dialog)
  
  observeEvent(input$tab_login.login, {
    # Read the CSV every time login is attempted
    user.access <- read.csv(credentials_file)
    
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
        
        # Display logged-in user
        output$logged_in_user <- renderText(paste("Logged in as:", username)) 
        
        # Show logout button
        shinyjs::show("logout_button")
        
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
  
  observeEvent(input$tab_login.register, {
    showModal(register_dialog)
  })
  
  observeEvent(input$tab_login.create_account, {
    new_username <- input$tab_login.new_username
    new_password <- input$tab_login.new_password
    
    if (nchar(new_username) > 0 && nchar(new_password) > 0) {
      
      # Read the CSV before creating a new account
      user.access <- read.csv(credentials_file)
      
      if (new_username %in% user.access$Name) {
        output$tab_login.register_msg <- renderText('Username already exists.')
        shinyjs::show('tab_login.register_msg')
        shinyjs::delay(1000, shinyjs::hide('tab_login.register_msg'))
      } else {
        new_row <- data.frame(Name = new_username, Password = new_password)
        user.access <- rbind(user.access, new_row)  # Update user.access in the server scope
        write.csv(user.access, credentials_file, row.names = FALSE)  # Write back to CSV
        
        # Re-read the CSV after creating a new account
        user.access <- read.csv(credentials_file)
        
        removeModal()
        showModal(login_dialog)
      }
    } else {
      output$tab_login.register_msg <- renderText('Please enter a username and password.')
      shinyjs::show('tab_login.register_msg')
      shinyjs::delay(1000, shinyjs::hide('tab_login.register_msg'))
    }
  })
  
  observeEvent(input$tab_login.delete_account, {
    showModal(delete_dialog)
  })
  
  observeEvent(input$tab_login.confirm_delete, {
    delete_username <- input$tab_login.delete_username
    delete_password <- input$tab_login.delete_password
    
    if (nchar(delete_username) > 0 && nchar(delete_password) > 0) {
      
      # Read the CSV before deleting an account
      user.access <- read.csv(credentials_file)
      row.names(user.access) <- user.access$Name
      
      if (delete_username %in% user.access[, 'Name']) {
        if (delete_password == user.access[delete_username, 'Password']) {
          
          # Delete user directory and contents
          user_dir <- file.path(base_dir, delete_username)
          if (dir.exists(user_dir)) {
            unlink(user_dir, recursive = TRUE)
          }
          
          # Remove user from CSV
          user.access <- user.access[user.access$Name != delete_username, ]
          write.csv(user.access, credentials_file, row.names = FALSE)
          
          # Re-read the CSV after deleting an account
          user.access <- read.csv(credentials_file)
          
          removeModal()
          showModal(login_dialog)
        } else {
          output$tab_login.delete_msg <- renderText('Incorrect Password')
          shinyjs::show('tab_login.delete_msg')
          shinyjs::delay(1000, shinyjs::hide('tab_login.delete_msg'))
        }
      } else {
        output$tab_login.delete_msg <- renderText('Username Not Found')
        shinyjs::show('tab_login.delete_msg')
        shinyjs::delay(1000, shinyjs::hide('tab_login.delete_msg'))
      }
    } else {
      output$tab_login.delete_msg <- renderText('Please enter a username and password.')
      shinyjs::show('tab_login.delete_msg')
      shinyjs::delay(1000, shinyjs::hide('tab_login.delete_msg'))
    }
  })
  
  observeEvent(input$logout_button, {
    # Reset logged-in user text
    output$logged_in_user <- renderText("")
    
    # Hide logout button
    shinyjs::hide("logout_button")
    
    output$ggplotheatmap <- renderPlot({
      NULL # Render an empty plot
    })
    
    # Show login dialog
    showModal(login_dialog)
  })
  
  # Hide logout button initially
  shinyjs::hide("logout_button")
}