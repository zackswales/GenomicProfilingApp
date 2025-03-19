user.access <- read.csv('~/Research Project/register.csv')

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
  tags$div(class = 'warn-text',textOutput('tab_login.login_msg'))
)

tab_login$server <- function(input, output) {
  
  # show login dialog box when initiated
  showModal(login_dialog)
  
  observeEvent(input$tab_login.login, {
    username <- input$tab_login.username
    password <- input$tab_login.password
    row.names(user.access) <- user.access$Name
    
    # validate login credentials
    if(username %in% user.access[,'Name']) {
      if(password == user.access[username,'Password']) {
        # succesfully log in
        removeModal() # remove login dialog
        output$tab_login.welcome_text <- renderText(glue('welcome, {username}'))
        shinyjs::show('tab_login.welcome_div') # show welcome message
      } else {
        # password incorrect, show the warning message
        # warning message disappear in 1 sec
        output$tab_login.login_msg <- renderText('Incorrect Password')
        shinyjs::show('tab_login.login_msg')
        shinyjs::delay(1000, hide('tab_login.login_msg'))
      }
    } else {
      # username not found, show the warning message
      # warning message disappear in 1 sec
      output$tab_login.login_msg <- renderText('Username Not Found')
      shinyjs::show('tab_login.login_msg')
      shinyjs::delay(1000, hide('tab_login.login_msg'))
    }
  })
}