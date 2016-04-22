#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("Fasta File Reader"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        # sliderInput("bins",
        #             "Number of bins:",
        #             min = 1,
        #             max = 50,
        #             value = 30),
        fileInput("file", "Choose FASTA file",accept=c('fasta'))
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
     
     output$contents<-renderTable({
      inFile<-input$file
      print("hello!")
      if(is.null(inFile)){
        return (NULL)
      }
      read.csv(inFile$datapath,header=TRUE,sep=',',
               quote=input$quote)
      
   })
      
      
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

