#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if(require("Biobase") == FALSE){
	source("https://bioconductor.org/biocLite.R")
	
	biocLite("Biobase")
}

if(require("Biostrings") == FALSE){
	source("https://bioconductor.org/biocLite.R")
	
	biocLite("Biostrings")
}

library(shiny)
library(seqinr)
library(Biobase)
library(Biostrings)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    # Application title
    titlePanel("Welcome to FASTAcrobat!"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Choose FASTA file", accept=c('.fasta')),
            selectInput("seqtype", 
					    label=h5("What is the sequence type of this file?"), 
                        choices = list("DNA" = 1, 
                                       "RNA" = 2, 
                                       "Protein" = 3),),
					   
            helpText("You may not select DNA or RNA if the fasta file
                     contains a peptide sequence.")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            h4("To get started, choose a fasta file on the left."),
            uiOutput("fileTab"),
            uiOutput("sum"),
			uiOutput("DNAtext"),
            plotOutput("distPlot"),
			uiOutput("RNAtext"),
			plotOutput("convertDNA_RNA")
        )
    )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
	
	##read in FASTA file
    fasta <- reactive({
        inFile <- input$file
        if(is.null(inFile)){
            return ()
        }
        read.fasta(file = inFile$datapath)
    })

	##check to see if anything is in the file -- do not think this is necessary?
    output$fileTab <- renderTable({
        if(is.null(fasta())){
           return()
        }
        input$file
    })
	
	##generates a summary of the file
    output$sum <- renderTable({
        if(is.null(fasta())){
            return ()
        }
        summary(fasta())
    })
	
	##print RNA sequences
	output$DNAtext <- renderText({
		
		if(is.null(fasta()) == F){
			
			f <- fasta()[[1]][1:length(fasta()[[1]])]
			x <- paste(f, collapse = "")
			dna <- DNAString(x)
			dnaBases <- strsplit(as.character(dna), split = "")[[1]]
			paste(dnaBases)
		}
	})
	
	##generate a bar plot of the DNA sequences
    output$distPlot <- renderPlot({
        if(is.null(fasta()) == F){
           
	        f <- fasta()[[1]][1:length(fasta()[[1]])]
			x <- paste(f, collapse = "")
			dna <- DNAString(x)
			dnaBases <- strsplit(as.character(dna), split = "")[[1]]
			
	        barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA Composition of each nucleotide" )
		}
    })

	##print RNA sequences
	output$RNAtext <- renderText({
			
			if(is.null(fasta()) == F){
				
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				rna <- RNAString(dna)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				paste(rnaBases)
			}
		})
	
	##generate a bar plot of the RNA sequences
	output$convertDNA_RNA <- renderPlot({
			if(is.null(fasta()) == F){
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				rna <- RNAString(dna)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA Composition of each nucleotide" )
			}
	})
})

# Run the application 
shinyApp(ui = ui, server = server)

