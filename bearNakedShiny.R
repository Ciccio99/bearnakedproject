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
            h3("Get started:"),
            fileInput("file", "Choose FASTA file", accept=c('.fasta')),
            selectInput("seqType", 
					    label=h5("What is the sequence type of this file?"), 
                        choices = list("DNA" = 1, 
                                       "RNA" = 2, 
                                       "Protein" = 3),),
            helpText("You may not select DNA or RNA if the fasta file
                     contains a peptide sequence."),
            uiOutput('seqSelectInput')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            h5("File information:"),
            uiOutput("fileTab"),
            # h5("Fasta Summary:"),
            # uiOutput("sum"),
            h4("Sequence Name:"),
            textOutput("seqName"),
            h4("DNA/RNA Sequence:"),
			uiOutput("seqText"),
            plotOutput("plotMM")
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
        #read.fasta(file = inFile$datapath)
        readDNAStringSet(file=inFile$datapath, format="fasta")
    })

	##check to see if anything is in the file -- do not think this is necessary?
    output$fileTab <- renderTable({
        if(is.null(fasta())){
           return()
        }
        input$file
    })
	
	# ##generates a summary of the file
 #    output$sum <- renderTable({
 #        if(is.null(fasta())){
 #            return ()
 #        }
 #        summary(fasta()[1])
 #    })

    output$seqSelectInput <- renderUI({
        seqNums <- 1:length(fasta())
        selectInput("seqSelectNum", label=h5("Choose a sequence:"), choices=seqNums)
    })
	
    output$seqName <- renderText({
        if (is.null(fasta()) || is.null(input$seqSelectNum)){
            return ()
        }
        num <- as.numeric(input$seqSelectNum)
        names(fasta()[num])
    })
    

	##print DNA, RNA, or protein sequences
	output$seqText <- renderText({
		
		if(is.null(fasta()) == F){
	        seqNum <- as.numeric(input$seqSelectNum)
			if(input$seqType == 1){ ##DNA
				
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				paste(dnaBases)
				
			} else if(input$seqType == 2){ ##RNA
				
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				rna <- RNAString(dna)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				paste(rnaBases)
				
			} else if(input$seqType == 3){ ##protein
				
				
			}
		}
	})
	
	##generate a bar plot of the DNA, RNA, or protein sequences
    output$plotMM <- renderPlot({
        if(is.null(fasta()) == F){
            seqNum <- as.numeric(input$seqSelectNum)
			if(input$seqType == 1){ ##DNA
				
		        f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
		        barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA composition of each nucleotide" )
				
			} else if(input$seqType == 2){ ##RNA
				
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				rna <- RNAString(dna)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA composition of each nucleotide" )
				
			} else if(input$seqType == 3){ ##protein
				
				
			}
		}
    })

})

# Run the application 
shinyApp(ui = ui, server = server)

