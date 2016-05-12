#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if(require("shiny") == FALSE){
	
	install.packages("shiny", dependencies = TRUE)
}

if(require("Biobase") == FALSE){
	source("https://bioconductor.org/biocLite.R")
	
	biocLite("Biobase")
}

if(require("Biostrings") == FALSE){
	source("https://bioconductor.org/biocLite.R")
	
	biocLite("Biostrings")
}

library(shiny)
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
            selectInput("seqType", 
					    label=h5("What is the sequence type of this file?"), 
                        choices = list("DNA" = 1, 
                                       "RNA" = 2, 
                                       "Protein" = 3),),
			
		   selectInput("convert", 
				   		label=h5("Sequence Conversion:"), 
				   		choices = list("DNA" = 1, 
						   			   "RNA" = 2, 
						    	   	   "Protein" = 3),),
					   
            helpText("Tip:"),
			helpText("You cannot convert a protein sequence to DNA or RNA.")
			
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            h4("To get started, choose a fasta file on the left."),
			uiOutput("seqText"),
            plotOutput("plotMM"),
			uiOutput("conversion")
        )
    )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
	
	check = 0
	
	##read in FASTA file
    fasta <- reactive({
        inFile <- input$file
        if(is.null(inFile) == F){
            
			if(input$seqType == 1){ ##DNA fasta file
				
				readDNAStringSet(inFile$name)
				
			} else if(input$seqType == 2){ ##RNA fasta file
				
				readRNAStringSet(inFile$name)
			} else if(input$seqType == 3){ ##Protein fasta file
				
				readAAStringSet(inFile$name)
			}
        }
        #readFASTA(file = inFile$datapath)
    })
	
	##generates a summary of the file
    output$sum <- renderDataTable({
        if(is.null(fasta()) == F){
			seqName <- names(fasta())
			start <- 1
			end <- length(fasta()[[1]])
			nAAs <- suppressWarnings(length(translate(fasta()[[1]])))
			
			df <-data.frame("seqName" = seqName, "start" = start, "end" = end, "nBases" = end, "Num of AAs" = nAAs)
			
        }
    })
	
	##print DNA, RNA, or protein sequences
	output$seqText <- renderText({
		
		if(is.null(fasta()) == F){
			
			if(input$seqType == 1){ ##DNA
				
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				
				total <- length(which(dnaBases == "T"))
				if(total > 0){
					
					paste(dnaBases)
				}
				
			} else if(input$seqType == 2){ ##RNA
				
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				rna <- RNAString(x)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				
				total <- length(which(rnaBases == "U"))
				if(total > 0){
					
					paste(rnaBases)
				}
				
			} else if(input$seqType == 3){ ##protein
			
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				peptide <- AAString(x)
				aa <- strsplit(as.character(peptide), split = "")[[1]]
				
				total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
				if(total != length(aa)){ 
					
					paste(aa)
				}
			}
		}
	})
	
	##generate a bar plot of the DNA, RNA, or protein sequences in the file
    output$plotMM <- renderPlot({
				
        if(is.null(fasta()) == F){
			
			if(input$seqType == 1){ ##DNA
				
		        f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				
				total <- length(which(dnaBases == "T"))
				if(total > 0){
					
		        	barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA composition of each nucleotide" )
				}
				
			} else if(input$seqType == 2){ ##RNA
				
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				rna <- RNAString(x)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				
				total <- length(which(rnaBases == "U"))
				if(total > 0){
					
					barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA composition of each nucleotide" )
				}
				
			} else if(input$seqType == 3){ ##protein
				
				f <- fasta()[[1]][1:length(fasta()[[1]])]
				x <- paste(f, collapse = "")
				peptide <- AAString(x)
				aa <- strsplit(as.character(peptide), split = "")[[1]]
				
				total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
				if(total != length(aa)){ 
					
					barplot(table(aa), xlab = "Amino acids", ylab = "Number of amino acids",  main = "Peptide composition of each amino acid" )
				}
			}
		}
    })

	##converts the sequences in the following manner
	##1) DNA to RNA
	##2) DNA to Protein
	##3) RNA to DNA
	##4) RNA to Protein
	output$conversion <- renderText({
			
		if(is.null(fasta()) == F){
			
			
			##1) check to see if seq type and convert are not the same
			##2) check to see if seq type is not protein
			if((input$convert != input$seqType) & input$seqType != 3){
				
				##convert seqtype to user's choice
				if((input$seqType == 2) & (input$convert == 1)){ ##convert RNA to DNA
					
					f <- fasta()[[1]][1:length(fasta()[[1]])]
					x <- paste(f, collapse = "")
					rna <- RNAString(x)
					dna <- DNAString(rna)
					paste("Converted DNA sequence: ", dna)
				} else if((input$seqType == 2) & (input$convert == 3)){ ##convert RNA to protein
				
					f <- fasta()[[1]][1:length(fasta()[[1]])]
					x <- paste(f, collapse = "")
					rna <- RNAString(x)
					peptide <- translate(rna)
					paste("Converted protein sequence: ", peptide)	
				} else if((input$seqType == 1) & (input$convert == 2)){ ##convert DNA to RNA
					
					f <- fasta()[[1]][1:length(fasta()[[1]])]
					x <- paste(f, collapse = "")
					dna <- DNAString(x)
					rna <- RNAString(dna)
					paste("Converted RNA sequence: ", rna)	
				} else if((input$seqType == 1) & (input$convert == 3)){ ##convert DNA to protein
					
					f <- fasta()[[1]][1:length(fasta()[[1]])]
					x <- paste(f, collapse = "")
					dna <- DNAString(x)
					peptide <- translate(dna)
					paste("Converted protein sequence: ", peptide)
				}
			}
		
		}
				
	})

})

# Run the application 
shinyApp(ui = ui, server = server)

