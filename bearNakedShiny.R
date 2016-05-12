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
    titlePanel("FASTAcrobat"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h3("Get started:"),
            fileInput("file", "Choose FASTA file", accept=c('.fasta')),
            selectInput("seqType", 
					    label=h4("What is the sequence type of this file?"), 
                        choices = list("Select Sequence Type" = 0,
                                        "DNA" = 1, 
                                       "RNA" = 2, 
                                       "Protein" = 3)),
            uiOutput('seqSelectInput')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            h4("Sequence Name:"),
            textOutput("seqName"),
            hr(),
            tabsetPanel(type="tabs",
                tabPanel("DNA",
                    p(textOutput("dnaSequence"))
                ),
                tabPanel("RNA",
                    p(textOutput("rnaSequence"))
                ),
                tabPanel("Protein",
                    p(textOutput("proteinSequence"))
                ),
                tabPanel("DNA Composition",
                    plotOutput("plotDNA")
                ),
				tabPanel("RNA Composition",
						plotOutput("plotRNA")
				),
				tabPanel("Protein Composition",
                        helpText("Fun Tip:"),
                        helpText('"X" Amino Acid = Ambiguous'),
						plotOutput("plotProtein")
				)
            )
        )
    )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {	
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

    # Instantiates sequence choices once the file is loaded
    output$seqSelectInput <- renderUI({
        if (is.null(fasta())) return ()

        seqNums <- 1:length(fasta())
        selectInput("seqSelectNum", label=h4("Choose a sequence:"), choices=seqNums)
    })

    # Sets name of the sequence
    output$seqName <- renderText({
        if (is.null(fasta()) || is.null(input$seqSelectNum)) return ()

        num <- as.numeric(input$seqSelectNum)
        names(fasta()[num])
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
	seqText <- renderText({
		if(is.null(fasta()) == F){
			seqNum <- as.numeric(input$seqSelectNum)

			if(input$seqType == 1){ ##DNA
				
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				dna <- DNAString(x)
				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
				
				total <- length(which(dnaBases == "T"))
				if(total > 0){
					
					paste(dnaBases)
				}
				
			} else if(input$seqType == 2){ ##RNA
				
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
				x <- paste(f, collapse = "")
				rna <- RNAString(x)
				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
				
				total <- length(which(rnaBases == "U"))
				if(total > 0){
					
					paste(rnaBases)
				}
				
			} else if(input$seqType == 3){ ##protein
			
				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
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
	output$plotDNA <- renderPlot({
		if (is.null(fasta())) return ()

		seqNum <- as.numeric(input$seqSelectNum)
		if (input$seqType == 1) { ##DNA	
	        f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
			x <- paste(f, collapse = "")
			dna <- DNAString(x)
			dnaBases <- strsplit(as.character(dna), split = "")[[1]]
			
			total <- length(which(dnaBases == "T"))
			if(total > 0){
				
	        	barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA composition of each nucleotide" )
			}
		} else if (input$seqType == 2) { ##DNA    
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            rna <- RNAString(x)
            dnaBases <- strsplit(as.character(DNAString(rna)), split = "")[[1]]
            
            total <- length(which(dnaBases == "T"))
            if(total > 0){
                
                barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA composition of each nucleotide" )
            }
        }
	})

	output$plotRNA <- renderPlot({
		if (is.null(fasta())) return ()

		seqNum <- as.numeric(input$seqSelectNum)
        if (input$seqType == 1) {
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            dna <- DNAString(x)
            rnaBases <- strsplit(as.character(RNAString(dna)), split = "")[[1]]
            total <- length(which(rnaBases == "U"))
            if (total > 0) {
                barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA composition of each nucleotide" )
            }
        } else if (input$seqType == 2) { ##RNA
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            rna <- RNAString(x)
            rnaBases <- strsplit(as.character(rna), split = "")[[1]]
            total <- length(which(rnaBases == "U"))
            total <- length(which(rnaBases == "U"))
            if (total > 0) {
                barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA composition of each nucleotide" )
            }
		}
	})

	output$plotProtein <- renderPlot({
		if(is.null(fasta())) return ()

		seqNum <- as.numeric(input$seqSelectNum)
        if(input$seqType == 1){ ##protein   
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            peptide <- AAString(translate(DNAString(x),  if.fuzzy.codon="solve"))
            aa <- strsplit(as.character(peptide), split = "")[[1]]
            
            total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
            if(total != length(aa)){ 
                barplot(table(aa), col="aquamarine2", xlab = "Amino acids", ylab = "Number of amino acids",  main = "Peptide composition of each amino acid" )
            }
        } else if(input$seqType == 2){ ##protein   
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            peptide <- AAString(translate(RNAString(x), if.fuzzy.codon="solve"))
            aa <- strsplit(as.character(peptide), split = "")[[1]]
            
            total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
            if(total != length(aa)){ 
                barplot(table(aa), col="aquamarine2", xlab = "Amino acids", ylab = "Number of amino acids",  main = "Peptide composition of each amino acid" )
            }
        } else if(input$seqType == 3){ ##protein	
			f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
			x <- paste(f, collapse = "")
			peptide <- AAString(x)
			aa <- strsplit(as.character(peptide), split = "")[[1]]
			
			total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
			if(total != length(aa)){ 
				barplot(table(aa), col="aquamarine2", xlab = "Amino acids", ylab = "Number of amino acids",  main = "Peptide composition of each amino acid" )
			}
		}
	})

#    output$plotMM <- renderPlot({
#				
#        if(is.null(fasta()) == F){
#            seqNum <- as.numeric(input$seqSelectNum)
#			if(input$seqType == 1){ ##DNA
#				
#		        f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
#				x <- paste(f, collapse = "")
#				dna <- DNAString(x)
#				dnaBases <- strsplit(as.character(dna), split = "")[[1]]
#				
#				total <- length(which(dnaBases == "T"))
#				if(total > 0){
#					
#		        	barplot(table(dnaBases), xlab = "Base", ylab = "Number of bases",  main = "DNA composition of each nucleotide" )
#				}
#				
#			} else if(input$seqType == 2){ ##RNA
#				
#				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
#				x <- paste(f, collapse = "")
#				rna <- RNAString(x)
#				rnaBases <- strsplit(as.character(rna), split = "")[[1]]
#				
#				total <- length(which(rnaBases == "U"))
#				if(total > 0){
#					
#					barplot(table(rnaBases), xlab = "Base", ylab = "Number of bases",  main = "RNA composition of each nucleotide" )
#				}
#				
#			} else if(input$seqType == 3){ ##protein
#				
#				f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
#				x <- paste(f, collapse = "")
#				peptide <- AAString(x)
#				aa <- strsplit(as.character(peptide), split = "")[[1]]
#				
#				total <- length(which(aa == "A")) + length(which(aa == "T")) + length(which(aa == "G")) + length(which(aa == "C")) + length(which(aa == "U"))
#				if(total != length(aa)){ 
#					
#					barplot(table(aa), xlab = "Amino acids", ylab = "Number of amino acids",  main = "Peptide composition of each amino acid" )
#				}
#			}
#		}
#    })

    # Instantiates DNA sequence
    output$dnaSequence <- renderText ({
        if (is.null(fasta())) return()
        seqNum <- as.numeric(input$seqSelectNum)
        if (input$seqType == 2) {
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            rna <- RNAString(x)
            strsplit(as.character(DNAString(rna)), split = "")[[1]]
        } else if (input$seqType == 1) {
            seqText()
        }   
    })

    # Instantiates RNA sequence 
    output$rnaSequence <- renderText ({
        if (is.null(fasta())) return()
        seqNum <- as.numeric(input$seqSelectNum)
        if (input$seqType == 1) {
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            dna <- DNAString(x)
            strsplit(as.character(RNAString(dna)), split = "")[[1]]
        } else if (input$seqType == 2) {
            seqText()
        }   
    })

    # Converts DNA/RNA to protein and Instantiates it
    output$proteinSequence <- renderText ({
        if (is.null(fasta())) return()

        seqNum <- as.numeric(input$seqSelectNum)
        if (input$seqType == 1) {
            #If DNA -> Protein
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            dna <- DNAString(x)
            strsplit(as.character(translate(dna)), split = "")[[1]]
        } else if (input$seqType == 2) {
            # If RNA -> Protein
            f <- fasta()[[seqNum]][1:length(fasta()[[seqNum]])]
            x <- paste(f, collapse = "")
            rna <- RNAString(x)
            strsplit(as.character(translate(rna)), split = "")[[1]]
        } else if(input$seqType == 3){
			seqText()
		}
    })
})

# Run the application 
shinyApp(ui = ui, server = server)

