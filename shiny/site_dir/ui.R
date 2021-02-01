library(shiny)
library(shinyjs)

# set up the menu choices for the plot selector
# first 30 are replicates of main t-SNE
plot_indices = 1:30
plot_names = unlist(lapply(plot_indices, function(x) paste("rep",x)))
plots = plot_indices
names(plots) = plot_names
# last one is the S-only map
plots = append(plots, c("S only (No An. coluzzii)"=31))
# the perplexity variants
plots = append(plots, c("t-SNE perplexity 50"=32))
plots = append(plots, c("t-SNE perplexity 100"=33))
plots = append(plots, c("t-SNE perplexity 250"=34))
plots = append(plots, c("t-SNE perplexity 1000"=35))

# and now the PCA plots
plots = append(plots, c("PCA PC1 vs PC2"=36))
plots = append(plots, c("PCA PC1 vs PC3"=37))
plots = append(plots, c("PCA PC2 vs PC3"=38))
plots = append(plots, c("PCA PC3 vs PC4"=39))

ui <- fluidPage(

  useShinyjs(),  # Include shinyjs
  
  # Some custom CSS for a smaller font for preformatted text
  tags$head(
    tags$style(HTML("
                    pre, table.table {
                    font-size: 14px;
                    }
                    "))
  ),

  fluidRow(
    h3("All-gene t-SNE plot"),
    column(width = 3, 
      wellPanel(
        selectizeInput(inputId = "plot_number",
	               label = "t-SNE plot",
		       choices = plots
                      )		       
      ),

      wellPanel(
        selectizeInput(inputId = "plot_type",
                       label = "Plot colouring",
                       choices = list(
                         `Gene info` = c("Chromosome and inversions" = "chr.factor", 
                                         "log(#SNPs)" = "log.snps",
                                         "#SNPs (rank lo-hi)" = "snp.rank",
                                         "log(1+#non-synonymous)" = "log.non.snps",
                                         "%non-synonymous SNPs" = "percent.nonsyn"
                                        ),
                         `Chromosomal location` = c("2L centro->telo" = "gene.start_2L",
			                     "2R telo->centro" = "gene.start_2R",
					     "3L centro->telo" = "gene.start_3L",
					     "3R telo->centro" = "gene.start_3R",
					     "X telo->centro" = "gene.start_X"),

                         `Consistency between 30 replicate plots` = c("mean variance" = "mean_variance",
			 	   	         "median variance" = "median_variance",
						 "# unique nearest 5 neighbours" = "nun5",
						 "# unique nearest 20 neighbours" = "nun20",
						 "# unique nearest 50 neighbours" = "nun50"
                                                ),

                         `Fst` = c("General" = "fst.general",
                                   "General (rank hi-lo)" = "fst.general.rank",
                                   "East-West" = "fst.east_west",
                                   "M-S" = "fst.m_s",
                                   "M-S (rank hi-lo)" = "fst.m_s.rank",
				   "2La" = "fst.2l",
				   "2Rb" = "fst.2r"),

                         `Fst population vs rest` =
                           c("Angola" = "aom", "Burkina Faso M" = "bfm", "Burkina Faso S" = "bfs",
                             "Cameroon" = "cms", "Gabon" = "gas", "Guinea" = "gns",
                             "Guinea-Bissau" = "gwa", "Kenya" = "kes", "Uganda" = "ugs"),
                         `Expression max` =
                           c("Tissues (Baker)" = "tissues_max",
                             "Blood meal (Marinotti)" = "bloodmeal_max",
                             "Circadian heads L-D (Rund)" = "circadian_max",
                             "Embryonic devel. (Goltsev)" = "embryonic_max",
                             "Dessication stress (Wang)" = "dessication_max"),
                         `Expression range` =
                           c("Tissues (Baker)" = "tissues_range",
                             "Blood meal (Marinotti)" = "bloodmeal_range",
                             "Circadian heads L-D (Rund)" = "circadian_range",
                             "Embryonic devel. (Goltsev)" = "embryonic_range",
                             "Dessication stress (Wang)" = "dessication_range")
                       )
        ) 
      ),
      wellPanel(
        textAreaInput(inputId = "gene_id_search",
                      label = "Highlight genes by ID",
                      height = "160px",
		      placeholder = "AGAP000001, AGAP000002\n(any delimiter will do)"),
        
        fluidRow(
          column(width=9, verbatimTextOutput("gene_id_search_result")),
          column(width=2, shinyjs::hidden(actionButton("reset_search", "Clear")))
        )
     )     
    ),

    
    column(width = 6,
           plotOutput("plot1", height = 500,
                      # Equivalent to: click = clickOpts(id = "plot_click")
                      click = "plot_click",
                      dblclick = "plot1_dblclick",
                      brush = brushOpts(
                        id = "plot1_brush",
#                        delay = 500,
#                        delayType = "throttle",
                        resetOnNew = FALSE)
           )),
    column(width = 3,
           h4("Help"),
           tags$ul(tags$li("Gene info: Click on a point to see more details below."),
                   tags$li("Gene list: Click/drag a rectangle to see IDs below."),
                   tags$li("Zoom in: Click/drag a rectangle and double click it."),
                   tags$li("Clear rectangle: Single click outside it."),
                   tags$li("Reset zoom: Double click in plot area.")),
           tags$br(), # spacer
           shinyjs::hidden(h4(id="karyogram-h4", "Genomic location of highlighted genes")),
           plotOutput("karyo", height = 300)
    )
  ),

  fluidRow(
    column(width = 6,
        h3("Gene details"),
        verbatimTextOutput("click_info")
    ),
    
    
#    column(width = 5,
#           h3("Gene-wise t-SNE plot"),
#           # In a plotOutput, passing values for click, dblclick, hover, or brush
#           # will enable those interactions.
#           # plotOutput("click_plot", height = 350)
#          tags$p("TBC")
#    )

    column(width = 6,
      h3("Genes in selection"),
      shinyjs::hidden(actionButton("highlight_genes", "Highlight these genes")),
      verbatimTextOutput("brush_genes"),
      h3("Gene neighbours"),
      numericInput("NN", "Number of nearest neighbours:", 20, min = 1, max = 500),        
      verbatimTextOutput("neighbours")
    )
  )    
)

