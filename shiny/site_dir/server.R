library(shiny)
library(shinyjs)
library(ggplot2)

tsne_data_file <- "data/bigfile.txt"
tsne_data <- read.csv(tsne_data_file, sep="\t", stringsAsFactors = F, header=TRUE)
gene_ID_colname <- "gene_name"
row.names(tsne_data) <- tsne_gene_IDs <- tsne_data[,gene_ID_colname]

#
# testing adding ranked measures
#
tsne_data$fst.m_s.rank = order(order(-tsne_data$fst.m_s))
tsne_data$fst.general.rank = order(order(-tsne_data$fst.general))
tsne_data$snp.rank = order(order(tsne_data$snp))

# now load in all 30 tnse coordinates
tsne30_file <- "data/super30tsne.txt"
tsne30_data <- read.table(tsne30_file)
# and copy the columns over to tsne_data
for (column in colnames(tsne30_data)) {
  tsne_data[[column]] <- tsne30_data[,column]
}
rm(tsne30_data)
# colnames(tsne_data)
# remove these because we just loaded them in again
tsne_data$tsne1 <- NULL
tsne_data$tsne2 <- NULL

# now load in the S-only map coordinates
S_only_file <- "data/output_S_bob.txt"
S_only_data <- read.csv(S_only_file, sep="\t", stringsAsFactors = F, header=TRUE, row.names=1)
tsne_data$tsne31_1 <- S_only_data[tsne_gene_IDs,"tsne1_S"]
tsne_data$tsne31_2 <- S_only_data[tsne_gene_IDs,"tsne2_S"]
#tsne_data[["tsne31_1"]] <- S_only_data[,"tsne1_S"]
#tsne_data[["tsne31_2"]] <- S_only_data[,"tsne2_S"]

# now load the different perplexity t-SNE plots
perplexity_file <- "data/tsne-perplexity.txt"
perplexity_data <- read.table(perplexity_file)
tsne_data$tsne32_1 <- perplexity_data[tsne_gene_IDs, "tsne_p50_1"]
tsne_data$tsne32_2 <- perplexity_data[tsne_gene_IDs, "tsne_p50_2"]
tsne_data$tsne33_1 <- perplexity_data[tsne_gene_IDs, "tsne_p100_1"]
tsne_data$tsne33_2 <- perplexity_data[tsne_gene_IDs, "tsne_p100_2"]
tsne_data$tsne34_1 <- perplexity_data[tsne_gene_IDs, "tsne_p250_1"]
tsne_data$tsne34_2 <- perplexity_data[tsne_gene_IDs, "tsne_p250_2"]
tsne_data$tsne35_1 <- perplexity_data[tsne_gene_IDs, "tsne_p1000_1"]
tsne_data$tsne35_2 <- perplexity_data[tsne_gene_IDs, "tsne_p1000_2"]
rm(perplexity_data)


# and now the PCA plots
pcs_file = "data/pcs1-4.txt"
pcs_data = read.table(pcs_file)
# PC1 vs PC2
tsne_data$tsne36_1 <- pcs_data$PC1
tsne_data$tsne36_2 <- pcs_data$PC2
# PC1 vs PC3
tsne_data$tsne37_1 <- pcs_data$PC1
tsne_data$tsne37_2 <- pcs_data$PC3
# PC2 vs PC3
tsne_data$tsne38_1 <- pcs_data$PC2
tsne_data$tsne38_2 <- pcs_data$PC3
# PC3 vs PC4
tsne_data$tsne39_1 <- pcs_data$PC3
tsne_data$tsne39_2 <- pcs_data$PC4





chr_cols <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","gold","#80b1d3",
          "#999999","#a65628","#f781bf")
tsne_data$chr.factor <- as.factor(tsne_data$chro3)
log.snps <- log(tsne_data$snp)
log.non.snps <- log(tsne_data$non + 1)
percent.nonsyn <- 100.0 * tsne_data$non / (tsne_data$non + tsne_data$syn)

gene.start_2L <- tsne_data$gene_start
gene.start_2L[tsne_data$chro != "2L"] <- NA
gene.start_2R <- tsne_data$gene_start
gene.start_2R[tsne_data$chro != "2R"] <- NA
gene.start_3L <- tsne_data$gene_start
gene.start_3L[tsne_data$chro != "3L"] <- NA
gene.start_3R <- tsne_data$gene_start
gene.start_3R[tsne_data$chro != "3R"] <- NA
gene.start_X <- tsne_data$gene_start
gene.start_X[tsne_data$chro != "X"] <- NA

gene_data_file <- "data/anopheles_gambiae_genes_VB_search.csv"
gene_data <- read.csv(gene_data_file, stringsAsFactors=FALSE, header=TRUE)
row.names(gene_data) <- gene_data[,"accession"]

# get the lengths of all the chromosomes
chr_summary <- aggregate(seq_region_end ~ seq_region_name, data=gene_data, max)
# ignore Mt and UNKN 'chromosomes'
chr_summary <- subset(chr_summary, !seq_region_name %in% c("Mt", "UNKN"))

chr_summary$chr_index <- 1:nrow(chr_summary)
# rename max aggregated seq_region_end column as chr_size
names(chr_summary)[names(chr_summary)=="seq_region_end"] <- "chr_size"

# copy the chromosome/inversion factor into gene_data for colouring karyograms
gene_data[tsne_gene_IDs,"chr.factor"] = tsne_data[tsne_gene_IDs,"chr.factor"]

expression_ranges_file = 'data/expression-ranges.tsv'
expression_ranges = read.table(expression_ranges_file, sep="\t")

# copy the data over into the main data frame, using the gene IDs
# of the t-SNE
tsne_data$circadian_range = expression_ranges[tsne_gene_IDs,"circadian"]
tsne_data$tissues_range = expression_ranges[tsne_gene_IDs,"tissues"]
tsne_data$embryonic_range = expression_ranges[tsne_gene_IDs,"embryonic"]
tsne_data$bloodmeal_range = expression_ranges[tsne_gene_IDs,"bloodmeal"]
tsne_data$dessication_range = expression_ranges[tsne_gene_IDs,"dessication"]

expression_max_file = 'data/expression-maximums.tsv'
expression_max = read.table(expression_max_file, sep="\t")

# copy the data over into the main data frame, using the gene IDs
# of the t-SNE
tsne_data$circadian_max = expression_max[tsne_gene_IDs,"circadian"]
tsne_data$tissues_max = expression_max[tsne_gene_IDs,"tissues"]
tsne_data$embryonic_max = expression_max[tsne_gene_IDs,"embryonic"]
tsne_data$bloodmeal_max = expression_max[tsne_gene_IDs,"bloodmeal"]
tsne_data$dessication_max = expression_max[tsne_gene_IDs,"dessication"]

distance_data_file = "data/gene-stats-local-w-headers.txt"
distance_data = read.table(distance_data_file)
tsne_data$mean_variance = distance_data$mean_variance
tsne_data$median_variance = distance_data$median_variance
tsne_data$nun5 = distance_data$nun5
tsne_data$nun20 = distance_data$nun20
tsne_data$nun50 = distance_data$nun50
rm(distance_data)


myNearPoints <- function(df, queryCoord, xvar = NULL, yvar = NULL, maxpoints = 20) {
  sq_dists = (df[[xvar]]-queryCoord$x)^2 + (df[[yvar]]-queryCoord$y)^2

  sort_idx = order(sq_dists)
  sort_idx = sort_idx[seq_len(maxpoints)]

  df[sort_idx, ]
}

server <- (function(input, output, session) {

  genes <- reactiveValues(data = NULL)
 
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot1 <- renderPlot({

    # set up special colour palettes if needed
    if (input$plot_type == "chr.factor") {
      scale = scale_colour_manual(name="Key", values=chr_cols)
    } else {
      scale = scale_colour_gradient(name="Key", low="yellow", high="blue")
    }

    tsne_x = sprintf("tsne%02d_1", as.numeric(input$plot_number))
    tsne_y = sprintf("tsne%02d_2", as.numeric(input$plot_number))

    if (!is.null(input$plot_type)) {
      gg <- ggplot(tsne_data, aes(get(tsne_x), get(tsne_y))) +
                   geom_point(aes(colour=get(input$plot_type), alpha=1.0-is.na(get(input$plot_type))*0.9), size=2) +
                   labs(x=format("t-SNE1", digits=2), y=format("t-SNE2", digits=2)) +
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                   coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE) + scale + scale_alpha(guide="none")
    }

    if (!is.null(input$gene_id_search) && input$gene_id_search != '') {
      gene_ids = unlist(strsplit(input$gene_id_search, "\\W+", perl = TRUE))
      hits <- subset(tsne_data, row.names(tsne_data) %in% gene_ids)
      output$gene_id_search_result <- renderText({ sprintf("%d genes highlighted", nrow(hits)) })
      gg <- gg + geom_point(data=hits, size=3, colour="#101010", shape=22, stroke=2)
      shinyjs::show("reset_search")

      gene_data_hits <- subset(gene_data, row.names(gene_data) %in% gene_ids)
      gene_data_hits <- merge(gene_data_hits, chr_summary, by = "seq_region_name", all = FALSE)

      # figure out which chromosome/inversion colours we need

      chr_all = levels(gene_data$chr.factor)
      chr_needed = unique(gene_data_hits$chr.factor)
      chr_cols_needed <- chr_cols[which(chr_all %in% chr_needed)]

      karyo_gg = ggplot() +
                 geom_segment(data = chr_summary,
                              aes(x = seq_region_name, xend = seq_region_name, y = 0, yend = chr_size),
                              lineend = "round", color = "lightgrey", size = 5) +
                 geom_segment(data = gene_data_hits,
                              aes(x = chr_index - 0.1, xend = chr_index + 0.1,
                              y = seq_region_start, yend = seq_region_end, color = chr.factor),
                              size = 1) +
                 scale_colour_manual(name="Key", values=chr_cols_needed) +
                 geom_text(data = chr_summary,
                           aes(x = seq_region_name, y = -3000000, label=seq_region_name)) +
                              theme_void()

      output$karyo <- renderPlot(karyo_gg)
      shinyjs::show("karyogram-h4")
    } else {
      output$gene_search_num_hits <- NULL
      output$karyo <- NULL
      shinyjs::hide("karyogram-h4")
    }

    gg   
  })

  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    # shinyjs::logjs("observed plot1_dblclick")

    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      shinyjs::reset("plot1_brush")      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  

  output$click_info <- renderPrint({
    if (is.null(input$plot_click)){
      cat("Click on a point on the plot to show gene info...","\n")

      output$neighbours <- renderPrint({
        cat("For singly clicked genes only")
      })
    } else {
       tsne_x = sprintf("tsne%02d_1", as.numeric(input$plot_number))
       tsne_y = sprintf("tsne%02d_2", as.numeric(input$plot_number))

       location <- nearPoints(tsne_data, input$plot_click, xvar=tsne_x, yvar=tsne_y, addDist = FALSE)
       if (nrow(location) == 0){
       	   cat("Pick a gene by selecting a point on the plot","\n")
           output$neighbours <- renderPrint({
             cat("For singly clicked genes only")
           })
       } else {		
         gene_id = location[1,gene_ID_colname]
	 cat("Annotation:\n")
	 for (column in colnames(gene_data)) {
	   cat(sprintf("%-20s : %s\n", column, gene_data[gene_id,column]))
         }


         # all data from the data file
	 cat("\n\nData:\n")
         for (column in colnames(location)) {
	   # don't show the other 29 tsne coordinates:
	   if (column == tsne_x | column == tsne_y | !length(grep("tsne[0-3][0-9]_[12]", column, perl=TRUE))) {
  	     cat(sprintf("%-20s : %s\n", column, location[1,column]))
           }
         }

         output$neighbours <- renderPrint({
           cat("# number of replicate t-SNEs that gene IDs appear in nearest neighbours\n\n")
           counts <- new.env(hash=TRUE)

           countFun <- function(id) {
             if (is.null(counts[[id]])) { counts[[id]] = 0 }
             counts[[id]] = counts[[id]] + 1
           }

           for (plot_number in 1:30) {
             tsne_x = sprintf("tsne%02d_1", plot_number)
             tsne_y = sprintf("tsne%02d_2", plot_number)

             geneLoc = { }
             geneLoc$x = tsne_data[gene_id, tsne_x]
             geneLoc$y = tsne_data[gene_id, tsne_y]

             locs <- myNearPoints(tsne_data, geneLoc, xvar=tsne_x, yvar=tsne_y, maxpoints=input$NN)

             apply(locs[gene_ID_colname], 1, countFun)
           }
           # make a sortable data.frame with the hash count results
           countData = data.frame() 
           for (id in names(counts)) {
             countData[id,"count"] = counts[[id]]
             countData[id,"id"] = id
           }
 
           countSortIdx = order(countData$count, decreasing=TRUE)
           for (idx in countSortIdx) {
             cat(sprintf("%-2d %s\n", countData[idx,"count"], countData[idx,"id"]))
           }
         })

       }        
    }
  })
  
  output$brush_genes <- renderPrint({
    tsne_x = sprintf("tsne%02d_1", as.numeric(input$plot_number))
    tsne_y = sprintf("tsne%02d_2", as.numeric(input$plot_number))
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      # shinyjs::logjs("brush is not null")
      hits <- subset(tsne_data, get(tsne_x) < brush$xmax & get(tsne_x) > brush$xmin &
                                get(tsne_y) < brush$ymax & get(tsne_y) > brush$ymin);
      genes$brushed = row.names(hits)
      cat(genes$brushed, fill=TRUE)
      shinyjs::show("highlight_genes")
    } else {
     # shinyjs::logjs("brush is NULL")
     shinyjs::hide("highlight_genes")
     cat("Click-drag to select a rectangle in the plot area...")
    }
  })
  

  #
  # the following is TBC
  # 
  output$click_plot <- renderPlot({
    location <- nearPoints(tsne_data, input$plot_click, addDist = FALSE)
    if (nrow(location) == 0){
      cat("")
    }
    else{
    plot_name <- paste("/home/kw1016/MSc_Project2_Data_analysis/WebsiteApp/data_files/tsne_files/",rownames(location[1,]),".tsne.p30.txt",sep="")
    plot_data <- read.csv(plot_name, sep="\t", stringsAsFactors = F, header=FALSE)
    ggplot(plot_data, aes(V1, V2)) + geom_point(aes(color = pop_code)) +
      scale_colour_manual(values = c( "#ff7f00", #orange
                                      "#d8e544", #yellow
                                     
                                      "#984ea3", #purple  
                                      "#00e673", #green
                                      "#a65628", #brown
                                      "#999999", #grey
                                      "#00bfff", #blue
                                      "#f781bf", #pink
                                      "#e73032" #red
                          ),
                          labels = c("Angola", 
                                     "B.Faso(M)",
                                     "B.Faso(S)",
                                     "Camaroon",
                                     "Gabon",
                                     "Guinea",
                                     "G.Bissau",
                                     "Kenya",
                                     "Uganda"
                          ))+
      labs(x=format("t-SNE1", digits=2), y=format("t-SNE2", digits=2), colour=format("Population origin", digits=2)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
  })

  observeEvent(input$highlight_genes, {
      # shinyjs::logjs("observed highlight_genes")
      updateTextInput(session, "gene_id_search", value = genes$brushed)
  })

  observeEvent(input$reset_search, {
      # shinyjs::logjs("observed reset_search")
      updateTextInput(session, "gene_id_search", value = "")
      output$gene_id_search_result <- NULL
      shinyjs::hide("reset_search")
  })

})


