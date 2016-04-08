
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
require(ggplot2)
library(scales)
source("helper_functions.R")
library(DT)
load("gene_disease.Rdata")
#library(data.table)

shinyServer(function(input, output) {
	attr(input, "readonly") <- FALSE
	
	#input$depth <- reactive({ 
	#	if (isdepth <- 0
	#	})
	
	val <- read.csv("parameters.course.csv",header=TRUE,stringsAsFactors = FALSE)
	sampleNames <- as.character(val[,1])
	files <- as.character(val[,2])
	last_export <- 0
	
	
	data_loc <- paste0(getwd(),"/.table.Rdata")
	#data.unfiltered <- reactive({
	data.filtered <- reactive({
		if (file.exists(data_loc)) {
			#print("Using cached data (delete .table.Rdata if you want to reprocess the data)")
			load(data_loc)
		} else { 
			#print("Importing data")
			#print(data_loc)
			data.unfiltered <- readSamples(sampleNames=sampleNames,fileLocations=files)
	
			data.filtered <- cleanAndFilterVariants(data.unfiltered)
			save(data.unfiltered,data.filtered,file=data_loc)
		}
		data.filtered
	})

	#data.filtered <- head(data.filtered,n=1000)
	sampleList <- 1:length(sampleNames)
	names(sampleList) <- sampleNames
	#print(sampleList)
	
	output$choose_include_samples <- renderUI ({ 
		#print(head(data.filtered))
		   checkboxGroupInput("include_samples",label="",choices=sampleList,selected=1:length(sampleList))
	})
	output$choose_exclude_samples <- renderUI ({ 
		   checkboxGroupInput("exclude_samples",label="",choices=sampleList)
	})
	
	output$slider_quality <- renderUI ({ 
		sliderInput("quality", label = "Variant Quality Score", min = 0,max=floor(quantile(data.filtered()$score,probs=0.25)/10)*10,value=0,step=10)
	})
	output$slider_depth <- renderUI ({ 
		sliderInput("depth", label = "Minimum read depth", min = 0,max=floor(quantile(data.filtered()$DP,probs=0.9)/10)*10,value=0,step=10)
#		input$depth <- 0
	})
	
	to.process <- reactive({
		#time0 <- Sys.time()
		
		#print(time0)
		if (is.null(input$include_samples)) { 
			to.plot <- data.filtered()
		} else {
			to.plot <- data.filtered()[data.filtered()$sample %in% sampleNames[as.numeric(input$include_samples)],]
		}
	#	print(head(data.filtered()))
	#	print(head(input$include_samples))
		if (input$shared) {
			to.plot$sample <- droplevels(to.plot$sample)
			tab <- table(to.plot$site_id)
			to.plot <- to.plot[to.plot$site_id %in% names(tab[tab == length(levels(to.plot$sample))]),]
			rm(tab)
		}
		if (length(input$exclude_samples) > 0) {
			tmp3 <- data.filtered()[data.filtered()$sample %in% sampleNames[as.numeric(input$exclude_samples)],]
		#	print(head(tmp3$site_id))
			to.plot <- to.plot[!to.plot$site_id %in% tmp3$site_id,]
			rm(tmp3)
		}
		if (input$rare) { 
			to.plot <- to.plot[is.na(to.plot$G5),]
		} 
		if (input$dbNSFP_deleterious) { 
			to.plot <- to.plot[!is.na(to.plot$dbNSFP_MetaSVM_pred) & to.plot$dbNSFP_MetaSVM_pred == 'D',]
		} 
		if (input$CADD_phred > 0) { 
			to.plot <- to.plot[!is.na(to.plot$dbNSFP_CADD_phred) & as.numeric(to.plot$dbNSFP_CADD_phred) >= input$CADD_phred,]
		} 
		if (input$alt_af[1] > 0 | input$alt_af[2] < 1) { 
			to.plot <- to.plot[to.plot$ALT_RATIO >= input$alt_af[1] & to.plot$ALT_RATIO <= input$alt_af[2], ]
		}
		#if (is.null(input$depth)) { 
		#	print("here")
		#}
		if ((!is.null(input$depth)) && input$depth > 0) { 
		#if (input$depth > 0) { 
			to.plot <- to.plot[to.plot$DP >= input$depth,]
		}
		if ((!is.null(input$quality)) && input$quality > 0) { 
		#if (input$quality > 0) { 
			to.plot <- to.plot[to.plot$score >= input$quality,]
		}
		
		
		if (!is.null(input$Export) && input$Export) {
			#output$export_output <- renderPrint({ print("") })
			#input$ExportMemory<- input$Export
			#if (input$ExportMemory) {
			require(xlsx)
			to.output <- data.filtered()[rownames(data.filtered()) %in% rownames(to.plot),]
			for (i in ncol(to.output):1) {
				if (sum(is.na(to.output[,i])) == nrow(to.output)) {
					to.output <- to.output[-i]
				}
			}
			to.output <- to.output[order(to.output$Gene_ID,to.output$site, as.numeric(to.output$Annotation_Impact), to.output$transcript_length,to.output$Feature_ID,decreasing=TRUE),]
			
			#WriteXLS("to.output","variants.xlsx")
			#write.table()
			write.xlsx(x=to.output,file="variants.xlsx")
			output$export_output <- renderPrint({ print("Exported successfully") })
			input$Export<- 0
			Sys.sleep(3)
			output$export_output <- renderPrint({ print("") })
		}
		to.plot <- to.plot[,colSums(is.na(to.plot)) < nrow(to.plot)]
		#to.plot <- droplevels(to.plot)
		#to.plot2 <- to.plot[,apply(to.plot,2,function(x) { length(table(x)) }) > 1]
	#	print(paste("Time here:",Sys.time() - time0))
	#	time0 <- Sys.time()
		#if (!("Annotation_Impact" %in% colnames(to.plot2))) { to.plot2$Annotation_Impact <- to.plot$Annotation_Impact }
		#if (!("sample" %in% colnames(to.plot2))) { to.plot2$sample <- to.plot$sample }
	#	to.plot2 <- to.process()
	
		#print("Done")
		#to.plot <- to.plot2
		
		#print(head(to.plot))
		#rm(to.plot2)
		to.plot
	})
	
	output$main_plot <- renderPlot ({ 
		df <- to.process()
		#print(head(df))
		if (nrow(df) > 0) { 
			p <- ggplot(df,aes(x=sample,fill=Annotation_Impact))  + geom_bar(stat="bin")  + ylab("Number of variants")  + scale_y_continuous(labels=comma)  + scale_fill_brewer(palette="Set2")
			print(p)
		}
			  
	})
		
	#output$variantColumns <- reactive ({
	#	tmp <- droplevels(to.process())
	#	tmp <- tmp[,apply(tmp,2,function(x) { length(table(x)) > 1 })]
	#	tmp
	#})
	
	output$choose_columns <- renderUI ({ 
		#print("here")
		tmp <- droplevels(to.process())
		tmp <- tmp[,apply(tmp,2,function(x) { length(table(x)) > 1 })]
		#print("here")
		checkboxGroupInput("table_columns",label="Columns to Display",choices=colnames(tmp),selected=colnames(tmp)[1:5])
	})
	
	output$choose_pubmed <- renderUI ({ 
		sliderInput("gd_pubmed", label = "Number of Pubmed Citations", min=0,max=quantile(geneDisease$NumberOfPubmeds,probs=0.99),value=5,step=1)
	})
	
	output$variants_table <- DT::renderDataTable ({
		df <- to.process()
		drawColumns <- input$table_columns
		if (length(drawColumns) == 0) { 
			drawColumns <- c("Gene_Name",colnames(df)[1:5])
		}  else if (!"Gene_Name" %in% drawColumns){
			drawColumns <- c("Gene_Name",drawColumns)
		}  
		
		df <- df[,colnames(df) %in% drawColumns]
		if (!is.null(input$gene_disease_table_rows_selected)) {
			geneSymbols <-  as.character(geneDisease[input$gene_disease_table_rows_selected,]$geneSymbol)
		} else { geneSymbols <- c() }
		
		if (length(geneSymbols) > 0) {
			#print(geneSymbols)
			df <- df[df$Gene_Name %in% geneSymbols,]
		}
		rownames(df) <- 1:nrow(df)
		df
	}, options = list(lengthMenu = c(10, 30, 50,100), pageLength = 10,autowidth=TRUE),escape=FALSE)
	
	output$gene_disease_table <- DT::renderDataTable ({
		#	pubmedMax <- input$gd_pubmed[2]
		#	if (pubmedMax == quantile(geneDisease$NumberOfPubmeds,probs=0.99)) { 
		#		pubmedMax = max(geneDisease$NumberOfPubmeds)
		#	}
		#eval(
		geneDisease[geneDisease$NumberOfPubmeds >= input$gd_pubmed[1] &  
						# geneDisease$NumberOfPubmeds <= pubmedMax & 
						geneDisease$score <= input$gd_score[2] & 
						geneDisease$score >= input$gd_score[1] & 
						geneDisease$geneSymbol %in% to.process()$Gene_Name ,
					]
		#)
		
	}, options = list(lengthMenu = c(10, 30, 50,100), pageLength = 10,autowidth=TRUE),escape=FALSE)
	
})
