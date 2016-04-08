
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)


shinyUI(fluidPage(
	
	titlePanel("Interactive report for prioritizing variant interpretation based on parameters below."),
	
	plotOutput("main_plot"),
	
	
	tabsetPanel(
		id = 'variant_table',
		tabPanel('hide table'),
		tabPanel('variants table',
			 DT::dataTableOutput('variants_table')
		),
		tabPanel('select table columns', 
			uiOutput("choose_columns")
		),
		tabPanel('Gene/Disease Filters',
		 	column(12,
			   	tags$div(title="Quality score in database of gene/disease association",
		 	 		sliderInput("gd_score", label = "Quality Score", min = 0,max=1,value=c(0.1,1),step=.01)
				),
				tags$div(title="Minimum number of PubMed Citations associating gene with disease",
					uiOutput("choose_pubmed")
				),
				DT::dataTableOutput('gene_disease_table')
		 	)
		)
	),

	fluidRow(
		column(4,
				h4("Show variants in samples:"),
				uiOutput("choose_include_samples"),
			   	checkboxInput("shared",label="Only include shared variants (intersection)",value=FALSE)
		),
		column(4,
			   	h4("Exclude variants in:"),
				uiOutput("choose_exclude_samples")
			   	#checkboxGroupInput("exclude_samples",label="",choices=mylist)
		),
		column(4,
			   h4("Filter options"),
			   tags$div(title="A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 5%",
			   		 checkboxInput("rare",label="Exclude common variants",value=FALSE)
			   ),
			   tags$div(title="Prediction of dbNSFP SVM based ensemble prediction score, 'T(olerated)' or'D(amaging'. The score cutoff between 'D' and 'T' is 0. The rankscore cutoff between 'D' and 'T' is 0.82268.",
			   		 checkboxInput("dbNSFP_deleterious",label="Only dbNSFP deleterious mutations",value=FALSE)
			   )
		)
	),
	fluidRow(
		column(4,
			   tags$div(title="Variant Quality score",
				uiOutput("slider_quality")
			   ),
			   tags$div(title="CADD phred-like score. 10 corresponds to top 10% predicted most damaging variants, 20 to 1%, and so on. This is phred-like rank score based on whole genome CADD raw scores. Please refer to Kircher et al. (2014) Nature Genetics 46(3):310-5  for details. The larger the score the more likely the SNP has damaging effect.", 
			   		 sliderInput("CADD_phred", label = "CADD phred score", min = 0, max = 40, value = 0,step=5)
			   ),
			actionButton("Export",label="Export variants"),
			textOutput("export_success")
		),
		column(4,
			tags$div(title="Minimum read depth of coverage.", 
				uiOutput("slider_depth")
			)
		),
		column(4,
			   tags$div(title="Fraction of reads supporting the alternate (non-reference) allele.",
			   		 sliderInput("alt_af",label="Alternate allele frequency",min=0,max=1,value=c(0,1),step=0.05)
			   )
		)
	)
		
	
	# Adapted solution from https://github.com/rstudio/shiny/issues/167
	
))
