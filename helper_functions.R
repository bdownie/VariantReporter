
cleanClinvar <- function(df) {
	df$CLNDBN <- gsub("\\|.*","",as.character(df$CLNDBN))
	df$CLNDBN <- gsub(",.*","",as.character(df$CLNDBN))
	df$CLNDBN <- gsub("\\\\x2c","",as.character(df$CLNDBN))
	df <- df[df$CLNDBN != "not_specified" & df$CLNDBN != "not_provided",]
	tmp <- strsplit(as.character(df$CLNSIG),"\\||,")
	to.output <- numeric(length=length(tmp))
	for (i in 1:length(tmp)) {
		v <- tmp[[i]]
		for (j in 1:length(v)) {
			#	if (regexpr("\\D",v[j],perl=TRUE)
			
			x <- as.numeric(v[j])
			if (!is.na(x) & x > 3 & x < 8) {
				to.output[i] <- x
			}
		}
	}
	df$CLNSIG <- to.output
	df <- df[df$CLNSIG > 0,]
	df
}

cleanAnnotation <- function(df) {
	if (!"splice_acceptor_variant" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"splice_acceptor_variant") }
	if (!"splice_donor_variant" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"splice_donor_variant") }
	if (!"frameshift_variant" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"frameshift_variant") }
	if (!"start_lost" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"start_lost") }
	if (!"stop_gained" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"stop_gained") }
	if (!"stop_lost" %in% levels(df$Annotation)) { levels(df$Annotation) <- c(levels(df$Annotation),"stop_lost") }
	
	df <- droplevels(df)
	
	df$Annotation <- sub("&.*","",as.character(df$Annotation),perl=TRUE)
	#df[df$Annotation == "frameshift_variant&splice_region_variant",]$Annotation <- "frameshift_variant"
	#df[df$Annotation == "frameshift_variant&start_lost",]$Annotation <- "start_lost"
	#df[df$Annotation == "frameshift_variant&stop_gained",]$Annotation <- "stop_gained"
	#df[df$Annotation == "splice_acceptor_variant&intron_variant",]$Annotation <- "splice_acceptor_variant"
	#df[df$Annotation == "splice_acceptor_variant&splice_donor_variant&intron_variant",]$Annotation <- "splice_acceptor_variant"
	#df[df$Annotation == "splice_acceptor_variant&splice_region_variant&intron_variant",]$Annotation <- "splice_acceptor_variant"
	#df[df$Annotation == "splice_donor_variant&intron_variant",]$Annotation <- "splice_donor_variant"
	#df[df$Annotation == "splice_acceptor_variant&splice_donor_variant&splice_region_varia",]$Annotation <- "splice_acceptor_variant"
	#df[df$Annotation == "frameshift_variant&stop_lost",]$Annotation <- "stop_lost"
	
	df <- droplevels(df)
}


plotAnnotation <- function(dfs,names="Sample") {
	library(ggplot2)
	if (length(names) != length(dfs)) {
		stop("Bad parameters")
	}
	tmp <- list()
	if (is.list(dfs)) {
		for (i in 1:length(dfs)) {
			df <- dfs[[i]]
			tmp[[i]] <- data.frame(Annotation=table(df$Annotation))
			tmp[[i]]$sample <- names[[i]]
		}
		to.plot <- do.call(rbind,tmp)
	}
	else { 
		to.plot <- data.frame(Annotation=table(dfs$Annotation_Impact))
		to.plot$sample <- names
	}
	colnames(to.plot) <- c("Annotation","Frequency","Sample")
	print(ggplot(to.plot,aes(x=Sample,fill=Annotation,y=Frequency)) + geom_bar(stat="identity") + ylab("Number of variants") + xlab("Sample name") + ggtitle("Distribution of high impact variants by class"))
}


plotDbSNP <- function(dfs,names="Sample") {
	library(ggplot2)
	if (length(names) != length(dfs)) {
		stop("Bad parameters")
	}
	tmp <- list()
	if (is.list(dfs)) {
		for (i in 1:length(dfs)) {
			df <- dfs[[i]]
			if (is.factor(df$COMMON)) {
				levels(df$COMMON) <- c(levels(df$COMMON),"0")
				df[!is.na(df$COMMON),] <- "1"
				df[is.na(df$COMMON),] <- "0"
				df$COMMON <- droplevels(df$COMMON)
				df$COMMON <- as.factor(as.numeric(as.character(df$COMMON)))
				levels(df$COMMON) <- c("No","Yes")
			}
			tmp[[i]] <- data.frame(CommonSNP=df$COMMON)
			tmp[[i]]$sample <- names[[i]]
		}
		to.plot <- do.call(rbind,tmp)
	}
	else { 
		to.plot <- data.frame(Annotation=table(dfs$Annotation_Impact))
		to.plot$sample <- names
	}
	colnames(to.plot) <- c("CommonSNP","Sample")
	print(ggplot(to.plot,aes(x=Sample,fill=CommonSNP)) + geom_bar(stat="bin") + ylab("Number of variants") + xlab("Sample name") + ggtitle("Distribution of moderate/high impact variants by population frequency"))
}

reorderColumns <- function(df) {
	first_columns <- c("IGV","chr","site","ID","ref","alt","score")
	second_columns <- c("Gene_Name","Gene_ID","Feature_ID","Transcript_BioType","Annotation","Annotation_Impact","cDNA.pos.cDNA.length",
						"CDS.pos.CDS.length","AA.pos.AA.length","DP","HGVS.c","HGVS.p","sample","AF")
	col_order <- c(first_columns,second_columns)
	tmp <- df[,!colnames(df) %in% col_order]
	tmp <- tmp[,order(apply(tmp,2,function(x) { sum(is.na(x)) } ),colnames(tmp))]
	o <- pmatch(col_order,colnames(df))
	o <- o[!is.na(o)]
	df <- data.frame(df[,o],tmp)
	df
}

processVariants <- function(df) { 
	
	df <- df[!is.na(df$site),]
	df$unique_id <- paste(df$chr,df$site,df$Feature_ID,sep="_")
	df$site_id <- paste(df$chr,df$site,df$alt,sep="_")
	
	df <- df[!duplicated(df$unique_id),]
	rownames(df) <- df$unique_id
	df <- df[,order(colnames(df))]
	#df$site_code <- paste(df$chr,df$site,sep="_")
	df <- df[is.na(df$ERRORS),]
	
	df <- df[!is.na(df$Annotation_Impact),]
	df$Annotation_Impact <- factor(df$Annotation_Impact,levels=c("MODIFIER","LOW","MODERATE","HIGH"))
	df$rank <- as.numeric(df$Annotation_Impact)
	df <- df[order(df$rank,df$site_id,decreasing=TRUE),]
	#df <- df[!duplicated(df$site_id),]
	

	df <- df[,colSums(is.na(df)) < nrow(df)]
	#df <- df[,apply(df,2,function(x) { if (x == "sample") { TRUE } else { length(table(x)) > 1 } })]
	#for (i in ncol(df):1) {
	#	if (sum(is.na(df[,i])) == nrow(df)) {
	#		df <- df[-i]
	#	}
	#}
	
	df$transcript_length <- as.numeric(sub(".*\\/","",df$CDS.pos.CDS.length))
	df[is.na(df$transcript_length),]$transcript_length <- 0
	
	df <- reorderColumns(df)
	
	df
}



plotVariantImpact <- function(dfs,names="Sample") {
	library(ggplot2)
	if (length(names) != length(dfs)) {
		stop("Bad parameters")
	}
	tmp <- list()
	if (is.list(dfs)) {
		for (i in 1:length(dfs)) {
			df <- dfs[[i]]
			tmp[[i]] <- data.frame(Annotation=table(dfs[[i]]$Annotation_Impact))
			tmp[[i]]$sample <- names[[i]]
		}
		to.plot <- do.call(rbind,tmp)
	}
	else { 
		to.plot <- data.frame(Annotation=table(dfs$Annotation_Impact))
		to.plot$sample <- names
	}
	colnames(to.plot) <- c("Impact","Frequency","Sample")
	print(ggplot(to.plot,aes(y=Frequency,fill=Impact,x=Sample)) + geom_bar(stat="identity") + ylab("Number of variants") + xlab("Sample name") + ggtitle("Distribution of variants by impact"))
}

# Taken from http://stackoverflow.com/questions/13123638/there-is-pmin-and-pmax-each-taking-na-rm-why-no-psum
psum <- function(...,na.rm=FALSE) { 
	rowSums(do.call(cbind,list(...)),na.rm=na.rm) } 

# Taken from http://stackoverflow.com/questions/28117556/clickable-links-in-shiny-datatable
createLink <- function(val) {
	sprintf('<a href="http://localhost:60151/goto?locus=%s" target="trash" class="btn btn-primary">IGV</a>',val)
}

subsetDataTable <- function(df,columnsToShow=1:5,genesToShow=NULL) {
	if (length(columnsToShow) == 0) { 
		df <- df[,1:5]
	} else { 
		df <- df[,colnames(df) %in% columnsToShow]
	}
	if (length(genesToShow) > 0) {
		df <- df[df$Gene_name %in% genesToShow,]
	}
	df
}

cleanAndFilterVariants <- function(df) {
	#df <- data.unfiltered.reactive()
	#df <- df[sample.int(nrow(df),size=1000),]
		
		
	
	#print("0")
	#print(head(df))
	df <- df[!is.na(df$rank),]
	df <- df[!is.na(df$transcript_length),]

	df <- df[order(df$rank,decreasing=TRUE),]
	df.list <- list()
	for (i in 1:length(levels(df$sample))) {
		tmp <- df[df$sample == levels(df$sample)[i],]
		df.list[[i]] <- tmp[!duplicated(tmp$site_id),]
	}
	#print("2")
	df <- do.call(rbind,df.list)
	df <- df[df$gt != "./.",]
	df$REF_DP <- sub(":.*","",sub(".*?:","",df$gt,perl=TRUE),perl=TRUE)
	#df[is.na(df$REF_DP),]$REF_DP 
	df$ALT_DP <- as.numeric(unlist(lapply(strsplit(df$REF_DP,split=","),function(x) { sum(as.numeric(x[-1])) })))
	df$REF_DP <- as.numeric(unlist(lapply(strsplit(df$REF_DP,split=","),function(x) { as.numeric(x[1]) })))
	df$ALT_RATIO <- df$ALT_DP/psum(df$ALT_DP,df$REF_DP)
	df[is.na(df$ALT_RATIO),]$ALT_RATIO <- as.numeric(sub(",.*","",df[is.na(df$ALT_RATIO),]$AF),perl=TRUE)
	df$IGV <- sapply(paste(df$chr,df$site,sep=":"), function(x) { createLink(x) })
	df <- df[,colSums(is.na(df)) < nrow(df)]
	df <- reorderColumns(df)
	df
}

readSamples <- function(sampleNames,fileLocations) {
	tables <- list()
	for (i in 1:length(fileLocations)) {
		tables[[i]] <- processVariants(data.frame(read.table(fileLocations[i],header=TRUE),sample=sampleNames[i]))
	}
	for (i in 1:length(tables)) {
		for (j in 1:length(tables)) {
			if (i != j) {
				x <- tables[[i]]
				y <- tables[[j]]
				tables[[i]] <- x[colnames(x) %in% colnames(y)]
			}
		}
	}
	#if (exists(x) && exists(y)) {
	#	rm(x)
	#	rm(y)
	#}
	df <- do.call(rbind,tables)
	df
	
			
#	save(data.unfiltered,file=data_loc)
}

