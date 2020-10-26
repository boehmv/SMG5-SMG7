#!/usr/bin/env Rscript

# Title: IRFinder
# Objective: Run DESeq2-like analysis downstream of IRFinder pipeline to identify log2FC of retained introns
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(tximport)
library(readr)
library(tximportData)
library(DESeq2)
library(Glimma)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(tibble)
library(gprofiler2)
library(data.table)
library(plotly)
library(htmlwidgets)
library(scales)
suppressPackageStartupMessages(library(circlize))
library(extrafont)
source("/opt/IRFinder-1.2.6/bin/DESeq2Constructor.R")

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition string is the second, gets converted to character
cond=arguments$args[2]
cond <- unlist(strsplit(cond,","))

# Create "Plots" folder if it does not exists
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)

# Create "DESeq2" folder if it does not exists
dir.create(file.path(mydir, "Plots", "IRFinder"), showWarnings = FALSE)


	# Get samples and check if all files are present
	results <- read.table(file.path(mydir, "IRFinder", "filePaths.txt"))
	paths = as.vector(results$V1)
	experiment <- read.table(file.path(mydir, "IRFinder", "experiment.txt"), header = TRUE)
	experiment$Condition=factor(experiment$Condition,levels=cond)
	rownames(experiment)=NULL
	
	# Generate a meta list containing three slots
		# First slot is a DESeq2 Object that can be directly pass to DESeq2 analysis.  
		# Second slot is a matrix for trimmed mean of intron depth
		# Third slot  is a matrix for correct splicing depth flanking introns
		# Fourth slot is a matrix for maximum splicing reads at either ends of introns
	metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
	# Extract DESeq2 Object with normalization factors ready
	dds = metaList$DESeq2Object
	# Build a formula of GLM. Read below for more details.
	design(dds) = ~Condition + Condition:IRFinder
	# Estimate parameters and fit to model
	dds = DESeq(dds)
	# Set working directory
	setwd(file.path(mydir, "IRFinder"))
	# Save R session to file for later	
	save.image(file='All_IRFinder_session.RData')

# Build for-loop to go through all conditions vs control
mycount=2
for (i in 2:length(cond)){
	# Tests if the number of IR reads are significantly different from normal spliced reads, in the control samples.
	res.control = results(dds, name = "Conditioncontrol.IRFinderIR")
	# Get the value of (intronic reads/normal spliced reads)
	control.IR_vs_Splice=2^res.control$log2FoldChange
	# Convert the above value to IR ratio
	IRratio.control = control.IR_vs_Splice/(1+control.IR_vs_Splice)
	#Same for condition
	res.condition = results(dds, name = paste0("Condition",paste(cond[mycount]),".IRFinderIR"))
	condition.IR_vs_Splice=2^res.condition$log2FoldChange
	IRratio.condition = condition.IR_vs_Splice/(1+condition.IR_vs_Splice)
	# Test the difference of (intronic reads/normal spliced reads) ratio between control and condition
	res.diff = results(dds, contrast=list("Conditioncontrol.IRFinderIR", paste0("Condition",paste(cond[mycount]),".IRFinderIR")))		
	# Convert results file to data.frame
	res.diff.DF <- as.data.frame(res.diff)
	# Rownames as column
	res.diff.DF <- rownames_to_column(res.diff.DF, var = "IDs")
	# Split rownames
	res.diff.DF <- res.diff.DF %>% separate(IDs, c("symbol", "gene_ID","status", "coordinates"), sep = "/")
	# Create condition-specific subfolder
	dir.create(file.path(mydir, "IRFinder",paste(cond[mycount])), showWarnings = FALSE)	
	# Generate output csv file
	write.csv(res.diff.DF, file=(file.path(mydir, "IRFinder", cond[mycount] ,paste0(paste(cond[mycount]),"_vs_control_IRFinder_results.csv"))))
	# Make Combined folder
	dir.create(file.path(mydir, "IRFinder", "Combined"), showWarnings = FALSE)
	# Create symbolic link if needed
	if(file.exists(file.path(mydir, "IRFinder", "Combined", paste0(paste(cond[mycount]),"_vs_control_IRFinder_results.csv")))){
	print("Symlink already established")
	} else {
	file.symlink(file.path(mydir, "IRFinder", cond[mycount] ,paste0(paste(cond[mycount]),"_vs_control_IRFinder_results.csv")), file.path(mydir, "IRFinder", "Combined", paste0(paste(cond[mycount]),"_vs_control_IRFinder_results.csv")))
	}
	
	# Generate Glimma MD plot
	#status <- as.numeric(res.diff.$padj < .05)
	#anno <- data.frame(GeneID=res.diff.DF$gene_ID, symbol=res.diff$symbol, coordinates=res.diff$coordinates)
	#glMDPlot(res.diff, status=status, counts=counts(dds,normalized=TRUE),
         #groups=dds$condition, transform=TRUE,
         #samples=colnames(dds), anno=anno,
         #path=(file.path(mydir, paste(cond[mycount]))), folder="glimma-MD", launch=FALSE)

	dir.create(file.path(mydir, "Plots", "IRFinder", paste(cond[mycount])), showWarnings = FALSE)


	# Filter results by padj
	res.diff.DF <- res.diff.DF[order(res.diff.DF$padj),]
	res.diff.DF_filt <- res.diff.DF %>% filter(padj < 0.05)

	# Generate volcano plot with plotly
	p <- plot_ly(data = res.diff.DF_filt, x = ~log2FoldChange, y = -log10(res.diff.DF_filt$padj), type = 'scatter', mode = 'markers', hoverinfo = 'text', color = -log10(res.diff.DF_filt$padj), colors=c("orange","red"), text = ~paste('</br> Gene: ', symbol, '</br> Gene ID: ', gene_ID, '</br> BaseMean: ', baseMean, '</br> padj: ', padj, '</br> log2FC: ', log2FoldChange, '</br> Coordinates: ', coordinates), marker = list(opacity = 0.6, size = log2(res.diff.DF_filt$baseMean))) %>%  layout(title = paste(cond[mycount]), yaxis = list(zeroline = FALSE, title = "-log10 padjust"), xaxis = list(zeroline = FALSE, title = "log2FoldChange"))

	# Set working directory
	setwd(file.path(mydir, "Plots", "IRFinder", paste0(paste(cond[mycount]))))
	# Publish the interactive volcano plot
	htmlwidgets::saveWidget(widget=p, (paste0(paste(cond[mycount]), "_vs_control_volcano.html")), title = paste(cond[mycount]))

	# Create unfiltered volcano plot

	res.diff.DF_all <- res.diff.DF %>% mutate(padj = replace(padj, padj == 0, 1e-320))

	myplot <- ggplot(data = res.diff.DF_all, aes(x=log2FoldChange, y=-log10(padj))) + 
	theme_classic() + 	
  	theme(legend.position="top", legend.background = element_rect(size=0.5, linetype="solid", colour ="black"), axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) +
	coord_cartesian(xlim = c(-25, 25), ylim = c(0, 330), clip = 'off') +
	geom_point(alpha = 0.5, aes(color=log2FoldChange)) + 
	scale_color_gradientn(limits = c(-25,25), colors = c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"), values = rescale(c(-25,-3, -1.5, 0, 1.5, 3, 25), to = c(0, 1)), name = "log2 FoldChange") +	
 	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	geom_hline(aes(yintercept=320), linetype="dashed", color="gray") +
  	geom_vline(xintercept = 1, linetype="dashed") +
  	geom_vline(xintercept = -1, linetype="dashed") +
 	labs(title = paste(cond[mycount]), subtitle = "IRFinder IR analysis", x = "log2 FoldChange", y = "-log10 p.adjust", caption = "IRFinder (1.2.6)")
	
	myplot <- myplot + annotate("text", x=20, y=340, label="Max p.adjust", size =3, color = "gray")

	myplot

	# Save plot
	ggsave(file.path(mydir, "Plots", "IRFinder", paste(cond[mycount]), paste0(paste(cond[mycount]), "_unfiltered_Volcano.pdf")), width = 4, height=4, device=cairo_pdf)

	# Save R session to file for later	
	save.image(file=paste0(paste(cond[mycount]),'IRFinder_session.RData'))

	# Counter increment
	mycount = mycount + 1
}

