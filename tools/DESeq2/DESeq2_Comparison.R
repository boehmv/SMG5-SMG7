#!/usr/bin/env Rscript

# Title: DESeq2_Comparison
# Objective: Comparison DESeq2 pipeline, creating volcano, heatmap and barcode plots
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(dplyr)
library(limma)
library(scales)
library(ggplot2)
library(ggrepel)
library(ggplotify)
suppressPackageStartupMessages(library(ComplexHeatmap))
library("RColorBrewer")
suppressPackageStartupMessages(library(circlize))
library(extrafont)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition table as second argument, condition followed by absolute path to condition-specific DESeq2 csv output file, seperated by tab
myComp=arguments$args[2]

# Create "Plots" folder if it does not exists
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)

# Create "DESeq2" folder if it does not exists
dir.create(file.path(mydir, "Plots", "DESeq2"), showWarnings = FALSE)

# Read in file 
df <- read.table(file = myComp, sep = '\t', header = TRUE)

# Generate empty data frame as reference
combined <- data.frame(matrix(ncol = 2, nrow = 0))
x <- c("geneID", "symbol")
colnames(combined) <- x

# Generate empty data frame for lfc
lfc_df <- data.frame(matrix(ncol = 1, nrow = 0))
y <- c("geneID")
colnames(lfc_df) <- y

# Generate empty data frame for concatenated file
DGE_cat <- data.frame(matrix(ncol = 6, nrow = 0))
y <- c("geneID","baseMean","log2FoldChange","padj","symbol", "condition_2")
colnames(DGE_cat) <- y

# Iterate over conditions and generate data frames, generate combined df by merge
for (row in 1:nrow(df)) {
	cond=df[row, "cond"]
	path=df[row, "path"]

	d <- read.table(file = file.path(path), sep = ',', header = TRUE)
	d_short <- subset(d, abs(log2FoldChange) > 1 & padj < 0.05, select = c("geneID","log2FoldChange","padj","symbol"))
	d_short <- subset(d_short, !grepl("^SIRV", geneID))
	d_short <- subset(d_short, !grepl("^ERCC", geneID))
	d_long <- subset(d, !grepl("^SIRV", geneID) & !grepl("^ERCC", geneID), select = c("geneID","log2FoldChange","padj","symbol"))

	d_cat <- subset(d, !grepl("^SIRV", geneID) & !grepl("^ERCC", geneID), select = c("geneID", "baseMean","log2FoldChange","padj","symbol"))
	d_cat$condition_2 <- paste0(cond)
	DGE_cat <- rbind(DGE_cat, d_cat)

	# Create condition-specific folder if it does not exists
	dir.create(file.path(mydir, "Plots", "DESeq2", paste0(cond)), showWarnings = FALSE)

	# Create volcano filtered plot
	ggplot(data = d_short, aes(x=log2FoldChange, y=-log10(padj))) + 
	theme_classic() + 	
	  theme(strip.background = element_blank(), 
	        plot.background = element_blank(),  
	        panel.background = element_blank(),
	        legend.position="top",
	        legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
	        axis.text=element_text(size=6, color="black"),
	        axis.title=element_text(size=6),
	        legend.title = element_text(size = 6),
	        legend.text = element_text(size = 6),
	        plot.title = element_text(size = 6),
	        plot.subtitle = element_text(size = 6),
	        plot.caption = element_text(size = 5),
	        text=element_text(family="Arial"), 
	        axis.line = element_line(colour = 'black', size = 0.1), 
	        axis.ticks = element_line(colour = "black", size = 0.1)) +
	geom_point(alpha = 0.5, aes(color=log2FoldChange)) + 
	scale_color_gradientn(colors = c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"), values = rescale(c(min(d_short$log2FoldChange),-3, -1.5, 0, 1.5, 3, max(d_short$log2FoldChange)), to = c(0, 1)), name = "log2 FoldChange") +	
 	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  	geom_vline(xintercept = 1, linetype="dashed") +
  	geom_vline(xintercept = -1, linetype="dashed") +
 	labs(title = paste0(cond), subtitle = "DESeq2 DEG analysis", x = "log2 FoldChange", y = "-log10 p.adjust", caption = "DESeq2")
	
	# Save plot
	ggsave(file.path(mydir, "Plots", "DESeq2", paste0(cond), paste0(paste(cond), "_Volcano.pdf")), width = 8, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

	# Create unfiltered volcano plot

	d_long <- d_long %>% mutate(padj = replace(padj, padj == 0, 1e-320))

	myplot <- ggplot(data = d_long, aes(x=log2FoldChange, y=-log10(padj))) + 
	theme_classic() + 	
  	theme(strip.background = element_blank(), 
  	      plot.background = element_blank(), 
  	      panel.background = element_blank(),
  	      legend.position="top",
  	      legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
  	      axis.text=element_text(size=6, color="black"),
  	      axis.title=element_text(size=6),
  	      legend.title = element_text(size = 6),
  	      legend.text = element_text(size = 6),
  	      plot.title = element_text(size = 6),
  	      plot.subtitle = element_text(size = 6),
  	      plot.caption = element_text(size = 5),
  	      text=element_text(family="Arial"),
  	      axis.line = element_line(colour = 'black', size = 0.1), 
  	      axis.ticks = element_line(colour = "black", size = 0.1)) +
	coord_cartesian(xlim = c(-25, 25), ylim = c(0, 330), clip = 'off') +
	geom_point(alpha = 0.5, aes(color=log2FoldChange)) + 
	scale_color_gradientn(limits = c(-25,25), colors = c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"), values = rescale(c(-25,-3, -1.5, 0, 1.5, 3, 25), to = c(0, 1)), name = "log2 FoldChange") +	
 	geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	geom_hline(aes(yintercept=320), linetype="dashed", color="gray") +
  	geom_vline(xintercept = 1, linetype="dashed") +
  	geom_vline(xintercept = -1, linetype="dashed") +
 	labs(title = paste0(cond), subtitle = "DESeq2 DEG analysis", x = "log2 FoldChange", y = "-log10 p.adjust", caption = "DESeq2")
	
	myplot <- myplot + annotate("text", x=20, y=340, label="Max p.adjust", size =3, color = "gray")

	# Save plot
	ggsave(file.path(mydir, "Plots", "DESeq2", paste0(cond), paste0(paste(cond), "_unfiltered_Volcano.pdf")), width = 8, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

	d_sig <- subset(d, padj < 0.05, select = c("geneID","log2FoldChange"))
	names(d_sig)[names(d_sig) == "log2FoldChange"] <- paste0("log2FC_", paste(cond))	
	d_sig <- subset(d_sig, !grepl("^SIRV", geneID))
	d_sig <- subset(d_sig, !grepl("^ERCC", geneID))
	lfc_df <- merge(x = lfc_df, y = d_sig, by = c("geneID"), all = TRUE)

	names(d_short)[names(d_short) == "log2FoldChange"] <- paste0("log2FC_", paste(cond))
	names(d_short)[names(d_short) == "padj"] <- paste0("padj_", paste(cond))	
	combined <- merge(x = combined, y = d_short, by = c("geneID", "symbol"), all = TRUE)
	assign(paste(cond), d)
}

write.csv(DGE_cat, file=file.path(mydir, "DESeq2", "DGE_cat.csv"))

# Select for heatmap generation only log2foldchange columns
combined_heat <- combined %>% dplyr:: select(starts_with("log2FC"))
colnames(combined_heat) <- df[, "cond"]


# Remove NA (replace with 0)
combined_heat[is.na(combined_heat)] <- 0

# Order data frame
HeatcolSums <- colSums(combined_heat != 0)

HeatcolSums <- HeatcolSums[order(-HeatcolSums)]

combined_heat <- combined_heat[,names(HeatcolSums)]

i <- length(HeatcolSums)

while (i>0) {
combined_heat <- combined_heat[order(-combined_heat[,i]),]
i=i-1
}

# Get Annotations of Up-/Downregulated genes
UpReg <- colSums(combined_heat > 1)
UpAnno = HeatmapAnnotation("Genes Up" = anno_barplot(UpReg, gp = gpar(fill="#253494"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

DownReg <- colSums(combined_heat < -1)
DownAnno = HeatmapAnnotation("Genes Down" = anno_barplot(DownReg, gp = gpar(fill="#980043"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

######
# Heatmap
######
	# Set Color for Heatmap
	col = colorRamp2(c(min(combined_heat),-3, -1.5, 0, 1.5, 3, max(combined_heat)), c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"))

	# Generate Heatmap
	heat <- Heatmap(as.matrix(combined_heat),
	name = "log2FC",
	heatmap_legend_param = list(title_position = "topleft", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
	col =  col,
	na_col = "white",
	cluster_rows = FALSE,
	cluster_columns = FALSE,	# Set to TRUE if clustering is desired
	column_dend_side = "top",
	border = TRUE,
	show_row_names = FALSE,
	show_row_dend = FALSE,
	column_dend_reorder = FALSE,
	row_dend_reorder = FALSE,
	column_names_gp = gpar(fontsize=6),
	top_annotation = UpAnno,
	bottom_annotation = DownAnno,
	width = unit((length(HeatcolSums)/2), "cm"),
	height = unit(5, "cm")
	)
	
	# Print to pdf
	pdf(file = file.path(mydir, "Plots", "DESeq2", "DEseq2_Heatmap.pdf"))

	draw(heat)

	dev.off()

######
# Barcode plots
######

# Get conditions
n <- length(df[,"cond"])
myCond <- as.vector(df[,"cond"])

# Loop over all conditions as reference
while (n>0) {
	# Set reference name
	reference <- as.character(myCond[n])
	referenceName <- paste0("log2FC_", paste(reference))

	# Create condition-specific folder if it does not exists
	dir.create(file.path(mydir, "Plots", "DESeq2", paste0(reference)), showWarnings = FALSE)
	
	# Order by lfc_df by reference
	lfc_df <- lfc_df[order(-lfc_df[,(n+1)]),]
	lfc_df_n <- subset(lfc_df, select = c("geneID",referenceName))
	lfc <- lfc_df_n[,referenceName]
	names(lfc) <- lfc_df_n$geneID
	
	# Get remaining conditions
	ConCond <- myCond[myCond != reference]

	# Loop over remaining conditions
	p <- length(ConCond)
	while (p>0) {
		Up_Index <- vector(,length(lfc))
		Down_Index <- vector(,length(lfc))
		condName <- paste0("log2FC_", paste(ConCond[p]))
		Up_Cond <- combined$geneID[which(combined[[condName]] > 1)]
		Down_Cond <- combined$geneID[which(combined[[condName]] < -1)]

		for (i in 1:length(lfc)) {
			if (names(lfc[i]) %in% Up_Cond) {
				Up_Index[i] <- TRUE
			} else { 
				Up_Index[i] <- FALSE	
			}
			if (names(lfc[i]) %in% Down_Cond) {
				Down_Index[i] <- TRUE
			} else {
				Down_Index[i] <- FALSE
			}
		}
		# Write barcode plot to pdf
		pdf(file = file.path(mydir, "Plots", "DESeq2", paste0(reference), paste0(paste(reference), "_vs_", paste(ConCond[p]), "_Barcode.pdf")), useDingbats=FALSE)

		myplot <- as.ggplot(~barcodeplot(lfc, 
                index = Up_Index, 
                index2 = Down_Index,  
                xlab = "log2FoldChange",
                labels = c("",""),
                worm = TRUE,
                par(ps = 6, cex = 1, cex.lab = 0.8, cex.axis = 1)
                ), scale = 1, hjust = 0, vjust = 0)
		
		myplot <- myplot + 
		labs(title = paste0(paste(reference), " vs ", paste(ConCond[p])), subtitle = "DESeq2 Barcode plot", caption = "DESeq2 (1.24.0)")

		ggsave(file.path(mydir, "Plots", "DESeq2", paste0(reference), paste0(paste(reference), "_vs_", paste(ConCond[p]), "_Barcode.pdf")), plot = myplot, width = 6, height=4, device=cairo_pdf)

		p=p-1	
	}
	n=n-1
}

writeLines(capture.output(sessionInfo()), paste0(mydir, "/DESeq2/DESeq2_comparison_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
