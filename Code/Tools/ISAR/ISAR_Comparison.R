#!/usr/bin/env Rscript

# Title: ISAR_comparison
# Objective: Comparison ISAR pipeline giving extra plots (volcano, heatmap, barcodes) and analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(extrafont)
library(cowplot)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition table as second argument, condition followed by absolute path, seperated by tab
myComp=arguments$args[2]

# Create "Plots" folder if it does not exists
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)

# Create "ISAR" folder if it does not exists
dir.create(file.path(mydir, "Plots", "ISAR"), showWarnings = FALSE)

# Read in file 
df <- read.table(file = myComp, sep = '\t', header = TRUE)

# Generate empty data frame as reference
combined <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("isoform_id", "gene_name", "PTC")
colnames(combined) <- x

# Generate empty data frame for dIF
dIF_df <- data.frame(matrix(ncol = 1, nrow = 0))
y <- c("isoform_id")
colnames(dIF_df) <- y

# Initiate empty vector for All conditions
AllCond <- c()

# Get IsoformSwitchAnalyzeR version
ISAR_version <- packageVersion("IsoformSwitchAnalyzeR")

# Iterate over conditions and generate data frames, generate combined df by merge
for (row in 1:nrow(df)) {
	cond=as.character(df[row, "cond"])
	cond <- as.list(strsplit(cond, " ")[[1]])
	path=df[row, "path"]
	
	# Subset the data (filter for dIF and q_value, select desired columns and remove SIRV/ERCC hits)
	d <- read.table(file = file.path(path), sep = ',', header = TRUE)
	d_sig <- subset(d, isoform_switch_q_value < 0.05, select = c("isoform_id", "gene_name", "condition_2", "dIF","isoform_switch_q_value","PTC"))
	d_sig <- subset(d_sig, !grepl("^SIRV", isoform_id))
	d_sig <- subset(d_sig, !grepl("^ERCC", isoform_id))
	d_all <- subset(d, !grepl("^SIRV", isoform_id) & !grepl("^ERCC", isoform_id), select = c("isoform_id", "gene_name", "condition_2", "dIF","isoform_switch_q_value","PTC"))

	
	
	# Generate seperate data frames for each condition
	for (i in 1:length(cond)) {
		
		# Generate vector of all conditions
		AllCond <- c(AllCond, cond[i])

		# Subset the ISAR output for condition-specific rows
		dIF <- d_sig[which(d_sig$condition_2 == cond[i]),]

		dIF_all <- d_all[which(d_all$condition_2 == cond[i]),]		

		dIF_all <- dIF_all %>% mutate(isoform_switch_q_value = replace(isoform_switch_q_value, isoform_switch_q_value == 0, 1e-320))

		dIF_sig <- subset(dIF_all, isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)

		# Create condition-specific folder if it does not exists
		dir.create(file.path(mydir, "Plots", "ISAR", paste0(cond[i])), showWarnings = FALSE)

		# Create unfiltered volcano plot
		myplot <- ggplot(data = dIF_all, aes(x=dIF, y=-log10(isoform_switch_q_value))) + 
	  theme_classic() + 	
	  theme(legend.position="top", strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), legend.background = element_rect(size=0.5, linetype="solid", colour ="black"), axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 10), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) +
	  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 330), clip = 'off') +
	  geom_point(data=filter(dIF_all, is.na(PTC)), aes(color="lightgray"), alpha = 0.5) + 
	  geom_point(data=filter(dIF_all, PTC == FALSE), aes(color="royalblue"), alpha = 0.5) +
	  geom_point(data=filter(dIF_all, PTC == TRUE), aes(color="red"), alpha = 0.5) +
	  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
	  geom_hline(aes(yintercept=320), linetype="dashed", color="gray") +
	  geom_vline(xintercept = 0.1, linetype="dashed") +
	  geom_vline(xintercept = -0.1, linetype="dashed") +		
	  labs(title = paste0("Control_vs_",paste0(cond[i])), subtitle = "ISAR DTU analysis", x = "delta isoform fraction (dIF)", y = "-log10 q-value") +
	  scale_color_identity(name = "PTC",
                       breaks = c("lightgray", "royalblue", "red"),
                       labels = c("NA", "FALSE", "TRUE"),
                       guide = "legend")

	myplot <- myplot + annotate("text", x=-0.85, y=340, label="Max q-value", size =3, color = "gray")
	
		# Create filtered density plot
		myplot2 <- ggplot(data = dIF_sig, aes(dIF)) + 
  theme_classic() + 	
  theme(legend.position="none", strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 12), legend.text = element_text(size = 10), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) +
  xlim(-1, 1) +
  geom_density(data=filter(dIF_sig, is.na(PTC)), aes(color="lightgray"), fill="lightgray", alpha = 0.2) +
  geom_density(data=filter(dIF_sig, PTC == FALSE), aes(color="royalblue"), fill="royalblue", alpha = 0.2) +
  geom_density(data=filter(dIF_sig, PTC == TRUE), aes(color="red"), fill="red", alpha = 0.2) +
  labs(x = "delta isoform fraction (dIF)", y = "filtered density", caption = paste0("Filter: |dIF| > 0.1 & q-value < 0.05 \nISAR (", ISAR_version, ")")) +
  scale_color_identity(name = "PTC",
                       breaks = c("lightgray", "royalblue", "red"),
                       labels = c("NA", "FALSE", "TRUE"),
                       guide = "legend")

		finalplot <- plot_grid(myplot, myplot2, ncol = 1, align = "v", rel_heights = c(2, 1))
		finalplot		


		# Save plot
		ggsave(file.path(mydir, "Plots", "ISAR", paste0(cond[i]), paste0(paste(cond[i]), "_Volcano_Density.pdf")), width = 10, height=15, units = "cm", device=cairo_pdf, bg = "transparent")

		# Remove unneccessary columns
		dIF_filt <- subset(dIF, select = c("isoform_id", "gene_name", "dIF","isoform_switch_q_value","PTC")) 

		# Seperate subset to get just dIF, merge for barcode-plot data frame (dIF_df)
		dIF_Barcode <- subset(dIF, select = c("isoform_id","dIF"))
		names(dIF_Barcode)[names(dIF_Barcode) == "dIF"] <- paste0("dIF_", paste(cond[i]))
		dIF_df <- merge(x = dIF_df, y = dIF_Barcode, by = c("isoform_id"), all = TRUE)

		# Subset / filter dIF_filt by absolute dIF > 0.1
		dIF_filt <- subset(dIF_filt, abs(dIF) > 0.1)

		# Rename columns to give them condition-specific names
		names(dIF_filt)[names(dIF_filt) == "dIF"] <- paste0("dIF_", paste(cond[i]))
		names(dIF_filt)[names(dIF_filt) == "isoform_switch_q_value"] <- paste0("qvalue_", paste(cond[i]))	
	
		# Merge q-value/dIF filtered data 
		combined <- merge(x = combined, y = dIF_filt, by = c("isoform_id", "gene_name", "PTC"), all = TRUE)

		
		assign(paste0("ISAR_", paste(cond[i])), dIF)
	}	
}

write.csv(combined, file = file.path(mydir, "Plots", "ISAR", "ISAR_combined.csv"))

# Select for heatmap generation only dIF columns
combined_heat <- combined %>% dplyr:: select(starts_with("dIF"))
colnames(combined_heat) <- AllCond


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

# Get Annotations of |dIF| > 0.1
UpReg <- colSums(combined_heat > 0.1)
UpAnno = HeatmapAnnotation("Isoforms Up" = anno_barplot(UpReg, gp = gpar(fill="#253494"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

DownReg <- colSums(combined_heat < -0.1)
DownAnno = HeatmapAnnotation("Isoforms Down" = anno_barplot(DownReg, gp = gpar(fill="#980043"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

# Make Heatmap
	# Set Color for Heatmap
	col = colorRamp2(c(min(combined_heat),-0.3, -0.1, 0, 0.1, 0.3, max(combined_heat)), c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"))

	# Generate Heatmap
	heat <- Heatmap(as.matrix(combined_heat),
	name = "dIF",
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
	pdf(file = file.path(mydir, "Plots", "ISAR", "ISAR_Heatmap.pdf"))

	draw(heat)

	dev.off()

# Make Volcano plots

# Barcode plots

# Get conditions
n <- length(AllCond)
myCond <- as.vector(AllCond)

# Loop over all conditions as reference
while (n>0) {
	# Set reference name
	reference <- as.character(myCond[n])
	referenceName <- paste0("dIF_", paste(reference))

	# Create condition-specific folder if it does not exists
	dir.create(file.path(mydir, "Plots", "ISAR", paste0(reference)), showWarnings = FALSE)
	
	# Order lfc_df by reference
	dIF_df <- dIF_df[order(-dIF_df[,(n+1)]),]
	dIF_df_n <- subset(dIF_df, select = c("isoform_id",referenceName))
	dIF <- dIF_df_n[,referenceName]
	names(dIF) <- dIF_df_n$isoform_id
	
	# Get remaining conditions
	ConCond <- myCond[myCond != reference]

	# Loop over remaining conditions
	p <- length(ConCond)
	while (p>0) {
		Up_Index <- vector(,length(dIF))
		Down_Index <- vector(,length(dIF))
		condName <- paste0("dIF_", paste(ConCond[p]))
		Up_Cond <- combined$isoform_id[which(combined[[condName]] > 0.1)]
		Down_Cond <- combined$isoform_id[which(combined[[condName]] < -0.1)]

		for (i in 1:length(dIF)) {
			if (names(dIF[i]) %in% Up_Cond) {
				Up_Index[i] <- TRUE
			} else { 
				Up_Index[i] <- FALSE	
			}
			if (names(dIF[i]) %in% Down_Cond) {
				Down_Index[i] <- TRUE
			} else {
				Down_Index[i] <- FALSE
			}
		}
		# Write barcode plot to pdf
		pdf(file = file.path(mydir, "Plots", "ISAR", paste0(reference), paste0(paste(reference), "_vs_", paste(ConCond[p]), "_Barcode.pdf")), useDingbats=FALSE)

		barcodeplot(dIF, 
			index = Up_Index, 
			index2 = Down_Index, 
			main=paste0(paste(reference), " vs ", paste(ConCond[p])), 
			xlab = "dIF", 
			par(ps = 6, cex = 1, cex.lab = 1, cex.axis = 1, cex.main = 1, lwd = 0.25, pin = c(1.9685, 0.9842)))
		dev.off()
		p=p-1	
	}
	n=n-1
}

writeLines(capture.output(sessionInfo()), paste0(mydir, "/ISAR/ISAR_comparison_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
