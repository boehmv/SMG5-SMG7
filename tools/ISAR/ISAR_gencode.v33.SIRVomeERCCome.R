#!/usr/bin/env Rscript


# Title: ISAR_gencode.v33.SIRVomeERCCome
# Objective: Standard IsoformSwitchAnalyseR (ISAR) pipeline giving differential transcript usage (DTU) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(IsoformSwitchAnalyzeR)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition string is the second, gets converted to character
cond=arguments$args[2]
cond <- unlist(strsplit(cond,","))

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Get design table
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)
myDesign <- data.frame(
    sampleID = samples$sample,
    condition = samples$condition
)

# Compile sampleVector from experiment file
myFiles <- file.path(mydir, "Salmon", samples$sample, "quant.sf")
names(myFiles) <- samples$sample

# Import quantifications
salmonQuant <- importIsoformExpression(
    sampleVector = myFiles,
    addIsofomIdAsColumn = TRUE
)

# Get comparison data frame
myComparison <- data.frame(
    condition_1 = "control",
    condition_2 = cond[-1]
)

# Generate switchlist
SwitchList <- importRdata(
		isoformCountMatrix   = salmonQuant$counts,
		isoformRepExpression = salmonQuant$abundance,
		designMatrix         = myDesign,
		addAnnotatedORFs     = TRUE,
		onlyConsiderFullORF = FALSE,
		removeNonConvensionalChr = TRUE,
		isoformExonAnnoation = file.path(ref_dir, "Gencode",  "gencode.v33.SIRVomeERCCome.annotation.gtf"),
		comparisonsToMake= myComparison,	    
		showProgress = TRUE,
		ignoreAfterBar = TRUE,
    		ignoreAfterSpace = TRUE,
    		ignoreAfterPeriod = FALSE,	# Set to FALSE for Gencode
    		removeTECgenes = TRUE		# If set to TRUE, spike_ins need "gene_name" and "gene_type"
)

# What is in the Switchlist
SwitchList

# Filtering
SwitchList_filt <- preFilter(
  SwitchList,
  geneExpressionCutoff = 1, # FPMK threshold
  isoformExpressionCutoff = 0, # FPMK threshold
  IFcutoff=0.01,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  quiet=FALSE
)

# Testing for Isoform Switches via DEXSeq
SwitchList_filt_Analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchList_filt,
    reduceToSwitchingGenes=TRUE
)

# Set WD
setwd((file.path(mydir, "ISAR")))

save.image(file='ISAR_session.RData')

# Summarize switching features
extractSwitchSummary(SwitchList_filt_Analyzed)

# Predicting Switch Consequences
SwitchList_filt_Analyzed <- analyzeSwitchConsequences(
    SwitchList_filt_Analyzed,
    consequencesToAnalyze = c('NMD_status'), 
    dIFcutoff = 0.1, 
    showProgress=TRUE
)

# Summarize switching features without consequences
extractSwitchSummary(SwitchList_filt_Analyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

# Summarize switching features with consequences
extractSwitchSummary(SwitchList_filt_Analyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

# Plot Top 10 Switches
switchPlotTopSwitches(
    switchAnalyzeRlist = SwitchList_filt_Analyzed, 
    n = 10,
    pathToOutput = file.path(mydir, "ISAR","Plots"),
    filterForConsequences = FALSE, 
    splitFunctionalConsequences = TRUE
)

# Summarize switching features with consequences
extractTopSwitches(
    SwitchList_filt_Analyzed, 
    filterForConsequences = TRUE, 
    n = 10, 
    sortByQvals = TRUE
)

# Volcano plot
pdf("Volcano_plot.pdf")
ggplot(data=SwitchList_filt_Analyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('gray','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

pdf("Volcano_plot_PTC.pdf")
ggplot(data=SwitchList_filt_Analyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=PTC & isoform_switch_q_value < 0.05 & abs(dIF) > 0.1 ), # default cutoff
        size=0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('PTC', values = c('gray','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

write.csv(SwitchList_filt_Analyzed$isoformFeatures, file="SwitchList_filt_Analyzed.csv")
save.image(file='ISAR_session.RData')

writeLines(capture.output(sessionInfo()), paste0(mydir, "/ISAR/ISAR_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
