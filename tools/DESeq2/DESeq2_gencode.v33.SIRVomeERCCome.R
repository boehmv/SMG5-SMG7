#!/usr/bin/env Rscript

# Title: DESeq2_gencode.v33.SIRVomeERCCome
# Objective: Standard DESeq2 pipeline giving differentially expressed genes (DGE) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######
suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(readr)
  library(DESeq2)
  library(Glimma)
  library(limma)
  library(pheatmap)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(gprofiler2)
  library(data.table)
  library(ggplot2)
})

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition string is the second, gets converted to character
cond=arguments$args[2]
cond <- unlist(strsplit(cond,","))

# Survey name
myname <- tail(unlist(strsplit(mydir,"/")), n=1)

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Create "Plots" folder if it does not exists
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)

# Create "DESeq2" folder if it does not exists
dir.create(file.path(mydir, "Plots", "DESeq2"), showWarnings = FALSE)

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
all(file.exists(files))

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v33.SIRV.ERCC.tsv"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene,  ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Filter out all low expressed genes with less than 10 counts
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep,]

# Perform the DESeq analysis
dds <- DESeq(ddsTxi)

######
# Sample-to-sample distance
######

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Generate heatmap of the sample-to-sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, filename=(file.path(mydir, "Plots", "DESeq2", "sampleDistanceHeatmap.pdf")))

# Principal component plot of the samples
plt <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(plt, "percentVar"))

ggplot(plt, aes(PC1, PC2, color=condition)) +
  theme_classic() +	
  geom_point(size=3) +		
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10)) +
  labs(title = paste(myname), subtitle = "DESeq2 PCA analysis", caption = "DESeq2 (1.24.0)", x = paste0("PC1: ",percentVar[1],"% variance"), y = paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggsave(file.path(mydir, "Plots", "DESeq2", "PCA_plot.pdf"), width = 6, height=4, device=cairo_pdf)

# Generate Glimma MDS plot
glMDSPlot(dds, groups=dds$condition, path=(file.path(mydir, "Plots", "DESeq2")), folder="glimma-MDS", launch=FALSE)

#######
# Individual analyses
#######

# Build for-loop to go through all conditions vs control
mycount=2
for (i in 2:length(cond)){
  # Perform DESeq2
  setwd((file.path(mydir, "DESeq2")))
  condition <- paste(cond[mycount])
  res <- results(dds, alpha=0.05, contrast=c("condition",paste(cond[mycount]),"control"))	
  
  # Get additional IDs (remove transcript version numbers to get other identifiers)
  ensIDs <- sub("\\..*", "", row.names(res))	
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys=ensIDs,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys=ensIDs,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$uniprot <- mapIds(org.Hs.eg.db,
                        keys=ensIDs,
                        column="UNIPROT",
                        keytype="ENSEMBL",
                        multiVals="first")
  # Convert results file to data.frame
  resDF <- as.data.frame(res)
  # Create condition-specific subfolder
  dir.create(file.path(mydir, "DESeq2", paste(cond[mycount])), showWarnings = FALSE)
  # Modify data frame column headers
  resDF_out <- setDT(resDF, keep.rownames = "geneID")[]
  resDF_out[,condition] <- NA
  # Generate output csv file
  message("Generate csv output file for ", paste(cond[mycount]))
  write.csv(resDF_out, 
            file=(file.path(mydir, "DESeq2", cond[mycount] ,paste0(paste(cond[mycount]),"_vs_control_DESeq2_results.csv"))))
  # Make Combined folder
  dir.create(file.path(mydir, "DESeq2", "Combined"), showWarnings = FALSE)
  # Create symbolic link if needed
  if(file.exists(file.path(mydir, "DESeq2", "Combined", paste0(paste(cond[mycount]),"_vs_control_DESeq2_results.csv")))){
    print("Symlink already established")
  } else {
    file.symlink(file.path(mydir, "DESeq2", cond[mycount] ,paste0(paste(cond[mycount]),"_vs_control_DESeq2_results.csv")), file.path(mydir, "DESeq2", "Combined", paste0(paste(cond[mycount]),"_vs_control_DESeq2_results.csv")))
  }
  # Generate Glimma MD plot
  message("Generate Glimma MD plots for ", paste(cond[mycount]))
  status <- as.numeric(res$padj < .05)
  anno <- data.frame(GeneID=rownames(res), symbol=res$symbol, UNIPROT=res$uniprot)
  glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE),
           groups=dds$condition, transform=TRUE,
           samples=colnames(dds), anno=anno,
           path=(file.path(mydir, "DESeq2", paste(cond[mycount]))), folder="glimma-MD", launch=FALSE)
  
  #######
  # Gene ontology analysis via g:profiler
  #######
  
  # Create condition-specific subfolder
  dir.create(file.path(mydir, "DESeq2", paste(cond[mycount]), "GO"), showWarnings = FALSE)
  setwd(file.path(mydir, "DESeq2", paste(cond[mycount]), "GO"))
  
  message("GO analysis for ", paste(cond[mycount]), " running")
  
  # Get background vector of ensembl gene IDs (basically all non-zero genes), strip transcript version (incompatible with g:profiler)
  background_genes <- sub("\\..*", "", resDF[['geneID']])
  
  # Filter DESeq2 results (padj < 0.05 and |log2FC| > 1)
  resDF_filt <- resDF %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  # Prepare ordered list (by padj) for g:profiler analysis
  resDF_ordered <- resDF_filt[order(resDF_filt$padj),]
  # Get ordered query vector of ensembl gene IDs
  padj_query_genes <- sub("\\..*", "", resDF_ordered[['geneID']])
  
  # Prepare ordered list (by upregulated) for g:profiler analysis
  resDF_ordered <- resDF_filt[order(resDF_filt$log2FoldChange, decreasing = TRUE),]
  resDF_ordered <- resDF_ordered %>% filter(padj < 0.05 & log2FoldChange > 1)
  # Get ordered query vector of ensembl gene IDs
  Up_query_genes <- sub("\\..*", "", resDF_ordered[['geneID']])
  
  # Prepare ordered list (by downregulated) for g:profiler analysis
  resDF_ordered <- resDF_filt[order(resDF_filt$log2FoldChange, decreasing = FALSE),]
  resDF_ordered <- resDF_ordered %>% filter(padj < 0.05 & log2FoldChange < -1)
  # Get ordered query vector of ensembl gene IDs
  Down_query_genes <- sub("\\..*", "", resDF_ordered[['geneID']])
  
  # Run the gene list functional enrichment from gprofiler2 using tryCatch (for robust error handling)
  tryCatch(
    
    expr = {
      
      # Perform the GO analysis
      gostres <- gost(query = list("padj" = padj_query_genes, "Up" = Up_query_genes, "Down" = Down_query_genes), multi_query = TRUE, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, exclude_iea = TRUE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "g_SCS", custom_bg = background_genes, sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "HP"))
      
      # Produce static plot
      static_plot <- gostplot(gostres, capped = FALSE, interactive = FALSE, pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG = "#dd4477", REAC = "#3366cc", HP = "#990099"))
      
      # Publish the static plot
      publish_gostplot(static_plot, highlight_terms = NULL, filename = (paste0(paste(cond[mycount]), "_vs_control_gprofiler.pdf")))
      
      # Produce interactive plot
      interactive_plot <- gostplot(gostres, capped = FALSE, interactive = TRUE, pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG = "#dd4477", REAC = "#3366cc", HP = "#990099"))
      
      # Publish the interactive plot
      htmlwidgets::saveWidget(widget=interactive_plot, (paste0(paste(cond[mycount]), "_vs_control_gprofiler.html")), title = paste(cond[mycount]))
      
      gprofiler_table <- subset(gostres$result, select = c("query", "significant", "p_value","term_id", "source", "term_name"))
      
      # Produce g:profiler output table
      #publish_gosttable(gostres, highlight_terms = NULL, use_colors = TRUE, show_columns = c("source", "term_name", "term_size", "intersection_size"), ggplot = TRUE)
      
      message("GO analysis for ", paste(cond[mycount]), " successful")
    },
    error = function(e){ 
      message("ERROR - ", "GO analysis for ", paste(cond[mycount]), " NOT successful!")
      # Do this if an error is caught...
    },
    warning = function(w){
      message("WARNING - ", "GO analysis for ", paste(cond[mycount]), " caught an warning!")
      # Do this if an warning is caught...
    },
    finally = {
      message("GO analysis for ", paste(cond[mycount]), " finished")
    }
  )
  
  # Save R session to file for later	
  message("Saved DESeq2 analysis RData for ", paste(cond[mycount]), " in ", file.path(mydir, "DESeq2", paste(cond[mycount])))
  setwd(file.path(mydir, "DESeq2", paste(cond[mycount])))
  save.image(file='DESeq2_session.RData')
  
  # Counter increment
  mycount = mycount + 1
}

writeLines(capture.output(sessionInfo()), paste0(mydir, "/DESeq2/DESeq2_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
