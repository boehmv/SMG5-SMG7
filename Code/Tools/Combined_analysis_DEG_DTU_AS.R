#!/usr/bin/env Rscript

# Combined analysis of DESeq2, ISAR and Leafcutter output

# Created by: boehmv (Volker Böhm)

######
# Load libraries
######

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("readxl"))
suppressMessages(library(extrafont))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 1)

# mydir is the first argument
mydir=arguments$args[1]

####
# Transcriptome-wide changes - relaxed cut-offs
####

# Read in concatenated DGE file
DGE <- read.table(file = file.path(mydir, "DESeq2", "DGE_cat.csv"), sep = ',', header = TRUE)

# Group and summarise by cut-offs: |log2foldchange| > 1 & padj < 0.05
DGE_filt <- DGE %>% 
  group_by(condition_2, geneID, baseMean) %>% 
  summarise(DGE = any(abs(log2FoldChange) > 1 & padj < 0.05))

# Replace NA and True with 0 and 1 respectively
DGE_filt$DGE[is.na(DGE_filt$DGE)] <- 0
DGE_filt$DGE[isTRUE(DGE_filt$DGE)] <- 1 

# Read in concatenated DGE file
DTU <- read.table(file = file.path(mydir, "ISAR", "SwitchList_filt_Analyzed.csv"), sep = ',', header = TRUE)

# Extract from switchlist, group by and summarise
SwitchListTbl <- as_tibble(DTU) %>% 
  group_by(condition_2, gene_id) %>% 
  summarise(DTU = any(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05))

# Leafcutter: read in ConditionTable
df <- read.table(file = file.path(mydir, "leafcutter", "ConditionTable.txt"), sep = '\t', header = TRUE)

# Generate empty data frame for concatenated file
AS_cat <- data.frame(matrix(ncol = 3, nrow = 0))
y <- c("gene_id", "deltapsi","p.adjust")
colnames(AS_cat) <- y

# Iterate over conditions and generate data frames, generate combined df by merge
for (row in 1:nrow(df)) {
  cond=df[row, "cond"]
  path=df[row, "path"]
  
  # Read in Excel files, subset for relevant genes and columns, concatenate into AS_cat
  d <- read_excel(file.path(path))
  d_cat <- subset(d, !grepl("^SIRV", gene_id) & !grepl("^ERCC", gene_id), select = c("gene_id", "deltapsi","p.adjust"))
  d_cat$condition_2 <- paste0(cond)
  AS_cat <- rbind(AS_cat, d_cat)
}

# Produce Tibble, group by and summarise 
AS_Tbl <- as_tibble(AS_cat) %>% 
  group_by(condition_2, gene_id) %>% 
  summarise(AS = any(abs(deltapsi) > 0.1 & p.adjust < 0.05))

# Join DGE and DTU tibble together
Total_Data <- left_join(DGE_filt, SwitchListTbl, by = c("geneID" = "gene_id", "condition_2" = "condition_2"))

# Join Total_Data and AS_Tbl together
Total_Data <- left_join(Total_Data, AS_Tbl, by = c("geneID" = "gene_id", "condition_2" = "condition_2"))

# Replace NA and True with 0 and 1 respectively
Total_Data$DTU[is.na(Total_Data$DTU)] <- 0
Total_Data$AS[is.na(Total_Data$AS)] <- 0
Total_Data$DTU[isTRUE(Total_Data$DTU)] <- 1
Total_Data$AS[isTRUE(Total_Data$AS)] <- 1

# Give numbers for DTU total and fractions of all non-0 genes
Total_Data %>% group_by(condition_2) %>%  summarise(sum(DTU))
Total_Data %>% group_by(condition_2) %>% summarise(sum(DTU)/length(DTU))

# Give numbers for AS total and fractions of all non-0 genes
Total_Data %>% group_by(condition_2) %>%  summarise(sum(AS))
Total_Data %>% group_by(condition_2) %>% summarise(sum(AS)/length(AS))

# Give numbers for DGE total and fractions of all non-0 genes
Total_Data %>% group_by(condition_2) %>% summarise(sum(DGE))
Total_Data %>% group_by(condition_2) %>% summarise(sum(DGE)/length(DGE))

total_number <- Total_Data %>% group_by(condition_2) %>% summarise(length(DGE))
total_number <- as.double(total_number[1,2])

# Get classifiers
Total_Data$Class <- interaction(Total_Data$DGE, Total_Data$DTU, Total_Data$AS, lex.order = TRUE)
#Total_Data$Identity <- ifelse(Total_Data$DTU == 0 & Total_Data$DGE == 0, "none", ifelse(Total_Data$DTU == 1 & Total_Data$DGE == 0, "DTU", ifelse(Total_Data$DTU == 0 & Total_Data$DGE == 1, "DGE", ifelse(Total_Data$DTU == 1 & Total_Data$DGE == 1, "both","")))) 

# Get final table of numbers and identifiers
Total_Data_Final <- Total_Data %>% 
  group_by(condition_2, Class) %>% 
  summarise(number = n(), percentage = n()/total_number)

Total_Data_Final <- Total_Data_Final %>% group_by(condition_2, Class) %>%
  mutate(pos = cumsum(percentage) - (0.5 * percentage))

Total_Data_Final <- Total_Data_Final %>% 
  mutate(y1 = (round(percentage, 2)))

Total_Data_Final <- Total_Data_Final %>% 
  mutate(y1 = replace(y1, y1<=0.05, ""))

geom.text.size = 6
theme.size = geom.text.size / ggplot2:::.pt

# Draw overview plot
ggplot(Total_Data_Final, aes(fill=Class, y=percentage, x=condition_2, label = y1)) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(), 
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
  geom_bar(colour="black", 
           stat="identity", 
           size=0.1) +
  geom_text(aes(label = y1), 
            color = "black", 
            size = theme.size, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_brewer() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Transcriptome alterations",
       subtitle = "DGE & DTU & AS",
       x = "",
       y = "Fraction of expressed genes",
       caption = "DESeq2: |log2FC| > 1 & p.adj < 0.05 \n ISAR: |dIF| > 0.1 & & p.adj < 0.05 \n leafcutter: |dPSI| > 0.1 & & p.adj < 0.05")

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes.pdf"), width = 8, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

p1 <- ggplot(Total_Data, aes(fill = factor(Class), factor(condition_2), baseMean)) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(),
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
  scale_fill_brewer() + 
  geom_boxplot(outlier.colour = "darkgray",
               outlier.size=0.1,
               outlier.shape = 1,
               outlier.alpha = 0.5,
               outlier.stroke = 0.1,
               size=0.1,
               colour="black") +
  scale_y_continuous(trans='log10') +
  labs(fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Transcriptome alterations", 
       subtitle = "DEG & DTU & AS",
       x = "",
       y = "Mean of normalized counts",
       caption = "DESeq2: |log2FC| > 1 & p.adj < 0.05 \n ISAR: |dIF| > 0.1 & & p.adj < 0.05 \n leafcutter: |dPSI| > 0.1 & & p.adj < 0.05")

p1

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Boxplot_baseMean.pdf"), width = 12, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")


# Prepare cumulative plot
dat <- Total_Data %>% 
  group_by(condition_2) %>% 
  arrange(baseMean, .by_group = TRUE) %>%
  mutate(x = seq(1, length(baseMean)))

ggplot(dat %>%
         arrange(Class), aes(y=baseMean, x = x, color = Class)) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(),
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
  scale_color_viridis(discrete=TRUE) +
  geom_point(alpha = 0.2) + 
  facet_wrap(~ condition_2) +
  scale_y_continuous(trans='log10') +
  labs(title = "Transcriptome alterations", 
       subtitle = "DGE & DTU & AS",
       x = "Index of expressed genes",
       y = "Mean of normalized counts", 
       caption = "DESeq2: |log2FC| > 1 & p.adj < 0.05 \n ISAR: |dIF| > 0.1 & & p.adj < 0.05 \n leafcutter: |dPSI| > 0.1 & & p.adj < 0.05")


# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Cumulative.pdf"), width = 12, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

####
# Transcriptome-wide changes - stringent cut-offs
####

# Group and summarise by cut-offs: |log2foldchange| > 1 & padj < 0.05
DGE_filt_stringent <- DGE %>% 
  group_by(condition_2, geneID, baseMean) %>% 
  summarise(DGE = any(abs(log2FoldChange) > 2 & padj < 0.01))

# Replace NA and True with 0 and 1 respectively
DGE_filt_stringent$DGE[is.na(DGE_filt_stringent$DGE)] <- 0
DGE_filt_stringent$DGE[isTRUE(DGE_filt_stringent$DGE)] <- 1 

# Extract from switchlist, group by and summarise
SwitchListTbl_stringent <- as_tibble(DTU) %>% 
  group_by(condition_2, gene_id) %>% 
  summarise(DTU = any(abs(dIF) > 0.2 & isoform_switch_q_value < 0.01))

# Generate empty data frame for concatenated file
AS_cat <- data.frame(matrix(ncol = 3, nrow = 0))
y <- c("gene_id", "deltapsi","p.adjust")
colnames(AS_cat) <- y

# Iterate over conditions and generate data frames, generate combined df by merge
for (row in 1:nrow(df)) {
  cond=df[row, "cond"]
  path=df[row, "path"]
  
  # Read in Excel files, subset for relevant genes and columns, concatenate into AS_cat
  d <- read_excel(file.path(path))
  d_cat <- subset(d, !grepl("^SIRV", gene_id) & !grepl("^ERCC", gene_id), select = c("gene_id", "deltapsi","p.adjust"))
  d_cat$condition_2 <- paste0(cond)
  AS_cat <- rbind(AS_cat, d_cat)
}

# Produce Tibble, group by and summarise 
AS_Tbl_stringent <- as_tibble(AS_cat) %>% 
  group_by(condition_2, gene_id) %>% 
  summarise(AS = any(abs(deltapsi) > 0.2 & p.adjust < 0.01))

# Join DGE and DTU tibble together
Total_Data_stringent <- left_join(DGE_filt_stringent, SwitchListTbl_stringent, by = c("geneID" = "gene_id", "condition_2" = "condition_2"))

# Join Total_Data and AS_Tbl together
Total_Data_stringent <- left_join(Total_Data_stringent, AS_Tbl_stringent, by = c("geneID" = "gene_id", "condition_2" = "condition_2"))

# Replace NA and True with 0 and 1 respectively
Total_Data_stringent$DTU[is.na(Total_Data_stringent$DTU)] <- 0
Total_Data_stringent$AS[is.na(Total_Data_stringent$AS)] <- 0
Total_Data_stringent$DTU[isTRUE(Total_Data_stringent$DTU)] <- 1
Total_Data_stringent$AS[isTRUE(Total_Data_stringent$AS)] <- 1

# Give numbers for DTU total and fractions of all non-0 genes
Total_Data_stringent %>% group_by(condition_2) %>%  summarise(sum(DTU))
Total_Data_stringent %>% group_by(condition_2) %>% summarise(sum(DTU)/length(DTU))

# Give numbers for AS total and fractions of all non-0 genes
Total_Data_stringent %>% group_by(condition_2) %>%  summarise(sum(AS))
Total_Data_stringent %>% group_by(condition_2) %>% summarise(sum(AS)/length(AS))

# Give numbers for DGE total and fractions of all non-0 genes
Total_Data_stringent %>% group_by(condition_2) %>% summarise(sum(DGE))
Total_Data_stringent %>% group_by(condition_2) %>% summarise(sum(DGE)/length(DGE))

total_number_stringent <- Total_Data_stringent %>% group_by(condition_2) %>% summarise(length(DGE))
total_number_stringent <- as.double(total_number_stringent[1,2])

# Get classifiers
Total_Data_stringent$Class <- interaction(Total_Data_stringent$DGE, Total_Data_stringent$DTU, Total_Data_stringent$AS, lex.order = TRUE)
#Total_Data$Identity <- ifelse(Total_Data$DTU == 0 & Total_Data$DGE == 0, "none", ifelse(Total_Data$DTU == 1 & Total_Data$DGE == 0, "DTU", ifelse(Total_Data$DTU == 0 & Total_Data$DGE == 1, "DGE", ifelse(Total_Data$DTU == 1 & Total_Data$DGE == 1, "both","")))) 

# Get final table of numbers and identifiers
Total_Data_stringent_Final <- Total_Data_stringent %>% 
  group_by(condition_2, Class) %>% 
  summarise(number = n(), percentage = n()/total_number)

Total_Data_stringent_Final <- Total_Data_stringent_Final %>% group_by(condition_2, Class) %>%
  mutate(pos = cumsum(percentage) - (0.5 * percentage))

Total_Data_stringent_Final <- Total_Data_stringent_Final %>% 
  mutate(y1 = (round(percentage, 2)))

Total_Data_Final <- Total_Data_Final %>% 
  mutate(y1 = replace(y1, y1<=0.05, ""))

geom.text.size = 6
theme.size = geom.text.size / (14/5)

# Draw overview plot
ggplot(Total_Data_stringent_Final, aes(fill=Class, y=percentage, x=condition_2, label = y1)) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(), 
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
  geom_bar(colour="black", 
           stat="identity", 
           size=0.1) +
  geom_text(aes(label = y1), 
            color = "black", 
            size = theme.size, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_brewer() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Transcriptome alterations",
       subtitle = "DGE & DTU & AS",
       x = "",
       y = "Fraction of expressed genes",
       caption = "DESeq2: |log2FC| > 2 & p.adj < 0.01 \n ISAR: |dIF| > 0.2 & & p.adj < 0.01 \n leafcutter: |dPSI| > 0.2 & & p.adj < 0.01")

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_stringent.pdf"), width = 8, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")


p2 <- ggplot(Total_Data_stringent, aes(fill = factor(Class), factor(condition_2), baseMean)) + 
  theme_classic() +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(),
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
  scale_fill_brewer() + 
  geom_boxplot(outlier.colour = "darkgray",
               outlier.size=0.1,
               outlier.shape = 1,
               outlier.alpha = 0.5,
               outlier.stroke = 0.1,
               size=0.1,
               colour="black") +
  scale_y_continuous(trans='log10') +
  labs(fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Transcriptome alterations", 
       subtitle = "DEG & DTU & AS",
       x = "",
       y = "Mean of normalized counts",
       caption = "DESeq2: |log2FC| > 2 & p.adj < 0.01 \n ISAR: |dIF| > 0.2 & & p.adj < 0.01 \n leafcutter: |dPSI| > 0.2 & & p.adj < 0.01")

p2

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Boxplot_baseMean_stringent.pdf"), width = 12, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

# Pheatmap
Total_Data_stringent_Final_Matrix <- Total_Data_stringent_Final %>% 
  select(condition_2, Class, y1) %>%
  pivot_wider(names_from = condition_2, values_from = y1) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames(var = "Class")

pheatmap(as.matrix(Total_Data_stringent_Final_Matrix), color = inferno(50), display_numbers = as.matrix(Total_Data_stringent_Final_Matrix))

# Prepare cumulative plot
dat_stringent <- Total_Data_stringent %>% 
  group_by(condition_2) %>% 
  arrange(baseMean, .by_group = TRUE) %>%
  mutate(x = seq(1, length(baseMean)))

dat_stringent %>% 
  summarise(quantile = scales::percent(c(0, 0.25, 0.5, 0.75, 1)),
            baseMean = quantile(baseMean, c(0, 0.25, 0.5, 0.75, 1)))

ggplot(dat_stringent %>%
         arrange(Class), aes(y=baseMean, x = x, color = Class)) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(),
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
  scale_color_viridis(discrete=TRUE) +
  geom_point(alpha = 0.2) + 
  facet_wrap(~ condition_2) +
  scale_y_continuous(trans='log10') +
  labs(title = "Transcriptome alterations", 
       subtitle = "DGE & DTU & AS",
       x = "Index of expressed genes",
       y = "Mean of normalized counts",
       caption = "DESeq2: |log2FC| > 2 & p.adj < 0.01 \n ISAR: |dIF| > 0.2 & & p.adj < 0.01 \n leafcutter: |dPSI| > 0.2 & & p.adj < 0.01")

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Cumulative_Stringent.pdf"), width = 12, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

# Arrange in grid
g <- grid.arrange(p1, p2, nrow = 2)

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Boxplot_combined.pdf"), g, width = 12, height=14, device=cairo_pdf, units = c("cm"), bg = "transparent")

# Violin plot of baseMean
p3 <- ggplot(Total_Data_stringent, aes(factor(condition_2), baseMean))+
  theme_classic() +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        panel.background = element_blank(),
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
  scale_fill_brewer() +
  geom_violin(draw_quantiles = c(0, 0.33, 0.66, 1), scale = "count") +
  scale_y_continuous(trans='log10') +
  labs(fill = "Class") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Transcriptome alterations", 
       subtitle = "DGE & DTU & AS",
       x = "", 
       y = "Mean of normalized counts",
       caption = "DESeq2: |log2FC| > 2 & p.adj < 0.01 \n ISAR: |dIF| > 0.2 & & p.adj < 0.01 \n leafcutter: |dPSI| > 0.2 & & p.adj < 0.01")

# Save to disk
ggsave(file.path(mydir, "Plots", "Transcriptome_Changes_Mean.pdf"), p3, width = 12, height=8, device=cairo_pdf, units = c("cm"), bg = "transparent")

writeLines(capture.output(sessionInfo()), paste0(mydir, "/Plots/Combined_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
