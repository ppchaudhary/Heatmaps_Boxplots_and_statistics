# Load required libraries
library(ggplot2)
library(pheatmap)
library(dplyr)

# Read in count matrix and metadata
counts <- read.csv("Counts data_for use.csv", row.names = 1)
metadata <- read.csv("Malemetadata.csv")

# Ensure 'Week' is a factor and matches the correct order
metadata$Week <- factor(metadata$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Ensure metadata rows match count matrix columns
metadata <- metadata[match(colnames(counts), metadata$SampleName),]

# Load gene symbols mapping
gene_symbols <- read.csv("GeneSymbols.csv")
rownames(gene_symbols) <- gene_symbols$GeneID

# Define genes of interest
genes_of_interest <- c("Trpa1", "Trpc1", "Trpc2", "Trpc3", "Trpc4", "Trpc5", "Trpc6", "Trpc7",
                       "Trpm1", "Trpm2", "Trpm3", "Trpm4", "Trpm5", "Trpm6", "Trpm7", "Trpm8",
                       "Trpv1", "Trpv2", "Trpv3", "Trpv4", "Trpv5", "Trpv6", "Il4ra", "Il13",
                       "Il17a", "Il17b", "Il17c", "Il17d", "Il17f", "Il17ra", "Il17rb", "Il17rc",
                       "Il17rd", "Il17re", "Ifng", "Tnf")

# Map to Gene IDs
gene_ids_of_interest <- gene_symbols %>%
  filter(GeneSymbol %in% genes_of_interest) %>%
  pull(GeneID)

# Filter the count matrix to include only genes of interest
filtered_counts <- counts[rownames(counts) %in% gene_ids_of_interest, ]

# Scale data by row (z-score normalization)
heatmap_data <- t(scale(t(as.matrix(filtered_counts))))

# Remove rows with NA/NaN/Inf (usually from genes with no variance)
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

# Replace rownames with gene symbols
row_names_with_symbols <- gene_symbols[rownames(heatmap_data), "GeneSymbol"]
rownames(heatmap_data) <- row_names_with_symbols

# Prepare column annotations (must match column names of heatmap_data)
annotation_col <- data.frame(Week = metadata$Week)
rownames(annotation_col) <- metadata$SampleName

# Plot heatmap and save to file
pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         fontsize = 8,
         legend = TRUE,
         legend_title = "Z-score Expression",
         filename = "heatmap_selected_genes.png")


# Save scaled expression data to CSV
write.csv(heatmap_data, file = "heatmap_selected_genes_data.csv")




###Box plots ,kruskal wallis and dunns posthoc test

# Load libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(FSA)  # for dunnTest
library(rstatix)  # optional, for pairwise comparisons

# Read in data
counts <- read.csv("Counts data_for use.csv", row.names = 1)
metadata <- read.csv("Malemetadata.csv")
gene_symbols <- read.csv("GeneSymbols.csv")
rownames(gene_symbols) <- gene_symbols$GeneID

# Ensure metadata ordering and factor conversion
metadata <- metadata[match(colnames(counts), metadata$SampleName), ]
metadata$Week <- factor(metadata$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Create DESeq2 object and normalize
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Week)
dds <- estimateSizeFactors(dds)
expr_data <- as.data.frame(counts(dds, normalized = TRUE))

# Genes of interest
genes_of_interest <- c("Trpa1", "Trpc1", "Trpc2", "Trpc3", "Trpc4", "Trpc5", "Trpc6", "Trpc7",
                       "Trpm1", "Trpm2", "Trpm3", "Trpm4", "Trpm5", "Trpm6", "Trpm7", "Trpm8",
                       "Trpv1", "Trpv2", "Trpv3", "Trpv4", "Trpv5", "Trpv6", "Il4ra", "Il13",
                       "Il17a", "Il17b", "Il17c", "Il17d", "Il17f", "Il17ra", "Il17rb", "Il17rc",
                       "Il17rd", "Il17re", "Ifng", "Tnf")

# Map symbols to IDs
gene_ids_of_interest <- gene_symbols %>%
  filter(GeneSymbol %in% genes_of_interest) %>%
  pull(GeneID)

# Filter expression data
expr_filtered <- expr_data[rownames(expr_data) %in% gene_ids_of_interest, ]
expr_filtered$GeneSymbol <- gene_symbols[rownames(expr_filtered), "GeneSymbol"]

# Reshape for plotting
expr_long <- pivot_longer(expr_filtered, cols = -GeneSymbol, names_to = "Sample", values_to = "Expression")
expr_long <- left_join(expr_long, metadata[, c("SampleName", "Week")], by = c("Sample" = "SampleName"))
expr_long$Week <- factor(expr_long$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Directory for saving plots
dir.create("Gene_Boxplots", showWarnings = FALSE)

# Initialize p-value summary
pval_summary <- data.frame()

# Loop through each gene
for (gene in unique(expr_long$GeneSymbol)) {
  gene_data <- expr_long %>% filter(GeneSymbol == gene)
  
  # Check if Expression has any valid data (NA or Inf handling)
  if (all(is.na(gene_data$Expression)) || all(gene_data$Expression == 0)) {
    message(paste("Skipping gene due to empty expression:", gene))
    next  # Skip this gene
  }
  
  # Calculate y_max for the boxplot (ensure valid value)
  y_max <- max(gene_data$Expression, na.rm = TRUE)
  
  if (is.finite(y_max) && !is.na(y_max)) {
    y_label_pos <- y_max * 1.1
  } else {
    y_label_pos <- 1  # Fallback position if max is NA or Inf
  }
  
  # Create the boxplot
  p <- ggboxplot(gene_data, x = "Week", y = "Expression", fill = "Week",
                 palette = "jco", add = "jitter", title = paste("Expression of", gene)) +
    stat_compare_means(method = "kruskal.test", label.y = y_label_pos)
  
  # Save boxplot
  ggsave(filename = paste0("Gene_Boxplots/", gene, "_boxplot.png"), plot = p, width = 6, height = 4)
  
  # Save Dunn's test results (if applicable)
  if (exists("dunn_df") && !is.null(dunn_df) && nrow(dunn_df) > 0) {
    write.csv(dunn_df, file = paste0("Gene_Boxplots/", gene, "_dunn_test.csv"), row.names = FALSE)
  }
}

####Selective boxplots and with posthoc test values

# Load libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Read in data
counts <- read.csv("Counts data_for use.csv", row.names = 1)
metadata <- read.csv("Malemetadata.csv")
gene_symbols <- read.csv("GeneSymbols.csv")
rownames(gene_symbols) <- gene_symbols$GeneID

# Ensure metadata ordering and factor conversion
metadata <- metadata[match(colnames(counts), metadata$SampleName), ]
metadata$Week <- factor(metadata$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Create DESeq2 object and normalize
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Week)
dds <- estimateSizeFactors(dds)
expr_data <- as.data.frame(counts(dds, normalized = TRUE))

# Genes of interest
genes_of_interest <- c("Trpa1", "Trpc1", "Trpc2", "Trpc3", "Trpc4", "Trpc5", "Trpc6", "Trpc7",
                       "Trpm1", "Trpm2", "Trpm3", "Trpm4", "Trpm5", "Trpm6", "Trpm7", "Trpm8",
                       "Trpv1", "Trpv2", "Trpv3", "Trpv4", "Trpv5", "Trpv6", "Il4ra", "Il13",
                       "Il17a", "Il17b", "Il17c", "Il17d", "Il17f", "Il17ra", "Il17rb", "Il17rc",
                       "Il17rd", "Il17re", "Ifng", "Tnf")

# Map symbols to IDs
gene_ids_of_interest <- gene_symbols %>%
  filter(GeneSymbol %in% genes_of_interest) %>%
  pull(GeneID)

# Filter expression data
expr_filtered <- expr_data[rownames(expr_data) %in% gene_ids_of_interest, ]
expr_filtered$GeneSymbol <- gene_symbols[rownames(expr_filtered), "GeneSymbol"]

# Reshape for plotting
expr_long <- pivot_longer(expr_filtered, cols = -GeneSymbol, names_to = "Sample", values_to = "Expression")
expr_long <- left_join(expr_long, metadata[, c("SampleName", "Week")], by = c("Sample" = "SampleName"))
expr_long$Week <- factor(expr_long$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Directory for saving plots
dir.create("Gene_Boxplots", showWarnings = FALSE)

# Loop through each gene
for (gene in unique(expr_long$GeneSymbol)) {
  gene_data <- expr_long %>% filter(GeneSymbol == gene)
  
  # Check if Expression has any valid data (NA or Inf handling)
  if (all(is.na(gene_data$Expression)) || all(gene_data$Expression == 0)) {
    message(paste("Skipping gene due to empty expression:", gene))
    next  # Skip this gene
  }
  
  # Perform Kruskal-Wallis test
  kw_result <- kruskal.test(Expression ~ Week, data = gene_data)
  if (kw_result$p.value <= 0.05) {
    message(paste("Kruskal-Wallis result for gene", gene, "p-value:", kw_result$p.value))
    
    # Calculate y_max for the boxplot (ensure valid value)
    y_max <- max(gene_data$Expression, na.rm = TRUE)
    y_label_pos <- ifelse(is.finite(y_max), y_max * 1.1, 1)  # Fallback position if max is NA or Inf
    
    # Create the boxplot
    p <- ggboxplot(gene_data, x = "Week", y = "Expression", fill = "Week",
                   palette = "jco", add = "jitter", title = paste("Expression of", gene)) +
      stat_compare_means(method = "kruskal.test", label.y = y_label_pos)
    
    # Save boxplot
    ggsave(filename = paste0("Gene_Boxplots_posthoc/", gene, "_boxplot.png"), plot = p, width = 6, height = 4)
  } else {
    message(paste("Gene", gene, "did not pass Kruskal-Wallis test (p >", 0.05, "), skipping plot."))
  }
}

















