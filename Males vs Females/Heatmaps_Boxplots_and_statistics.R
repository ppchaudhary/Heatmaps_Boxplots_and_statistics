# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Create output directory if it doesn't exist
output_dir <- "better_annotations_and_stats"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read in count matrix and metadata
counts <- read.csv("Counts data_for use.csv", row.names = 1)
metadata <- read.csv("Metadata.csv")

# Ensure 'Week' is a factor with ordered levels
metadata$Week <- factor(metadata$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Ensure metadata rows are in the same order as count matrix columns
metadata <- metadata[match(colnames(counts), metadata$SampleName), ]

# Load gene symbols mapping and set rownames to GeneID
gene_symbols <- read.csv("GeneSymbols.csv")
rownames(gene_symbols) <- gene_symbols$GeneID

# Convert count matrix values to integers (round)
counts[] <- lapply(counts, function(x) as.integer(round(x)))

# Filter out lowly expressed genes (mean count < 10)
row_means <- rowMeans(counts)
counts <- counts[row_means >= 10, ]

# Define list of specific genes of interest
specific_genes <- c("Trpa1", "Trpc1", "Trpc2", "Trpc3", "Trpc4", "Trpc5", "Trpc6", "Trpc7",
                    "Trpm1", "Trpm2", "Trpm3", "Trpm4", "Trpm5", "Trpm6", "Trpm7", "Trpm8",
                    "Trpv1", "Trpv2", "Trpv3", "Trpv4", "Trpv5", "Trpv6", "Il4ra", "Il13",
                    "Il17a", "Il17b", "Il17c", "Il17d", "Il17f", "Il17ra", "Il17rb", "Il17rc",
                    "Il17rd", "Il17re", "Ifng", "Tnf")

# Map to Gene IDs
specific_gene_ids <- gene_symbols$GeneID[gene_symbols$GeneSymbol %in% specific_genes]

# Loop through each week to generate heatmaps
for (week in levels(metadata$Week)) {
  
  # Subset metadata and counts for current week
  metadata_week <- metadata[metadata$Week == week, ]
  counts_week <- counts[, colnames(counts) %in% metadata_week$SampleName, drop = FALSE]
  
  # Identify which of the specific genes are present
  specific_genes_present <- intersect(rownames(counts_week), specific_gene_ids)
  
  # Proceed only if at least one gene is found
  if (length(specific_genes_present) > 0) {
    
    # Extract and scale expression data for specific genes
    specific_heatmap_data <- counts_week[specific_genes_present, ]
    specific_heatmap_data <- t(scale(t(specific_heatmap_data)))  # Z-score normalization
    
    # Remove any rows with NA/NaN/Inf (usually due to no variance)
    specific_heatmap_data <- specific_heatmap_data[complete.cases(specific_heatmap_data), ]
    
    # Replace row names with gene symbols
    row_names_with_symbols <- gene_symbols[rownames(specific_heatmap_data), "GeneSymbol"]
    rownames(specific_heatmap_data) <- row_names_with_symbols
    
    # Create column annotation for the heatmap
    annotation_col <- data.frame(Gender = metadata_week$Gender)
    rownames(annotation_col) <- metadata_week$SampleName
    
    # Construct output filename
    file_name <- file.path(output_dir, paste0("Specific_Genes_Male_vs_Female_", week, "_heatmap.png"))
    
    # Plot and save the heatmap
    pheatmap(specific_heatmap_data,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             annotation_col = annotation_col,
             fontsize = 8,
             legend = TRUE,
             legend_title = "Z-score Expression",
             filename = file_name)
  }
}


#####Box plots and statistics

# Load libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(FSA)
library(rstatix)

# Read in data
counts <- read.csv("Counts data_for use.csv", row.names = 1)
metadata <- read.csv("Metadata.csv")
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
expr_long <- left_join(expr_long, metadata[, c("SampleName", "Week", "Gender")], by = c("Sample" = "SampleName"))
expr_long$Week <- factor(expr_long$Week, levels = c("Week2", "Week4", "Week8", "Week16", "Week30"))

# Create output directory
dir.create("better_annotations_and_stats/boxplots", recursive = TRUE, showWarnings = FALSE)

# Loop through genes
for (gene in unique(expr_long$GeneSymbol)) {
  gene_data <- expr_long %>% filter(GeneSymbol == gene)
  
  if (all(is.na(gene_data$Expression)) || all(gene_data$Expression == 0)) {
    message(paste("Skipping gene due to empty expression:", gene))
    next
  }
  
  # Run Wilcoxon tests per week only if both groups are valid
  wilcox_results <- data.frame()
  for (week in unique(gene_data$Week)) {
    week_data <- gene_data %>% filter(Week == week)
    if (nrow(week_data) < 3) next
    
    group_counts <- table(week_data$Gender)
    if (length(group_counts) == 2 && all(group_counts > 1)) {
      test_result <- tryCatch({
        week_data %>%
          wilcox_test(Expression ~ Gender) %>%
          mutate(Week = week)
      }, error = function(e) NULL)
      
      if (!is.null(test_result)) {
        wilcox_results <- bind_rows(wilcox_results, test_result)
      }
    }
  }
  
  # Adjust p-values
  if (nrow(wilcox_results) > 0) {
    wilcox_results <- wilcox_results %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
  }
  
  # Plot
  y_max <- max(gene_data$Expression, na.rm = TRUE)
  y_label_pos <- ifelse(is.finite(y_max), y_max * 1.1, 1)
  
  p <- ggboxplot(gene_data, x = "Gender", y = "Expression", fill = "Gender",
                 palette = "npg", add = "jitter", facet.by = "Week") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.y = y_label_pos) +
    labs(title = paste("Expression of", gene, "in Males vs Females across Weeks"),
         y = "Normalized Expression", x = "Gender") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = paste0("better_annotations_and_stats/boxplots/", gene, "_combined_boxplot.png"),
         plot = p, width = 12, height = 6)
  
  if (nrow(wilcox_results) > 0) {
    write.csv(wilcox_results,
              file = paste0("better_annotations_and_stats/boxplots/", gene, "_wilcox_results.csv"),
              row.names = FALSE)
  }
}
