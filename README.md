Males vs Females
Gene Expression Analysis and Visualization
This repository contains R scripts for analyzing RNA-seq count data, generating heatmaps for specific genes, and producing boxplots with statistical testing (Wilcoxon test) to compare gene expression between male and female samples across multiple time points.

📂 Project Structure
Input files required:

Counts data_for use.csv — Raw count matrix (genes × samples)

Metadata.csv — Sample metadata including SampleName, Week, and Gender

GeneSymbols.csv — Mapping of GeneID to GeneSymbol

Output:

Heatmaps in better_annotations_and_stats/

Boxplots and Wilcoxon test results in better_annotations_and_stats/boxplots/

📜 Analysis Workflow
Prepare Data

Load count matrix, metadata, and gene symbols.

Ensure sample ordering and proper factor levels for Week.

Filter lowly expressed genes (mean count < 10).

Generate Heatmaps

Focused on a curated list of ion channel and cytokine genes (e.g., TRP family, IL17 family, Ifng, Tnf).

For each week (Week2, Week4, Week8, Week16, Week30):

Scale gene expression data (Z-score normalization).

Generate non-clustered heatmaps comparing male vs female samples.

Statistical Testing and Boxplots

Normalize counts using DESeq2 size factor estimation.

For each gene:

Create gender-separated boxplots across weeks.

Perform Wilcoxon rank-sum tests per week.

Adjust p-values using Benjamini-Hochberg correction.

Annotate significant differences on boxplots.

📦 Required R Packages
Make sure you have the following R packages installed:

r

install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "ggpubr", "FSA", "rstatix"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
🚀 How to Run
Place your input files (Counts data_for use.csv, Metadata.csv, GeneSymbols.csv) in the working directory.

Open the R script(s).

Run the code.

Find outputs in the better_annotations_and_stats/ directory.

🧬 Genes of Interest
Focus is on genes related to ion channels and inflammatory cytokines, including:

TRP Channels: Trpa1, Trpv1, Trpm8, etc.

Cytokines: Il17a, Il17f, Il4ra, Ifng, Tnf, etc.

(Full list embedded in the code.)

📊 Outputs Examples
Heatmaps: Specific_Genes_Male_vs_Female_Week2_heatmap.png, etc.

Boxplots: Trpa1_combined_boxplot.png, Il17a_combined_boxplot.png, etc.

Statistics: Trpa1_wilcox_results.csv, etc.

✍️ Notes
The code automatically skips genes with missing or zero expression.

p-values are corrected for multiple comparisons (FDR method).

The analysis focuses on expression differences between males and females at different time points.

Male and Females seperately

This repository contains R scripts for analyzing RNA-seq count data, generating heatmaps for specific genes, and producing boxplots with statistical testing (Wilcoxon test) to compare gene expression of male and female samples across multiple time points seperately.

📂 Project Structure

Input files required:

Counts data_for use.csv — Raw count matrix (genes × samples)

Metadata.csv — Sample metadata including SampleName, Week, and Gender

GeneSymbols.csv — Mapping of GeneID to GeneSymbol

Output:

Heatmaps in better_annotations_and_stats/

Boxplots and Wilcoxon test results in better_annotations_and_stats/boxplots/

📜 Analysis Workflow

Prepare Data
Load count matrix, metadata, and gene symbols.

Ensure correct sample ordering and factor levels for Week.

Filter out lowly expressed genes (mean count < 10).

Generate Heatmaps
Focus on a curated list of ion channel and cytokine genes (e.g., TRP family, IL17 family, Ifng, Tnf).

For each week (Week2, Week4, Week8, Week16, Week30):

Scale gene expression data using Z-score normalization.

Generate non-clustered heatmaps comparing male vs female samples.

Statistical Testing and Boxplots
Normalize counts using DESeq2 size factor estimation.

For each gene:

Create gender-separated boxplots across weeks.

Perform Wilcoxon rank-sum test per week.

Adjust p-values using Benjamini-Hochberg correction.

Annotate significant differences on boxplots.

📦 Required R Packages

Make sure you have the following R packages installed:


install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "ggpubr", "FSA", "rstatix"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
🚀 How to Run

Place your input files (Counts data_for use.csv, Metadata.csv, GeneSymbols.csv) in the working directory.

Open the R script(s).

Run the code

Find outputs in the better_annotations_and_stats/ directory.

🧬 Genes of Interest

Focus is on genes related to ion channels and inflammatory cytokines, including:

TRP Channels: Trpa1, Trpv1, Trpm8, etc.

Cytokines: Il17a, Il17f, Il4ra, Ifng, Tnf, etc.

(Full list embedded in the code.)

📊 Outputs Examples

Heatmaps: Specific_Genes_Male_vs_Female_Week2_heatmap.png, etc.

Boxplots: Trpa1_combined_boxplot.png, Il17a_combined_boxplot.png, etc.

Statistics: Trpa1_wilcox_results.csv, etc.

✍️ Notes

Genes with missing or zero expression are automatically skipped.

p-values are corrected for multiple comparisons using the FDR method.

The analysis focuses on expression differences between males and females across multiple time points.
