# Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2 from Bioconductor
BiocManager::install("DESeq2")

# Load DESeq2 library
library(DESeq2)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# Function to map Ensembl gene IDs to gene symbols
add_gene_symbols <- function(data) {
  gene_ids <- rownames(data)
  annots <- select(org.Hs.eg.db, keys = gene_ids, columns = "SYMBOL", keytype = "ENSEMBL")
  data <- merge(data, annots, by.x = "row.names", by.y = "ENSEMBL", all.x = TRUE)
  rownames(data) <- data$Row.names
  data <- data[, -1]  # Remove redundant Row.names column
  return(data)
}

run_individual_tests <- function(count_matrix, treatments, colData, control_label_0M, control_label_0P) {
  # Ensure condition is a factor
  colData$condition <- as.factor(colData$condition)
  
  # Initialize DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = colData,
                                design = ~ condition)
  
  # filter out rows with low gene counts
  dds <- dds[rowSums(counts(dds)) > 1,]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Store results for each treatment
  results_list <- list()
  
  # Filter function without additional libraries
  filter_results <- function(res_df, lfc_thresh = 1, padj_thresh = 0.05) {
    # Remove rows with NA values in log2FoldChange or padj
    res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$padj), ]
    # filter by threshold values
    res_df <- res_df[abs(res_df$log2FoldChange) > lfc_thresh & res_df$padj < padj_thresh, ]
    return(res_df)
  }
  
  # Run comparisons for treatments vs. 0M control
  for (treatment in treatments) {
    res <- results(dds, contrast = c("condition", treatment, control_label_0M))
    res_df <- as.data.frame(res)  # Convert to data frame
    filtered_res <- filter_results(res_df)  # Filter results
    # Add gene symbols
    filtered_res <- add_gene_symbols(filtered_res)
    results_list[[paste0(treatment, "_vs_", control_label_0M)]] <- filtered_res
  }
  
  # Run comparisons for treatments vs. 0P control
  for (treatment in treatments) {
    res <- results(dds, contrast = c("condition", treatment, control_label_0P))
    res_df <- as.data.frame(res)  # Convert to data frame
    filtered_res <- filter_results(res_df)  # Filter results
    # Add gene symbols
    filtered_res <- add_gene_symbols(filtered_res)
    results_list[[paste0(treatment, "_vs_", control_label_0P)]] <- filtered_res
  }
  
  return(results_list)
}

# Load and transpose the count matrix
count_matrix <- read.csv('countMatrix.csv',
                         header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
count_matrix <- round(count_matrix)  # Ensure counts are integers

# Define metadata
colData <- data.frame(
  condition = c(rep("0M", 8), rep("0P", 8), rep("1A", 6), rep("2P", 8), rep("3L", 8), rep("4C", 8)),
  row.names = colnames(count_matrix)
)

all(colnames(count_matrix) %in% rownames (colData))

all(colnames(count_matrix) == rownames(colData))

# Define treatments
treatments <- c("1A", "2P", "3L", "4C")

output_dir <- "deseq2_filtered_results"

# run the tests
results <- run_individual_tests(
  count_matrix = count_matrix,
  treatments = treatments,
  colData = colData,
  control_label_0M = "0M",
  control_label_0P = "0P"
)

# Inspect results
for (name in names(results)) {
  cat("\nComparison:", name, "\n")
  print(results[[name]])
  # Write filtered results to CSV in the downloads directory
  write.csv(results[[name]], file = file.path(output_dir, paste0(name, "_filtered.csv")), row.names = TRUE)
}
