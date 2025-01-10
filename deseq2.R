# Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2 from Bioconductor
BiocManager::install("DESeq2")

# Load DESeq2 library
library(DESeq2)

# used for count matrix
library(dplyr)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# Function to map Ensembl gene IDs to gene symbols
add_gene_symbols <- function(data) {
  gene_ids <- rownames(data)
  annots <- select(org.Hs.eg.db,
                   keys = gene_ids,
                   columns = "SYMBOL",
                   keytype = "ENSEMBL")
  data <- merge(data,
                annots,
                by.x = "row.names",
                by.y = "ENSEMBL",
                all.x = TRUE)
  rownames(data) <- data$Row.names
  data <- data[, -1]  # Remove redundant Row.names column
  return(data)
}

run_individual_tests <- function(count_matrix,
                                 treatments,
                                 colData,
                                 control_label_0M,
                                 control_label_0P) {
  # Ensure condition is a factor
  colData$condition <- as.factor(colData$condition)
  
  # Initialize DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = colData,
                                design = ~ condition)
  
  # filter out rows with low gene counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Store results for each treatment
  results_list <- list()
  
  # Filter function without additional libraries
  filter_results <- function(res_df,
                             lfc_thresh = 1,
                             padj_thresh = 0.05) {
    # Remove rows with NA values in log2FoldChange or padj
    res_df <- res_df[!is.na(res_df$log2FoldChange) &
                       !is.na(res_df$padj), ]
    # filter by threshold values
    res_df <- res_df[abs(res_df$log2FoldChange) > lfc_thresh &
                       res_df$padj < padj_thresh, ]
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

# Define the function
load_count_matrix <- function() {
  # File paths for the counts files
  counts_files <- list(
    "0M" = "htseqModResults/0M_output.csv",
    "0P" = "htseqModResults/0P_output.csv",
    "1A" = "htseqModResults/1a_output.csv",
    "2P" = "htseqModResults/2P_output.csv",
    "3L" = "htseqModResults/3L_output.csv",
    "4C" = "htseqModResults/4C_output.csv"
  )
  
  # Read all files into a list of data frames
  dataframes <- lapply(names(counts_files), function(sample) {
    file_path <- counts_files[[sample]]
    df <- read.csv(
      file_path,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(df) <- c("geneID", paste0(sample, "_col", seq(1, ncol(df) - 1)))
    return(df)
  })
  
  # Merge all data frames on "geneID" using Reduce and full_join
  count_matrix <- Reduce(function(x, y) {
    full_join(x, y, by = "geneID")
  }, dataframes)
  
  # Replace NA values with zeros
  count_matrix[is.na(count_matrix)] <- 0
  
  # Filter out rows (genes) with zero counts across all samples
  count_matrix <- count_matrix[rowSums(count_matrix[, -1]) > 0, ]
  
  # Filter gene IDs that start with 'E'
  count_matrix <- count_matrix[grepl("^E", count_matrix$geneID), ]
  
  # Set geneID as rownames and remove the geneID column
  rownames(count_matrix) <- count_matrix$geneID
  count_matrix <- count_matrix[, -1]  # Drop the geneID column
  
  return(count_matrix)
}

# obtain count matrix
count_matrix <- load_count_matrix()

count_matrix <- round(count_matrix)  # Ensure counts are integers

# Define metadata
colData <- data.frame(
  condition = c(
    rep("0M", 8),
    rep("0P", 8),
    rep("1A", 6),
    rep("2P", 8),
    rep("3L", 8),
    rep("4C", 8)
  ),
  row.names = colnames(count_matrix)
)

all(colnames(count_matrix) %in% rownames (colData))

all(colnames(count_matrix) == rownames(colData))

# Define treatments
treatments <- c("1A", "2P", "3L", "4C")

output_dir <- "~/Downloads/filteredResults"

# run the tests
results <- run_individual_tests(
  count_matrix = count_matrix,
  treatments = treatments,
  colData = colData,
  control_label_0M = "0M",
  control_label_0P = "0P"
)

# Inspect results and print to files
for (name in names(results)) {
  cat("\nComparison:", name, "\n")
  print(results[[name]])
  # Write filtered results to CSV in the downloads directory
  write.csv(results[[name]],
            file = file.path(output_dir
                             , paste0(name, "_filtered.csv")),
            row.names = TRUE)
}
