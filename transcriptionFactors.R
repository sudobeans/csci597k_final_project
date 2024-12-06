library(readxl)
library(ggplot2)
# Note that dplyr masks some classic R functions like filter,
# lag, intersect, setdiff, setequal, union.
library(dplyr, warn.conflicts=FALSE)
library(ggrepel)
library(readr)

usage = "transcriptionFactors.R [chea3_output_tsv_file_name] [graph_output_file_name]"

args <- commandArgs(trailingOnly = TRUE)
cat(args)
if (length(args) != 2) {
  print(paste("Arguments: ", args))
  stop(paste("Not enough arguments! USAGE: ", usage))
}

input_file <- args[1]
output_file <- args[2]
print(paste("Reading from", input_file, "printing to", output_file))

df <- read_tsv(input_file)

graph_from_df <- function(df) {
  df$tf <- as.character(df$`TF`)
  df$log <- as.numeric(df$`Odds Ratio`)
  df$p <- as.numeric(df$`FET p-value`)
  df$rank <- as.numeric(df$`Rank`)
  
  # Graph top 50 transcription factors
  df <- df %>%
    mutate(
      significance = ifelse(rank <= 50 & (log > 0 | log < 0), "Significant", "Not Significant"),
      color = ifelse(p < .05 & log > 0, "#E20134", 
                     ifelse(p < .05 & log < -0, "#00C2F9", "black"))
    )
  
  
  # Select top 5 ranked TFs to highlight
  tfs = df$tf
  genes_to_highlight <- df$tf[df$rank <= 5]
  
  
  # Filter significant points based on gene names
  df_labeled <- df %>%
    filter(`TF`%in% genes_to_highlight)
  
  
  volcano_plot <- ggplot(df, aes(x = log, y = -log10(p))) +
    geom_point(aes(color = color), size = .8, alpha = 1) +   # Set points to one color and adjust size
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",linewidth =.5 ) +  # p-value threshold line
    #geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black",linewidth =.5 ) +   # fold change threshold lines
    geom_text_repel(data = df_labeled, aes(label = `TF`,color = color), size = 2,nudge_x = -0.03, nudge_y = 0.03) + 
    scale_color_identity() + # Label specified genes
    labs(x = "Odds Ratio",
         y = "FET p-value") +
    scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +  # Set x-axis limits and label frequency
    scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",  # Remove the legend
      panel.background = element_rect(fill = "transparent", color = NA),  # Transparent background
      plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot area
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = .5 ),  # Add black axis lines
      axis.title = element_text(color = "black", size = 10),  # Customize axis title text color and size
      axis.text = element_text(color = "black", size =8)    # Customize axis text color and size
    )
  
  return(volcano_plot)
}

volcano_plot = graph_from_df(df)
ggsave(output_file, plot = volcano_plot, width = 2, height = 2, dpi = 600)
