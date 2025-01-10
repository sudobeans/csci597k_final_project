library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

df <- read_xlsx("zGreenExperiment.xlsx", sheet = "TranscriptionFactorNetwork")

df$log <- as.numeric(df$`Odds Ratio`)
df$p <- as.numeric(df$padj)

df <- df %>%
  mutate(
    significance = ifelse(p < .05 & (log > 0 | log < 0), "Significant", "Not Significant"),
    color = ifelse(p < .05 & log > 0, "#E20134", 
                   ifelse(p < .05 & log < -0, "#00C2F9", "black"))
  )


# List of gene names to be highlighted
genes_to_highlight <- c("E2F4",
                        "STAT2")


# Filter significant points based on gene names
df_labeled <- df %>%
  filter(`Transcription Factor`%in% genes_to_highlight)


volcano_plot <- ggplot(df, aes(x = log, y = -log10(p))) +
  geom_point(aes(color = color), size = .8, alpha = 1) +   # Set points to one color and adjust size
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",linewidth =.5 ) +  # p-value threshold line
  #geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black",linewidth =.5 ) +   # fold change threshold lines
  geom_text_repel(data = df_labeled, aes(label = `Transcription Factor`,color = color), size = 2,nudge_x = -0.03, nudge_y = 0.03) + 
  scale_color_identity() + # Label specified genes
  labs(x = "Odds Ratio",
       y = "-Log10(adjusted p value)") +
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

# Print the volcano plot
#print(volcano_plot)

ggsave("zGreenTranscriptionFactors.png", plot = volcano_plot, width = 2, height = 2, dpi = 600)
