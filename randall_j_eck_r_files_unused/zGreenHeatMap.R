# Load necessary libraries
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr) 
library(heatmaply)

df <- read_xlsx("zGreenExperiment.xlsx", sheet = "NormalizedCountsforHeatMap")
df <- as.data.frame(df)
rownames(df) <- df$Gene
df <- df[,-1]
df_normalized <- t(apply(df, 1, function(x) x / mean(x)))
df_nom <- data.frame(rownames(df_normalized), df_normalized)
df_nom <- na.omit(df_nom)
# Define the list of genes you're interested in
genes_of_interest <- c(
  "Ano7",
  "Spp1",
  "Zfp959",
  "Stn1",
  "Ccl5",
  "Cryab",
  "Cxcl10",
  "Prpf31",
  "Ifi44l",
  "Stmn1",
  "Xaf1",
  "Cd72",
  "B3galt6",
  "Csrp1",
  "Zmat1",
  "Dipk2a",
  "Irf7",
  "2610528A11Rik",
  "Ifi202b",
  "Wtip",
  "Mrm3",
  "Pwp1",
  "Ttc12",
  "Myo1g")
# Subset the data to include only the genes of interest
data_subset <- df_nom[df_nom$rownames.df_normalized. %in% genes_of_interest, ]
melted_df <- melt(data_subset)
colnames(melted_df) <- c("Gene", "Sample", "Expression")


# Create the heatmap
heatmap_plot <- ggplot(melted_df, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile(color = "white", alpha = 1) +
  scale_fill_gradient(low = "#00C2F9", high = "#E20134") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot area
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = .5),  # Add black axis lines
    axis.title = element_text(color = "black", size = 10),  # Customize axis title text color and size
    axis.text = element_text(color = "black", size = 8),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black")# Customize axis text color and size
  ) 

# Print the heatmap
print(heatmap_plot)
output_file <- "zGreenHeatMap.png"  # Update with the desired output file path
ggsave(output_file, plot = heatmap_plot, width = 4, height = 4, dpi = 600)

