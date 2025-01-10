# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)
library(readxl)

# Read the input spreadsheet
go_data <- read_xlsx("mCherryExperiment.xlsx", sheet = "GeneSetEnrichmentAnalysis")

# Check the data
head(go_data)

# Filter and arrange the data for better visualization
# Select the top 10 terms based on p_value
top_go_data <- go_data %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>%
  head(10)

# Rename the top 10 GO terms
new_names <- c("Term 1", "Term 2", "Term 3", "Term 4", "Term 5", 
               "Term 6", "Term 7", "Term 8", "Term 9", "Term 10")
top_go_data$GO_term <- new_names

# Create a new column for ordering in descending order of p value
top_go_data <- top_go_data %>%
  mutate(ordered_go_term = factor(`Simplified Name`, levels = `Simplified Name`[order(-pvalue)]))

# Plot
dot_plot <- ggplot(top_go_data, aes(x = -log10(pvalue), y = ordered_go_term, size = Intercept, color = `Odds Ratio`)) +
  geom_point() +
  labs(x = "-log10(p-value)",
       y = "Microglia Gene Set",
       color = "Enrichment",
       size = "Count") +
  scale_color_gradient(low = "#00C2F9", high = "#E20134")+
  scale_x_continuous(limits = c(0, max(-log10(top_go_data$pvalue)) + 1)) +
  theme_minimal() +
  theme(
    text = element_text(color = "black"),  # Set all text to black
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.height    = unit(0.2, "cm"), 
    panel.background = element_rect(fill = "transparent", color = NA),  # Make background transparent
    plot.background = element_rect(fill = "transparent", color = NA),  # Make background transparent
    panel.grid.major = element_line(color = "black"),  # Major gridlines black
    panel.grid.minor = element_line(color = "black")  # Minor gridlines black
  ) 



# Export the plot to a PNG file
output_file <- "mCherryGSEA.png"  # Update with the desired output file path
ggsave(output_file, plot = dot_plot, width = 4, height = 2, dpi = 600)