# Heatmap Theme Functions for DGAT Immunology Analysis
# Provides consistent styling for heatmap visualizations

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Custom theme for heatmaps
theme_heatmap <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

# Color palette for heatmaps
heatmap_colors <- function(n = 100) {
  colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(n)
}

# Correlation heatmap colors
correlation_colors <- function(n = 100) {
  colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                     "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", 
                     "#D6604D", "#B2182B", "#67001F"))(n)
}

# Survival plot theme
theme_survival <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

# Forest plot theme
theme_forest <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90"),
      axis.line = element_line(color = "black", size = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Volcano plot theme
theme_volcano <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

# Box plot theme
theme_boxplot <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      strip.text = element_text(size = 11, face = "bold")
    )
}

# Scatter plot theme
theme_scatter <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

# Bar plot theme
theme_barplot <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

# Function to create publication-ready heatmap
create_publication_heatmap <- function(data, title = "Heatmap", 
                                      color_scheme = "RdBu",
                                      cluster_rows = TRUE, 
                                      cluster_cols = TRUE,
                                      show_rownames = TRUE,
                                      show_colnames = TRUE,
                                      fontsize = 10,
                                      width = 12,
                                      height = 8) {
  
  # Set up color palette
  if (color_scheme == "RdBu") {
    colors <- heatmap_colors()
  } else if (color_scheme == "correlation") {
    colors <- correlation_colors()
  } else {
    colors <- colorRampPalette(rev(brewer.pal(n = 11, name = color_scheme)))(100)
  }
  
  # Create heatmap
  pheatmap(data,
           color = colors,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           show_rownames = show_rownames,
           show_colnames = show_colnames,
           fontsize = fontsize,
           main = title,
           width = width,
           height = height,
           border_color = NA,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation")
}

# Function to save plots with consistent settings
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
}

# Export theme functions
cat("Heatmap theme functions loaded successfully!\n")
cat("Available functions:\n")
cat("- theme_heatmap()\n")
cat("- theme_survival()\n")
cat("- theme_forest()\n")
cat("- theme_volcano()\n")
cat("- theme_boxplot()\n")
cat("- theme_scatter()\n")
cat("- theme_barplot()\n")
cat("- create_publication_heatmap()\n")
cat("- save_plot()\n")
