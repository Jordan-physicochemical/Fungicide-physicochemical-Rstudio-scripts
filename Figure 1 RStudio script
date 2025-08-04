Figure 1 (physicochemical property spread)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Load your data - change file path accordingly
fungicide_data <- read.csv("20250624 Physicochemical properties for Figure 1 selection of fungicides.csv")

properties <- c(
  "MolLogP",
  "LogS.Chemicalize",
  "MolWt",
  "MolMR",
  "Maximum.projection.radius.Chemicalize",
  "TPSA",
  "Polarizability.Chemicalize",
  "BalabanJ"
)

if(!"Highlight" %in% colnames(fungicide_data)){
  stop("Column 'Highlight' not found in fungicide_data")
}


highlighted_fungicides <- unique(fungicide_data$Fungicide[fungicide_data$Highlight == "Y"])
num_colors <- length(highlighted_fungicides)
if(num_colors <= 8){
  palette_colors <- RColorBrewer::brewer.pal(num_colors, "Set1")
} else {
  palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_colors)
}
names(palette_colors) <- highlighted_fungicides

plot_property <- function(property_name) {
  if(!property_name %in% colnames(fungicide_data)){
    stop(paste("Property", property_name, "not found in fungicide_data"))
  }
  
  # Convert to numeric safely
  fungicide_data[[property_name]] <- as.numeric(as.character(fungicide_data[[property_name]]))
  
  min_val <- min(fungicide_data[[property_name]], na.rm = TRUE)
  max_val <- max(fungicide_data[[property_name]], na.rm = TRUE)
  
  fungicide_data$y_plot <- fungicide_data[[property_name]]
  fungicide_data$y_plot[fungicide_data$Highlight == "Y"] <- fungicide_data$y_plot[fungicide_data$Highlight == "Y"] + offset
  
  p_left <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(property_name, "\nMin: ", round(min_val, 2)),
             size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(5,5,5,5))
  
  p_middle <- ggplot(fungicide_data, aes(y = "", x = .data[[property_name]])) +
    geom_violin(fill = "lightblue") +
    geom_jitter(data = subset(fungicide_data, Highlight == "N"),
                aes(x = y_plot),
                height = 0.1, size = 1.5, alpha = 0.6, color = "black") +
    geom_jitter(data = subset(fungicide_data, Highlight == "Y"),
                aes(x = y_plot, fill = Fungicide),
                height = 0.1, size = 3.75, alpha = 0.9,
                shape = 21, color = "black", stroke = 0.8) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),    # Remove x axis tick labels
          axis.ticks.x = element_blank(),   # Remove x axis ticks
          plot.margin = margin(5,5,5,5)) +
    scale_fill_manual(values = palette_colors, name = "Highlighted Fungicide")
  
  p_right <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("Max: ", round(max_val, 2)),
             size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(5,5,5,5))
  
  combined <- p_left + p_middle + p_right + plot_layout(widths = c(1,3,1))
  return(combined)
}


plot_list <- lapply(properties, plot_property)

final_plot_violin1 <- wrap_plots(plot_list, ncol = 1, heights = rep(1, length(plot_list)))

print(final_plot_violin1)

ggsave("Figure 1 fungicide selection violin plots only.svg", plot = final_plot_violin1, width = 12, height = 8, units = "in", dpi = 300)




### boxplots for noncontinous data
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Load your data
fungicide_data <- read.csv("20250624 Physicochemical properties for Figure 1 selection of fungicides.csv")

properties_box <- c(
  "NumHAcceptors",
  "NumHDonors",
  "NOCount",
  "NHOHCount",
  "NumRotatableBonds",
  "NumValenceElectrons",
  "fr_amide",
  "NumAromaticRings",
  "NumSaturatedRings",
  "NumAliphaticRings"
)

if(!"Highlight" %in% colnames(fungicide_data)){
  stop("Column 'Highlight' not found in fungicide_data")
}

highlighted_fungicides <- unique(fungicide_data$Fungicide[fungicide_data$Highlight == "Y"])
num_colors <- length(highlighted_fungicides)
if(num_colors <= 8){
  palette_colors <- RColorBrewer::brewer.pal(num_colors, "Set1")
} else {
  palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_colors)
}
names(palette_colors) <- highlighted_fungicides

# Function to create a 3-panel (left-text, boxplot, right-text) layout
plot_box_property <- function(property_name) {
  if(!property_name %in% colnames(fungicide_data)){
    stop(paste("Property", property_name, "not found in fungicide_data"))
  }
  
  # Ensure numeric
  fungicide_data[[property_name]] <- as.numeric(as.character(fungicide_data[[property_name]]))
  
  # Calculate min and max
  min_val <- min(fungicide_data[[property_name]], na.rm = TRUE)
  max_val <- max(fungicide_data[[property_name]], na.rm = TRUE)
  
  # Left facet: property name + min
  p_left <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(property_name, "\nMin: ", round(min_val, 2)),
             size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Middle facet: horizontal boxplot + jittered points
  p_middle <- ggplot(fungicide_data, aes(y = "", x = .data[[property_name]])) +
    geom_boxplot(outlier.shape = NA, fill = "lightblue") +
    geom_jitter(data = subset(fungicide_data, Highlight == "N"),
                aes(x = .data[[property_name]]),
                width = 0.1, height = 0.3, size = 1.5, alpha = 0.6, color = "black") +
    geom_jitter(data = subset(fungicide_data, Highlight == "Y"),
                aes(x = .data[[property_name]], fill = Fungicide),
                width = 0.1, height = 0.3, size = 3.75, alpha = 0.9,
                shape = 21, color = "black", stroke = 0.8) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # remove x-axis tick labels
      axis.ticks.x = element_blank(), # remove x-axis ticks
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    scale_fill_manual(values = palette_colors, name = "Highlighted Fungicide")
  
  
  # Right facet: max value
  p_right <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("Max: ", round(max_val, 2)),
             size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Combine left + middle + right horizontally
  combined <- p_left + p_middle + p_right + plot_layout(widths = c(1, 3, 1))
  return(combined)
}

# Create list of 3-panel plots
plot_list_box <- lapply(properties_box, plot_box_property)

# Stack all vertically
final_plot_box1 <- wrap_plots(plot_list_box, ncol = 1)

# Show plot
print(final_plot_box1)
ggsave("Figure 1 fungicide selection box plots only.svg", plot = final_plot_box1, width = 12, height = 8, units = "in", dpi = 300)




###combined plots
final_plot <- final_plot_box1 + final_plot_violin1 + plot_layout(ncol = 1, heights = c(1, 1))

final_plot


library (svglite)
ggsave("Figure 1 fungicide selection combined violin box plots.svg", plot = final_plot, width = 12, height = 16, units = "in", dpi = 300)


### save legend
library(cowplot)
library(grid)

# Extract legend from the first boxplot (adjust index if needed)
legend_only <- get_legend(plot_list_box[[1]])

# Draw legend on a new page
grid.newpage()
grid.draw(legend_only)

library(ggplot2)

# Convert to ggdraw for saving
legend_plot <- cowplot::ggdraw(legend_only)

# Save as SVG
ggsave("boxplot_legend.svg", plot = legend_plot, width = 4, height = 2, device = "svg")

