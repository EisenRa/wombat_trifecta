##### Function for plotting dRep MASH outputs into a tree
##### Raphael Eisenhofer 3/2025

library(tidyverse)
library(ggdendro)
library(scales)

plot_mash_tree <- function(Mdb, Cdb, threshold = FALSE, plot_dir = NULL) {
  # Set ggplot theme
  theme_set(theme_minimal() +
              theme(panel.grid = element_blank(),
                    axis.line.x = element_line(),
                    axis.ticks.y = element_blank(),
                    plot.margin = margin(l = 100, r = 20)))
  
  # Convert to distance matrix
  dist_mat <- Mdb %>%
    select(genome1, genome2, dist) %>%
    pivot_wider(names_from = genome2, values_from = dist) %>%
    column_to_rownames("genome1") %>%
    as.matrix()
  
  # Make sure it's a proper distance matrix (symmetric with 0 diagonal)
  dist_mat <- as.dist(dist_mat)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_mat, method = "average")
  
  # Prepare cluster assignments
  name2cluster <- Cdb %>% 
    select(genome, primary_cluster) %>% 
    deframe()
  
  # Generate colors
  name2color <- gen_color_dictionary(hc$labels, name2cluster)
  
  # Convert hclust to dendrogram
  dend <- as.dendrogram(hc)
  
  # Prepare dendrogram data for plotting
  dend_data <- dendro_data(dend)
  
  # Prepare segment data
  segment_data <- segment(dend_data) %>%
    mutate(yend = ifelse(yend < 0, 0, yend))  # Ensure no negative values
  
  # Prepare label data with colors and clusters
  label_data <- label(dend_data) %>%
    mutate(color = unname(name2color),
           # Add secondary cluster info
           label_with_cluster = paste0(label, " (", 
                                       Cdb$secondary_cluster[match(label, Cdb$genome)], ")"))
  
  # Create the dendrogram plot
  p <- ggplot() +
    geom_segment(data = segment_data, 
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label_data,
              aes(x = x, y = -0.01, label = label_with_cluster, color = color),
              hjust = 0, size = 3) +
    scale_y_reverse(expand = expansion(mult = c(0.1, 0.1)),
                    labels = function(x) sprintf("%.2f", (1 - x) * 100)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    scale_color_identity() +
    labs(title = "MASH clustering",
         x = "MASH Average Nucleotide Identity (ANI)",
         y = "") +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10, margin = margin(t = 10)))
  
  # Add threshold line if specified
  if (!isFALSE(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed")
  }
  
  # Adjust plot size based on number of labels
  n_labels <- length(hc$labels)
  plot_height <- max(6, n_labels * 0.3)  # Minimum height of 6 inches
  
  # Save or display the plot
  if (!is.null(plot_dir)) {
    ggsave(file.path(plot_dir, "Primary_clustering_dendrogram.pdf"), 
           p, width = 10, height = plot_height, 
           device = "pdf", bg = "transparent")
  }
  
  print(p)
  return(invisible(p))
}

# Helper function to generate colors
gen_color_dictionary <- function(names, name2cluster) {
  clusters <- unique(name2cluster[names])
  n_clusters <- length(clusters)
  
  # Use a color palette that works well for clusters
  if (n_clusters <= 8) {
    pal <- RColorBrewer::brewer.pal(max(3, n_clusters), "Set2")
  } else if (n_clusters <= 12) {
    pal <- RColorBrewer::brewer.pal(n_clusters, "Set3")
  } else {
    pal <- scales::hue_pal()(n_clusters)
  }
  
  colors <- setNames(pal[1:n_clusters], clusters)
  return(colors[as.character(name2cluster[names])])
}