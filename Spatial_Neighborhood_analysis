# --- 1. User-Defined Parameters ---
# IMPORTANT: Adjust these parameters to match your specific dataset and analysis goals.

# 1.1. Path to your Seurat object RDS file.
# Example: seurat_object_path <- "/path/to/your/combined_seurat_object.rds"
seurat_object_path <- "ENTER_YOUR_SEURAT_OBJECT_PATH_HERE"

# 1.2. Name of the assay to use in your Seurat object (e.g., "Xenium", "RNA", "SCT").
# This specifies which data layer to use for spatial coordinates or cell identities.
assay_name <- "Xenium" # Default value, change if needed

# 1.3. Name of the target cell type for microenvironment analysis.
# This should exactly match a cell type label in your 'annotated_cell_type_column'.
# Example: target_cell_type_name <- "T/NK cell"
target_cell_type_name <- "ENTER_YOUR_TARGET_CELL_TYPE_NAME_HERE"

# 1.4. Name of the metadata column in your Seurat object containing annotated cell types.
# This column will be used to identify your 'target_cell_type_name' and other cell types.
# Example: annotated_cell_type_column <- "SingleR_full"
annotated_cell_type_column <- "ENTER_YOUR_CELL_TYPE_METADATA_COLUMN_HERE"

# 1.5. Maximum distance for defining neighbors (in spatial units, e.g., micrometers).
# Cells within this radius will be considered neighbors for your target cells.
# Adjust this value based on your spatial data's resolution and the biological scale of interactions.
neighbor_distance_um <- 80 # Default value, adjust based on your data

# 1.6. Number of clusters for K-means and NMF clustering.
# This determines how many distinct microenvironment clusters the target cells will be grouped into.
# You may need to experiment with different values to find optimal biological groupings.
num_clusters <- 5 # Default value, adjust as needed

# 1.7. Name of the metadata column containing sample or response information.
# This column will be used for visualizing cluster distributions across different groups (e.g., "Response", "patient_group", "treatment").
# Ensure this column exists in your Seurat object's metadata.
response_column_name <- "Response" # Default value, adjust if your column name differs

# --- 2. Load Required Packages ---
# These packages are essential for the script's functionality.
# If any are not installed, uncomment and run install.packages("package_name") commands in your R console.
# install.packages("Seurat")
# install.packages("dplyr")
# install.packages("SpatialExperiment")
# install.packages("SeuratObject")
# install.packages("sp")
# install.packages("spdep")
# install.packages("STutility")
# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("NMF")
# install.packages("RColorBrewer")
# install.packages("ggpubr")

library(Seurat)          # For single-cell and spatial omics data analysis
library(dplyr)           # For data manipulation and piping operations (%>%)
library(SpatialExperiment) # For spatial transcriptomics data structures (often used with spatial data)
library(SeuratObject)    # Core Seurat S4 object methods
library(sp)              # For spatial data classes and methods
library(spdep)           # For spatial dependence modeling, including neighbor finding (dnearneigh)
library(STutility)       # Utility functions for spatial transcriptomics analysis
library(tidyr)           # For data tidying, including reshaping (pivot_longer)
library(ggplot2)         # For creating high-quality data visualizations
library(NMF)             # For Non-negative Matrix Factorization
library(RColorBrewer)    # For color palettes in visualizations
library(ggpubr)          # For publication-ready plots based on ggplot2

# --- 3. Load Data and Initial Setup ---
# Load the Seurat object from the specified path provided in User-Defined Parameters.
Seurat_object <- readRDS(seurat_object_path) # Renamed to avoid conflict with Seurat package name

# Set options to increase the maximum allowed size for global variables in 'future' plans.
# This helps prevent memory issues, especially with large spatial datasets.
options(future.globals.maxSize = 80000 * 1024^2) # Set to 80GB (adjust if you have less RAM)

# Set the default assay using the user-defined parameter.
DefaultAssay(Seurat_object) <- assay_name

# --- 4. Add Spatial Coordinates to Metadata ---
# Extract spatial coordinates (centroids) from each image within the Seurat object.
# These coordinates represent the center of each cell's boundary and are crucial for spatial analysis.
message("Extracting spatial coordinates...")
all_coords <- lapply(names(Seurat_object@images), function(img_name) {
  coords <- Seurat_object@images[[img_name]]@boundaries$centroids@coords
  coords_df <- as.data.frame(coords)
  # Assign cell barcodes as row names for easy lookup and integration.
  rownames(coords_df) <- Seurat_object@images[[img_name]]$centroids@cells
  coords_df
})
# Combine coordinates from all images into a single data frame.
coords_df <- bind_rows(all_coords)

# Add the combined coordinate data frame as metadata to the Seurat object.
Seurat_object <- AddMetaData(Seurat_object, metadata = coords_df)
message("Spatial coordinates added to metadata.")

# --- 5. Define Cell Neighbors (Spatial Proximity Analysis) ---
# Identify neighboring cells for each cell within its respective image.
# The 'dnearneigh' function from 'spdep' is used to find neighbors within the specified distance range.
message("Defining cell neighbors...")
neighbor_results <- list()
for (img_name in names(Seurat_object@images)) {
  message(paste0("Processing image: ", img_name))
  coords <- Seurat_object@images[[img_name]]@boundaries$centroids@coords
  coords_df <- as.data.frame(coords)
  rownames(coords_df) <- Seurat_object@images[[img_name]]$centroids@cells
  coords_sp <- sp::SpatialPoints(coords_df) # Convert to SpatialPoints object for spatial analysis
  
  # Find neighbors within the user-defined distance (neighbor_distance_um).
  # 0 is the lower bound (no self-loops).
  neighbor_list <- dnearneigh(coords_sp, 0, neighbor_distance_um)
  neighbor_results[[img_name]] <- neighbor_list
}
message("Neighbor definitions complete.")

# --- 6. Calculate Neighboring Cell Type Proportions for Target Cells ---
# Focus on the user-defined target cell type and quantify the proportions of different cell types
# in their immediate neighborhood. This forms the basis for microenvironment characterization.
message(paste0("Calculating neighboring cell type proportions for '", target_cell_type_name, "' cells..."))
target_cells <- rownames(Seurat_object@meta.data)[Seurat_object@meta.data[[annotated_cell_type_column]] == target_cell_type_name]

# Get all unique cell types from the annotated metadata column.
all_cell_types <- unique(Seurat_object@meta.data[[annotated_cell_type_column]])

# Initialize a matrix to store cell type proportions.
# Rows correspond to target cells, columns to all unique cell types.
prop_matrix <- matrix(0, nrow = length(target_cells), ncol = length(all_cell_types))
rownames(prop_matrix) <- target_cells
colnames(prop_matrix) <- all_cell_types

# Iterate through each image and each target cell to calculate and store neighbor proportions.
for (img_name in names(neighbor_results)) {
  neighbor_list <- neighbor_results[[img_name]]
  cells <- Seurat_object@images[[img_name]]$centroids@cells
  target_idx <- which(cells %in% target_cells) # Get indices of target cells within the current image
  
  for (i in target_idx) {
    neighbors <- neighbor_list[[i]] # Get indices of neighbors for the current target cell
    if (length(neighbors) > 0) { # Proceed only if the target cell has neighbors
      barcode <- cells[i] # Barcode (unique identifier) of the current target cell
      neighbor_barcodes <- cells[neighbors] # Barcodes of its neighboring cells
      
      # Retrieve cell type labels for these neighboring cells using the user-defined column.
      neighbor_labels <- Seurat_object@meta.data[neighbor_barcodes, annotated_cell_type_column]
      
      # Calculate proportions of each cell type among the neighbors (relative frequency).
      tab <- prop.table(table(neighbor_labels))
      # Store these proportions in the main proportion matrix.
      prop_matrix[barcode, names(tab)] <- as.numeric(tab)
    }
  }
}

# Filter out target cells that had no neighbors (i.e., their row sum in prop_matrix is 0).
# These cells cannot be characterized by their neighborhood composition.
prop_matrix_filtered <- prop_matrix[rowSums(prop_matrix) > 0, ]
message("Neighboring cell type proportion calculation complete.")

# --- 7. Clustering: K-means and Non-negative Matrix Factorization (NMF) ---
# Cluster the target cells based on their neighboring cell type proportions.
# This step groups target cells with similar microenvironments.
message("Performing K-means and NMF clustering...")

# K-means clustering: a common method for partitioning data into k clusters.
set.seed(42) # Set a seed for reproducibility of random processes in clustering
km_res <- kmeans(prop_matrix_filtered, centers = num_clusters) # Use user-defined number of clusters

# Non-negative Matrix Factorization (NMF) clustering:
# NMF is particularly useful for parts-based representations (like cell type proportions).
set.seed(42) # For reproducibility
nmf_result <- nmf(prop_matrix_filtered, rank = num_clusters, method = "brunet", nrun = 10) # 'brunet' is a common NMF algorithm
nmf_basis <- basis(nmf_result) # Extract the basis matrix, which represents cluster profiles
# Assign each cell to the NMF cluster where it has the highest loading.
nmf_clusters <- apply(nmf_basis, 1, which.max)
message("Clustering complete.")

# --- 8. Integrate Clustering Results into Seurat Object Metadata and Save ---
# Add the K-means and NMF clustering results as new metadata columns to the Seurat object.
# These new columns will store the microenvironment cluster assignments for each target cell.
Seurat_object$cluster_kmeans_env <- NA
Seurat_object$cluster_kmeans_env[rownames(prop_matrix_filtered)] <- as.character(km_res$cluster)

Seurat_object$cluster_nmf_env <- NA
Seurat_object$cluster_nmf_env[rownames(prop_matrix_filtered)] <- as.character(nmf_clusters)

# Define output file names based on the target cell type for clear identification of results.
clean_target_name <- gsub(" ", "_", tolower(target_cell_type_name)) # Clean name for file paths
output_meta_file <- paste0(clean_target_name, "_cluster_meta_results.rds")
output_matrix_file <- paste0(clean_target_name, "_prop_matrix_filtered.rds")

# Save the updated Seurat object metadata. This includes the new cluster assignments.
saveRDS(Seurat_object@meta.data, output_meta_file)
message(paste0("Updated metadata saved to: ", output_meta_file))

# Convert the filtered proportion matrix to a data frame and save it.
# This matrix can be used for further downstream analysis or visualization.
subset_matrix_df <- as.data.frame(prop_matrix_filtered)
saveRDS(subset_matrix_df, output_matrix_file)
message(paste0("Filtered proportion matrix saved to: ", output_matrix_file))

# --- 9. Analyze and Visualize Cluster-specific Neighboring Cell Type Proportions ---
# Calculate and visualize the mean proportions of neighboring cell types for each cluster.
# This helps understand what defines each microenvironment cluster.
message("Generating cluster composition plots...")

# Add cluster assignments to the subset_matrix_df for plotting.
subset_matrix_df$cluster_kmeans <- Seurat_object$cluster_kmeans_env[rownames(subset_matrix_df)]
subset_matrix_df$cluster_nmf <- Seurat_object$cluster_nmf_env[rownames(subset_matrix_df)]

# Calculate mean cell type proportions for each K-means cluster.
cluster_means_kmeans <- subset_matrix_df %>%
  group_by(cluster_kmeans) %>%
  summarise(across(where(is.numeric), mean))

# Calculate mean cell type proportions for each NMF cluster.
cluster_means_nmf <- subset_matrix_df %>%
  group_by(cluster_nmf) %>%
  summarise(across(where(is.numeric), mean))

# Reshape data to 'long' format for ggplot2 visualization (suitable for bar plots/heatmaps).
cluster_means_kmeans_long <- pivot_longer(cluster_means_kmeans, -cluster_kmeans,
                                          names_to = "cell_type", values_to = "mean_prop")
cluster_means_nmf_long <- pivot_longer(cluster_means_nmf, -cluster_nmf,
                                       names_to = "cell_type", values_to = "mean_prop")

### Stacked Bar Plot of NMF Cluster Neighboring Cell Type Proportions
# Visualizes the average composition of neighboring cell types for each NMF-defined microenvironment cluster.
plot_stacked_bar <- ggplot(cluster_means_nmf_long, aes(x = factor(cluster_nmf), y = mean_prop, fill = cell_type)) +
  geom_col(position = "stack") + # Creates a stacked bar plot
  labs(x = paste0("NMF Cluster (", target_cell_type_name, ")"),
       y = "Mean Cell Type Proportion",
       fill = "Cell Type",
       title = paste0("NMF Cluster Composition for ", target_cell_type_name, " Microenvironment")) + # Clear labels
  theme_minimal(base_size = 14) + # Clean, minimal theme
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # Ensure x-axis labels are readable
ggsave(paste0(clean_target_name, "_nmf_cluster_stacked_barplot.pdf"), plot = plot_stacked_bar, width = 10, height = 7)

### Heatmap of NMF Cluster Neighboring Cell Type Proportions
# Provides a concise visual summary of which cell types are enriched or depleted in each NMF cluster.
plot_heatmap <- ggplot(cluster_means_nmf_long, aes(x = factor(cluster_nmf), y = cell_type, fill = mean_prop)) +
  geom_tile(color = "white") + # Creates a tile plot (heatmap) with white borders between tiles
  scale_fill_gradient(low = "white", high = "steelblue") + # Color scale from white (low) to steelblue (high proportion)
  labs(x = paste0("NMF Cluster (", target_cell_type_name, ")"),
       y = "Cell Type",
       fill = "Proportion",
       title = paste0("NMF Cluster Composition Heatmap for ", target_cell_type_name, " Microenvironment")) + # Clear labels
  theme_minimal(base_size = 14) # Clean, minimal theme
ggsave(paste0(clean_target_name, "_nmf_cluster_heatmap.pdf"), plot = plot_heatmap, width = 10, height = 7)
message("Cluster composition plots saved.")

# --- 10. Normalized Cell Count Visualization by Response and Cluster ---
# Visualize the normalized distribution of target cell clusters across different user-defined response groups.
# This helps investigate how different microenvironments are distributed across biological conditions.
message(paste0("Generating normalized cell count plot by '", response_column_name, "' and cluster..."))

# Check if the response_column_name exists in the Seurat object's metadata.
if (!response_column_name %in% colnames(Seurat_object@meta.data)) {
  warning(paste0("Column '", response_column_name, "' not found in Seurat object metadata. Skipping response-based visualization."))
} else {
  # Select relevant metadata columns (cluster, sample, response) and remove any rows with NA values.
  df <- Seurat_object@meta.data %>%
    select(cluster_nmf_env, sample_name, !!sym(response_column_name)) %>% # Use !!sym() for dynamic column selection
    na.omit() # Omit NA values from these specific columns to ensure valid counts

  # Calculate cell counts per Response group and NMF cluster.
  cell_counts <- df %>%
    group_by(!!sym(response_column_name), cluster_nmf_env) %>%
    summarise(cell_count = n(), .groups = "drop") # .groups = "drop" prevents grouped tibble output

  # Calculate the number of unique samples per Response group to normalize by.
  sample_counts <- df %>%
    group_by(!!sym(response_column_name)) %>%
    summarise(sample_n = n_distinct(sample_name), .groups = "drop")

  # Join counts and normalize cell counts by the number of samples in each response group.
  normalized_counts <- cell_counts %>%
    left_join(sample_counts, by = setNames(response_column_name, response_column_name)) %>% # Join by the response column
    mutate(count_per_sample = cell_count / sample_n) # Calculate normalized count

  # Create a 'group' variable for x-axis labeling using the response column's values.
  normalized_counts$group <- normalized_counts[[response_column_name]]

  # Convert 'cluster_nmf_env' to a factor for proper categorical plotting.
  normalized_counts$target_cell_cluster <- factor(normalized_counts$cluster_nmf_env)

  # Define a custom color palette for the stacked bar plot.
  set3_colors <- c(brewer.pal(n = 12, name = "Set3"),
                   '#D2691E', '#B39DDB', '#FFAB91', '#4DB6AC', '#F48FB1', '#C5E1A5', '#FFD54F') # Extended Set3 colors

  ### Stacked Bar Plot of Normalized Target Cell Counts
  # Visualizes the normalized counts of target cell clusters, stacked by Response group,
  # providing insight into cluster prevalence across conditions.
  plot_normalized_counts <- ggplot(normalized_counts, aes(x = group, y = count_per_sample, fill = target_cell_cluster)) +
    geom_bar(stat = "identity", position = "stack") + # Creates a stacked bar plot
    scale_fill_manual(values = set3_colors) + # Apply the custom color palette
    labs(x = response_column_name,
         y = "Normalized Cell Count per Sample",
         title = paste0("Normalized ", target_cell_type_name, " Cluster Counts by ", response_column_name),
         fill = paste0(target_cell_type_name, " Cluster")) + # Clear plot labels
    theme_minimal(base_size = 14) + # Clean, minimal theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis text for readability
  ggsave(paste0(clean_target_name, "_normalized_counts_by_", gsub(" ", "_", tolower(response_column_name)), ".pdf"), plot = plot_normalized_counts, width = 10, height = 7)
  message("Normalized cell count plot saved.")
}
message("Script execution complete.")
