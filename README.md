Spatial Target Cell Microenvironment Analysis
This R script provides a comprehensive workflow for analyzing the cellular microenvironment around Target cells in spatial transcriptomics data, leveraging the Seurat object framework. It involves identifying neighboring cells, calculating cell type proportions in the vicinity of Target cells, and subsequently clustering these Target cells based on their unique microenvironmental profiles. The results are then visualized to reveal distinct Target cell clusters and their spatial associations with other cell types.

Workflow Overview
The script performs the following key steps:

Package Loading & Data Preparation: Loads necessary R packages and prepares the Seurat object by adding spatial coordinates to its metadata.
Neighbor Definition: Identifies neighboring cells around each Target cell within the same image based on spatial proximity.
Cell Type Proportion Calculation: Quantifies the proportions of different cell types in the neighborhood of each Target cell.
Target Cell Clustering: Applies K-means and Non-negative Matrix Factorization (NMF) to cluster Target cells based on their calculated neighborhood profiles.
Result Integration & Saving: Integrates the clustering results back into the Seurat object's metadata and saves intermediate data.
Visualization: Generates insightful plots to visualize the composition of each Target cell cluster and their distribution across different sample conditions.
