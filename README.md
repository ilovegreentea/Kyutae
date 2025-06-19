Spatial Cell Microenvironment Analysis with Seurat
Unlocking the Secrets of Cellular Neighborhoods in Spatial Omics Data
This repository presents a robust and user-friendly R script designed to dissect the intricate cellular microenvironments within spatial transcriptomics data. Leveraging the powerful Seurat framework, this tool empowers researchers to move beyond single-cell analysis and explore the critical spatial interactions that shape biological function and disease.

Our script focuses on identifying and characterizing the unique cellular neighborhoods surrounding a target cell population (e.g., T/NK cells, B cells, Macrophages). By quantifying the proportional composition of neighboring cell types, we can cluster these target cells based on their distinct spatial contexts, revealing novel biological insights.

Key Features
Customizable Target Cell Analysis: Easily define your cell of interest to investigate its specific microenvironment.
Flexible Annotation Integration: Seamlessly integrate your existing cell type annotations.
Spatial Proximity Modeling: Precisely define cell-cell interactions based on adjustable spatial distances.
Advanced Clustering: Employ both K-means and Non-negative Matrix Factorization (NMF) to uncover hidden patterns in cellular neighborhoods.
Comprehensive Visualization: Generate intuitive plots that illustrate cluster composition and distribution across experimental conditions.
Reproducible Workflow: A clear, parameter-driven script ensures your analyses are consistent and shareable.
How It Works: The Analytical Journey
The script guides your data through a structured analytical pipeline:

Data Ingestion & Preparation: Your Seurat object is loaded, and precise spatial coordinates are integrated into the metadata, forming the foundation for spatial analysis.
Neighborhood Definition: We meticulously identify cells neighboring your target population within a user-defined spatial radius, capturing localized interactions.
Proportional Profiling: For each target cell, the script quantifies the relative abundance of all neighboring cell types, creating a unique "spatial signature."
Contextual Clustering: Target cells are then grouped using K-means and NMF based on these spatial signatures, revealing populations with distinct microenvironmental roles.
Insightful Visualization: The results are transformed into clear, publication-ready plots that illuminate the specific cell type compositions of each cluster and their variations across different conditions (e.g., disease states, treatment responses).
