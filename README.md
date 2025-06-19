# Spatial Cell Microenvironment Analysis with Seurat

---

### **Unlocking the Secrets of Cellular Neighborhoods in Spatial Omics Data**

This repository presents a **robust and user-friendly R script** designed to dissect the intricate **cellular microenvironments** within spatial transcriptomics data. Leveraging the powerful [Seurat](https://satijalab.org/seurat/) framework, this tool empowers researchers to move beyond single-cell analysis and explore the **critical spatial interactions** that shape biological function and disease.

Our script focuses on identifying and characterizing the **unique cellular neighborhoods** surrounding a **target cell population** (e.g., T/NK cells, B cells, Macrophages). By quantifying the **proportional composition of neighboring cell types**, we can cluster these target cells based on their **distinct spatial contexts**, revealing novel biological insights.

---

### **Key Features**

* **Customizable Target Cell Analysis**: Easily define your cell of interest to investigate its specific microenvironment.
* **Flexible Annotation Integration**: Seamlessly integrate your existing cell type annotations.
* **Spatial Proximity Modeling**: Precisely define cell-cell interactions based on adjustable spatial distances.
* **Advanced Clustering**: Employ both K-means and Non-negative Matrix Factorization (NMF) to uncover hidden patterns in cellular neighborhoods.
* **Comprehensive Visualization**: Generate intuitive plots that illustrate cluster composition and distribution across experimental conditions.
* **Reproducible Workflow**: A clear, parameter-driven script ensures your analyses are consistent and shareable.

---

### **How It Works: The Analytical Journey**

The script guides your data through a structured analytical pipeline:

* **Data Ingestion & Preparation**: Your Seurat object is loaded, and **precise spatial coordinates are integrated** into the metadata, forming the foundation for spatial analysis.
* **Neighborhood Definition**: We **meticulously identify cells neighboring your target population** within a user-defined spatial radius, capturing localized interactions.
* **Proportional Profiling**: For each target cell, the script **quantifies the relative abundance of all neighboring cell types**, creating a unique "spatial signature."
* **Contextual Clustering**: Target cells are then **grouped using K-means and NMF** based on these spatial signatures, revealing populations with distinct microenvironmental roles.
* **Insightful Visualization**: The results are transformed into **clear, publication-ready plots** that illuminate the specific cell type compositions of each cluster and their variations across different conditions (e.g., disease states, treatment responses).

---

### **Getting Started: Plug & Play Your Data**

To use this script, simply follow these steps:

1.  **Prerequisites**:
    * Ensure you have **R** and **RStudio** installed.
    * Install all necessary R packages by running:
        ```R
        install.packages(c("Seurat", "dplyr", "SpatialExperiment", "SeuratObject", "sp", "spdep", "STutility", "tidyr", "ggplot2", "NMF", "RColorBrewer", "ggpubr"))
        ```

2.  **Configuration**:
    * Open the `analyze_spatial_microenvironment.R` (or your chosen script name) file.
    * Navigate to the **"User-Defined Parameters"** section at the very top of the script.
    * **Modify the variables** to match your specific dataset:
        * `seurat_object_path`: Path to your Seurat object (`.rds` file).
        * `assay_name`: The Seurat assay containing your data (e.g., `"Xenium"`, `"RNA"`).
        * `target_cell_type_name`: The exact name of your target cell type (e.g., `"T/NK cell"`, `"B cell"`).
        * `annotated_cell_type_column`: The metadata column containing your cell type annotations (e.g., `"SingleR_full"`, `"predicted_labels"`).
        * `neighbor_distance_um`: The maximum distance (in spatial units, e.g., micrometers) to consider for neighboring cells. Tune this value based on your data's resolution and biological question.
        * `num_clusters`: The desired number of clusters for K-means and NMF. Experimentation may be required to find the optimal number.
        * `response_column_name`: The metadata column representing sample groups or response variables (e.g., `"Response"`, `"patient_group"`).

3.  **Run the Script**:
    * Execute the entire script in RStudio.

4.  **Explore Results**:
    * Check your working directory for newly generated `.rds` files (containing updated metadata and proportion matrices) and `.pdf` plots (visualizations of cluster compositions and distributions).
    * Messages in your R console will guide you through the analysis progress.

---

### **Beyond the Code: Contributing and Support**

We welcome contributions and feedback! If you find this tool useful or have suggestions for improvement, please feel free to open an issue or submit a pull request.



