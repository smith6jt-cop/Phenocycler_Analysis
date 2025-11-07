# Phenocycler Analysis R Functions
# R integration for advanced statistical analysis and visualization

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Load Phenocycler h5ad file into Seurat
#'
#' @param h5ad_path Path to h5ad file
#' @return Seurat object
load_phenocycler_h5ad <- function(h5ad_path) {
  # Check if anndata is available
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required. Install with: install.packages('anndata')")
  }
  
  # Read AnnData
  adata <- anndata::read_h5ad(h5ad_path)
  
  # Convert to Seurat
  counts <- t(adata$X)
  rownames(counts) <- adata$var_names
  colnames(counts) <- adata$obs_names
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = adata$obs
  )
  
  # Add spatial coordinates if available
  if ("spatial" %in% names(adata$obsm)) {
    spatial_coords <- adata$obsm$spatial
    colnames(spatial_coords) <- c("x", "y")
    seurat_obj@images <- list(
      coords = spatial_coords
    )
  }
  
  # Add dimensionality reductions
  if ("X_pca" %in% names(adata$obsm)) {
    seurat_obj[["pca"]] <- CreateDimReducObject(
      embeddings = adata$obsm$X_pca,
      key = "PC_"
    )
  }
  
  if ("X_umap" %in% names(adata$obsm)) {
    seurat_obj[["umap"]] <- CreateDimReducObject(
      embeddings = adata$obsm$X_umap,
      key = "UMAP_"
    )
  }
  
  return(seurat_obj)
}

#' Differential expression analysis
#'
#' @param seurat_obj Seurat object
#' @param group_by Metadata column for grouping
#' @param ident_1 First group
#' @param ident_2 Second group
#' @param test_method Test method (wilcox, t, MAST, etc.)
#' @return Data frame with DE results
run_de_analysis <- function(seurat_obj, group_by, ident_1, ident_2, 
                            test_method = "wilcox") {
  Idents(seurat_obj) <- group_by
  
  de_results <- FindMarkers(
    seurat_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_method,
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  de_results$gene <- rownames(de_results)
  de_results <- de_results %>%
    mutate(
      significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
      direction = ifelse(avg_log2FC > 0, "Up", "Down")
    )
  
  return(de_results)
}

#' Create volcano plot
#'
#' @param de_results DE results data frame
#' @param title Plot title
#' @return ggplot object
plot_volcano <- function(de_results, title = "Volcano Plot") {
  de_results <- de_results %>%
    mutate(
      log_pval = -log10(p_val_adj + 1e-300),
      label = ifelse(
        significant & abs(avg_log2FC) > 1 & log_pval > 20,
        gene, ""
      )
    )
  
  p <- ggplot(de_results, aes(x = avg_log2FC, y = log_pval)) +
    geom_point(aes(color = significant), alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "red", alpha = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", 
               color = "blue", alpha = 0.5) +
    scale_color_manual(values = c("gray", "red")) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add labels for top genes
  if (sum(de_results$label != "") > 0) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20
    )
  }
  
  return(p)
}

#' Plot spatial features
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param features Features to plot
#' @return ggplot object
plot_spatial_features <- function(seurat_obj, features) {
  if (!"coords" %in% names(seurat_obj@images)) {
    stop("No spatial coordinates found in Seurat object")
  }
  
  coords <- seurat_obj@images$coords
  
  plots <- lapply(features, function(feat) {
    expr <- GetAssayData(seurat_obj, layer = "data")[feat, ]
    
    df <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      expression = expr
    )
    
    ggplot(df, aes(x = x, y = y, color = expression)) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_viridis_c() +
      labs(title = feat, color = "Expression") +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      )
  })
  
  return(wrap_plots(plots, ncol = min(3, length(features))))
}

#' Plot cell type composition
#'
#' @param seurat_obj Seurat object
#' @param celltype_col Cell type column
#' @param group_by Grouping variable
#' @return ggplot object
plot_composition <- function(seurat_obj, celltype_col, group_by = NULL) {
  meta <- seurat_obj@meta.data
  
  if (is.null(group_by)) {
    # Simple bar plot
    comp <- meta %>%
      group_by(!!sym(celltype_col)) %>%
      summarise(count = n()) %>%
      mutate(percentage = count / sum(count) * 100)
    
    p <- ggplot(comp, aes(x = reorder(!!sym(celltype_col), -percentage), 
                          y = percentage)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(
        title = "Cell Type Composition",
        x = "Cell Type",
        y = "Percentage"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    # Grouped bar plot
    comp <- meta %>%
      group_by(!!sym(group_by), !!sym(celltype_col)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(!!sym(group_by)) %>%
      mutate(percentage = count / sum(count) * 100)
    
    p <- ggplot(comp, aes(x = !!sym(celltype_col), y = percentage, 
                          fill = !!sym(group_by))) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(
        title = "Cell Type Composition by Group",
        x = "Cell Type",
        y = "Percentage",
        fill = group_by
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}

#' Statistical testing for composition differences
#'
#' @param seurat_obj Seurat object
#' @param celltype_col Cell type column
#' @param group_col Group column
#' @return Data frame with test results
test_composition_differences <- function(seurat_obj, celltype_col, group_col) {
  meta <- seurat_obj@meta.data
  
  # Create contingency table
  cont_table <- table(meta[[celltype_col]], meta[[group_col]])
  
  # Chi-square test
  chi_test <- chisq.test(cont_table)
  
  # Fisher's exact test (if small counts)
  fisher_test <- fisher.test(cont_table, simulate.p.value = TRUE)
  
  results <- list(
    chi_square = chi_test,
    fisher_exact = fisher_test,
    contingency_table = cont_table
  )
  
  return(results)
}

# Print package loading message
cat("Phenocycler R analysis functions loaded successfully!\n")
cat("Available functions:\n")
cat("  - load_phenocycler_h5ad()\n")
cat("  - run_de_analysis()\n")
cat("  - plot_volcano()\n")
cat("  - plot_spatial_features()\n")
cat("  - plot_composition()\n")
cat("  - test_composition_differences()\n")
