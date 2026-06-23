##################################################
# 01_scRNA_functions.R
# BCM pAML manuscript
#
# Single-cell preprocessing and integration
##################################################
source("00_helper.R")

#' Create Seurat Object from Count Matrix
#'
#' @param counts Raw gene expression count matrix
#' @param sample_id Sample or project identifier
#' @param min_cells Minimum number of cells required per feature
#' @param min_features Minimum number of features per cell
#'
#' @return A Seurat object
#' @export
create_sc_object <- function(counts,
                             sample_id,
                             min_cells = 3,
                             min_features = 200) {
  
  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = min_cells,
    min.features = min_features
  )
  
  return(seu)
}


############################################################

#' Basic Single-Cell Preprocessing Pipeline
#'
#' @param seu A Seurat object
#' @param npcs Number of principal components to compute
#' @param resolution Clustering resolution
#'
#' @return A processed Seurat object
#' @export
run_sc_pipeline <- function(seu,
                            npcs = 30,
                            resolution = 2.5) {
  
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seu,
    pattern = "^MT-"
  )
  
  seu <- Seurat::NormalizeData(seu)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000)
  seu <- Seurat::ScaleData(seu)
  seu <- Seurat::RunPCA(seu, npcs = npcs)
  seu <- Seurat::FindNeighbors(seu, dims = 1:npcs)
  seu <- Seurat::FindClusters(seu, resolution = resolution)
  seu <- Seurat::RunUMAP(seu, dims = 1:npcs)
  
  return(seu)
}


############################################################

#' Integrate Multiple Patient Samples (Seurat Integration)
#'
#' @param object_list List of Seurat objects
#'
#' @return Integrated Seurat object
#' @export
integrate_patient_samples <- function(object_list) {
  
  # Normalize + variable features per object
  object_list <- lapply(object_list, function(x) {
    x <- Seurat::NormalizeData(x)
    x <- Seurat::FindVariableFeatures(x)
    return(x)
  })
  
  # Select integration features
  features <- Seurat::SelectIntegrationFeatures(object_list)
  
  # Find anchors
  anchors <- Seurat::FindIntegrationAnchors(
    object_list,
    anchor.features = features
  )
  
  # Integrate data
  integrated <- Seurat::IntegrateData(anchors)
  
  return(integrated)
}


############################################################

#' Assign Cell Types Using SingleR Annotation
#'
#' @param seu A Seurat object
#' @param reference SingleR reference dataset
#'
#' @return Seurat object with added `celltype` metadata
#' @export
assign_cell_types <- function(seu,
                              reference) {
  
  pred <- SingleR::SingleR(
    test = Seurat::GetAssayData(seu, slot = "data"),
    ref = reference,
    labels = reference$label.main
  )
  
  seu$celltype <- pred$labels
  
  return(seu)
}


############################################################

#' Collapse Highly Correlated Clusters
#'
#' This function identifies highly correlated cluster groups and merges them
#' based on a correlation threshold using a graph-like grouping strategy.
#'
#' @param correlation.matrix A symmetric numeric correlation matrix
#' @param correlation.threshold Numeric. Threshold for merging clusters
#'
#' @return A data.frame mapping original clusters to new grouped labels
#'
#' @details
#' Clusters with correlation above the threshold are grouped together.
#' This is used to collapse redundant or highly similar cell states.
#'
#' @export
regroup.cor <- function(correlation.matrix, correlation.threshold = 0.95) {
  
  # -----------------------------
  # Identify high correlations
  # -----------------------------
  top.cor <- which(
    abs(correlation.matrix) >= correlation.threshold &
      row(correlation.matrix) < col(correlation.matrix),
    arr.ind = TRUE
  )
  
  if (nrow(top.cor) == 0) {
    return(data.frame(new.group = character(0)))
  }
  
  # -----------------------------
  # Build cluster pairs
  # -----------------------------
  high_cor <- matrix(colnames(correlation.matrix)[top.cor], ncol = 2)
  
  groups <- list()
  count <- 0
  
  # Build initial groups
  for (elem in unique(high_cor[,1])) {
    count <- count + 1
    items <- high_cor[high_cor[,1] %in% elem, ]
    items <- unique(as.vector(items))
    groups[[paste0("group_", count)]] <- items
  }
  
  for (elem in unique(high_cor[,2])) {
    count <- count + 1
    items <- high_cor[high_cor[,2] %in% elem, ]
    items <- unique(as.vector(items))
    groups[[paste0("group_", count)]] <- items
  }
  
  # -----------------------------
  # Resolve overlaps (keep largest groups first)
  # -----------------------------
  groups <- groups[order(sapply(groups, length), decreasing = TRUE)]
  
  new.groups <- list()
  used <- character(0)
  count <- 0
  
  for (i in names(groups)) {
    if (i %in% used) next
    
    count <- count + 1
    new.groups[[paste0("group_", count)]] <- unique(unlist(groups[[i]]))
    used <- c(used, i)
  }
  
  # -----------------------------
  # Convert to annotation table
  # -----------------------------
  annotation <- stack(new.groups)
  colnames(annotation) <- c("cluster", "new.group")
  rownames(annotation) <- annotation$cluster
  
  return(annotation[, "new.group", drop = FALSE])
}


#' Regroup Seurat Clusters Based on Correlation Structure
#'
#' This function aggregates cluster-level expression profiles,
#' computes correlation between clusters, and merges highly similar clusters.
#'
#' @param object A Seurat object
#' @param cluster Character. Metadata column containing cluster labels
#' @param patient Character. Patient/sample name (used for outputs)
#' @param threshold Numeric. Correlation threshold for merging
#' @param assay Character. Seurat assay name (default: "RNA")
#' @param slot Character. Expression slot (default: "data")
#' @param save_outputs Logical. Whether to save plots and tables
#'
#' @return A Seurat object with updated `cellType` metadata
#'
#' @export
regrouping <- function(object,
                       cluster = "ClusterID",
                       patient = "Sample",
                       threshold = 0.986,
                       assay = "RNA",
                       slot = "data",
                       save_outputs = TRUE) {
  
  # -----------------------------
  # Extract cluster labels
  # -----------------------------
  cellType <- object@meta.data[[cluster]]
  
  groups <- split(seq_along(cellType), cellType)
  
  # -----------------------------
  # Extract expression matrix safely
  # -----------------------------
  expr <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # -----------------------------
  # Compute cluster means
  # -----------------------------
  cluster_means <- lapply(groups, function(x) {
    Matrix::rowMeans(expr[, x, drop = FALSE])
  })
  
  cluster_matrix <- do.call(cbind, cluster_means)
  
  # -----------------------------
  # Correlation matrix
  # -----------------------------
  correlation.matrix <- cor(cluster_matrix)
  
  annotation <- regroup.cor(correlation.matrix, threshold)
  
  # -----------------------------
  # Save outputs (optional)
  # -----------------------------
  if (save_outputs) {
    
    pdf(paste0(patient, "_cluster_correlation.pdf"),
        width = max(7, ncol(cluster_matrix) / 2),
        height = max(7, ncol(cluster_matrix) / 2))
    
    pheatmap::pheatmap(
      correlation.matrix,
      annotation = annotation,
      scale = "none",
      display_numbers = TRUE,
      number_format = "%.2f"
    )
    
    dev.off()
    
    write.csv(correlation.matrix,
              paste0(patient, "_correlation_matrix.csv"))
  }
  
  # -----------------------------
  # Update Seurat metadata
  # -----------------------------
  object$cellType <- object[[cluster]]
  
  counter <- 1
  
  if (nrow(annotation) > 0) {
    for (grp in unique(annotation$new.group)) {
      cells <- names(annotation$new.group[annotation$new.group == grp])
      object$cellType[object$cellType %in% cells] <- paste0("Regrouped.", counter)
      counter <- counter + 1
    }
  }
  
  # -----------------------------
  # Final consistency check
  # -----------------------------
  stopifnot(length(object$cellType) == ncol(expr))
  
  # -----------------------------
  # Return object
  # -----------------------------
  return(object)
}
