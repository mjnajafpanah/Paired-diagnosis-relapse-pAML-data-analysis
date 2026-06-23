##################################################
# 02_pAML_identification_functions.R
#
# pAML cluster identification
##################################################
source("00_helper.R")

#' Calculate Cluster Fractions Across Samples
#'
#' This function computes the percentage composition of each cluster
#' within each sample group (e.g., diagnosis vs relapse).
#'
#' @param metadata A data.frame containing cell metadata.
#' @param cluster_col Column name for cluster labels.
#' @param sample_col Column name for sample/condition grouping.
#'
#' @return A matrix of cluster proportions (%) per sample.
#'
#' @export
calculate_cluster_fractions <- function(metadata,
                                        cluster_col = "seurat_clusters",
                                        sample_col = "condition") {
  
  # -----------------------------
  # Build contingency table
  # -----------------------------
  tab <- table(
    metadata[[sample_col]],
    metadata[[cluster_col]]
  )
  
  # -----------------------------
  # Convert to row-wise proportions
  # -----------------------------
  prop <- prop.table(tab, margin = 1)
  
  # Convert to percentage
  prop <- round(prop * 100, 2)
  
  return(prop)
}


############################################################

#' Classify pAML Cluster Dynamics
#'
#' This function classifies clusters into biological categories
#' based on changes in abundance between conditions (e.g., diagnosis vs relapse).
#'
#' @param cluster_table A matrix/data.frame of cluster percentages
#'   (output of `calculate_cluster_fractions`).
#' @param threshold Numeric cutoff (%) to define meaningful presence.
#'
#' @return A list containing:
#'   - expanded: clusters enriched in relapse
#'   - diminished: clusters enriched in diagnosis
#'   - stable: clusters present in both conditions
#'   - relapse_specific: clusters unique to relapse
#'
#' @export
classify_paml_clusters <- function(cluster_table,
                                   threshold = 10) {
  
  # -----------------------------
  # Expanded in relapse
  # -----------------------------
  expanded <- names(which(
    cluster_table["Relapse", ] > threshold &
      cluster_table["Diagnosis", ] <= threshold
  ))
  
  # -----------------------------
  # Diminished in relapse
  # -----------------------------
  diminished <- names(which(
    cluster_table["Diagnosis", ] > threshold &
      cluster_table["Relapse", ] <= threshold
  ))
  
  # -----------------------------
  # Stable clusters
  # -----------------------------
  stable <- names(which(
    cluster_table["Diagnosis", ] > threshold &
      cluster_table["Relapse", ] > threshold
  ))
  
  # -----------------------------
  # Relapse-specific clusters
  # -----------------------------
  relapse_specific <- names(which(
    cluster_table["Diagnosis", ] == 0 &
      cluster_table["Relapse", ] > threshold
  ))
  
  # -----------------------------
  # Return structured output
  # -----------------------------
  return(list(
    expanded = expanded,
    diminished = diminished,
    stable = stable,
    relapse_specific = relapse_specific
  ))
}

