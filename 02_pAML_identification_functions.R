##################################################
# 02_pAML_identification_functions.R
#
# pAML cluster identification
##################################################


calculate_cluster_fractions <- function(
    metadata,
    cluster_col="seurat_clusters",
    sample_col="condition"
){
  
  tab <-
    table(
      metadata[[sample_col]],
      metadata[[cluster_col]]
    )
  
  prop <- prop.table(
    tab,
    margin=1
  )
  
  prop <- round(
    prop*100,
    2
  )
  
  return(prop)
  
}


classify_paml_clusters <- function(
    cluster_table,
    threshold=10
){
  
  expanded <-
    names(
      which(
        cluster_table["Relapse",] > threshold &
          cluster_table["Diagnosis",] <= threshold
      )
    )
  
  diminished <-
    names(
      which(
        cluster_table["Diagnosis",] > threshold &
          cluster_table["Relapse",] <= threshold
      )
    )
  
  stable <-
    names(
      which(
        cluster_table["Diagnosis",] > threshold &
          cluster_table["Relapse",] > threshold
      )
    )
  
  relapse_specific <-
    names(
      which(
        cluster_table["Diagnosis",]==0 &
          cluster_table["Relapse",] > threshold
      )
    )
  
  list(
    
    expanded=expanded,
    
    diminished=diminished,
    
    stable=stable,
    
    relapse_specific=relapse_specific
    
  )
  
}