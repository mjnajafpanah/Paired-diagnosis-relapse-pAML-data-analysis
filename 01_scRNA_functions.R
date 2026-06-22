##################################################
# 01_scRNA_functions.R
# BCM pAML manuscript
#
# Single-cell preprocessing and integration
##################################################

#' Create Seurat object
create_sc_object <- function(
    counts,
    sample_id,
    min_cells=3,
    min_features=200
){
  
  seu <- Seurat::CreateSeuratObject(
    counts=counts,
    project=sample_id,
    min.cells=min_cells,
    min.features=min_features
  )
  
  return(seu)
  
}


#' Basic single-cell pipeline
run_sc_pipeline <- function(
    seu,
    npcs=30,
    resolution=2.5
){
  
  seu[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(
      seu,
      pattern="^MT-"
    )
  
  seu <- Seurat::NormalizeData(seu)
  
  seu <- Seurat::FindVariableFeatures(
    seu,
    nfeatures=2000
  )
  
  seu <- Seurat::ScaleData(seu)
  
  seu <- Seurat::RunPCA(
    seu,
    npcs=npcs
  )
  
  seu <- Seurat::FindNeighbors(
    seu,
    dims=1:npcs
  )
  
  seu <- Seurat::FindClusters(
    seu,
    resolution=resolution
  )
  
  seu <- Seurat::RunUMAP(
    seu,
    dims=1:npcs
  )
  
  return(seu)
  
}


#' Integrate diagnosis, relapse and normal BM
integrate_patient_samples <- function(
    object_list
){
  
  object_list <-
    lapply(
      object_list,
      function(x){
        
        x <- Seurat::NormalizeData(x)
        
        x <- Seurat::FindVariableFeatures(x)
        
      }
    )
  
  features <-
    Seurat::SelectIntegrationFeatures(
      object_list
    )
  
  anchors <-
    Seurat::FindIntegrationAnchors(
      object_list,
      anchor.features=features
    )
  
  integrated <-
    Seurat::IntegrateData(
      anchors
    )
  
  return(integrated)
  
}


#' Assign cell types
assign_cell_types <- function(
    seu,
    reference
){
  
  pred <-
    SingleR::SingleR(
      test=Seurat::GetAssayData(
        seu,
        slot="data"
      ),
      
      ref=reference,
      
      labels=reference$label.main
    )
  
  seu$celltype <-
    pred$labels
  
  return(seu)
  
}