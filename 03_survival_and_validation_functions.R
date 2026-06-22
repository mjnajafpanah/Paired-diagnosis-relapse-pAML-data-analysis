##################################################
# 03_survival_and_validation_functions.R
#
# Survival analyses and Bulk validation analyses
##################################################


run_cox_model <- function(
    data,
    formula
){
  
  fit <-
    
    survival::coxph(
      
      formula,
      
      data=data
      
    )
  
  return(
    
    summary(
      fit
    )
    
  )
  
}



run_kaplan_meier <- function(
    data,
    formula
){
  
  fit <-
    
    survival::survfit(
      
      formula,
      
      data=data
      
    )
  
  return(fit)
  
}



adjust_pvalues <- function(
    pvalues,
    method="bonferroni"
){
  
  p.adjust(
    pvalues,
    method=method
  )
  
}


combine_pvalues <- function(
    pvalues
){
  
  statistic <-
    -2*sum(
      log(
        pvalues
      ),
      na.rm=TRUE
    )
  
  df <- 2*length(
    pvalues
  )
  
  pchisq(
    statistic,
    df,
    lower.tail=FALSE
  )
  
}



extract_filename_info <- function(
    filename
){
  
  filename <-
    basename(
      filename
    )
  
  info <-
    stringr::str_match(
      filename,
      "(\\\\w+)\\\\.(C\\\\d+)vC\\\\d+"
    )
  
  data.frame(
    
    patient_id=
      info[,2],
    
    cluster_id=
      info[,3]
    
  )
  
}



merge_correlated_clusters <- function(
    correlation_matrix,
    threshold=0.95
){
  
  idx <-
    which(
      
      abs(
        correlation_matrix
      ) >= threshold &
        
        row(
          correlation_matrix
        ) <
        
        col(
          correlation_matrix
        ),
      
      arr.ind=TRUE
      
    )
  
  if(
    nrow(idx)==0
  ){
    
    return(NULL)
    
  }
  
  pairs <-
    
    data.frame(
      
      cluster1=
        colnames(
          correlation_matrix
        )[idx[,1]],
      
      cluster2=
        colnames(
          correlation_matrix
        )[idx[,2]]
      
    )
  
  return(pairs)
  
}