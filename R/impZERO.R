
#' @title Imputation by 0.
#' 
#' @description 
#' This function performs missing values imputation by 0.
#' 
#' @param dataSet.mvs expression matrix containing abundances with MVs 
#'                         (either peptides or proteins)
#'
#' @return dataset containing complete abundances
#'  
#' @export
#' 

impute.ZERO = function(dataSet.mvs){
  
  dataSet.imputed = dataSet.mvs
  dataSet.imputed[which(is.na(dataSet.mvs))] = 0
  
  return(dataSet.imputed)
  
}