
#' @title imputation using the EM algorithm
#' 
#' @description 
#' This function performs missing values imputation using the EM algorithm
#' 
#' @param dataSet.mvs expression matrix with MVs (either peptides or proteins)
#'
#' @return expression matrix with MVs imputed
#'  
#' @export
#' @importFrom norm prelim.norm em.norm rngseed imp.norm
#' 
impute.wrapper.MLE = function(dataSet.mvs){
  
  s <- prelim.norm (dataSet.mvs)
  thetahat <- em.norm (s, showits = FALSE ) # find the mle
  rngseed (1234567)
  dataSet.imputed <- imp.norm (s, thetahat , dataSet.mvs ) # impute missing data under the MLE
  
  return(dataSet.imputed)
  
}