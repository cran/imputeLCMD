#' @title imputation based on SVD algorithm
#' 
#' @description 
#' this function performs missing values imputation based on SVD algorithm
#' 
#' @param dataSet.mvs expression matrix with MVs (either peptides or proteins)
#' @param K the number of PCs
#'
#' @return expression matrix with MVs imputed
#'  
#' @export
#' @importFrom pcaMethods pca
#'
impute.wrapper.SVD = function(dataSet.mvs,K){
  
  resultSVD = pca(dataSet.mvs, method="svdImpute", nPcs = K)
  dataSet.imputed = resultSVD@completeObs
  
  return(dataSet.imputed)
  
}