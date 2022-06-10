#' @title Imputation with min value
#' 
#' @description this function performs missing values imputation by the minimum value observed
#' 
#' @param dataSet.mvs expression matrix with MVs (either peptides or proteins)
#' @param q the q quantile used to estimate the minimum 
#                               value observed for each sample
#' 
#' @return dataset containing complete abundances
#'  
#' @export
#' 
#' @import stats
#' 

impute.MinDet = function(dataSet.mvs,q = 0.01){
  
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  
  lowQuantile.samples = apply(dataSet.imputed,2,quantile,prob = q,na.rm = T)
  
  for (i in 1:(nSamples)){
    dataSet.imputed[which(is.na(dataSet.mvs[,i])),i] = lowQuantile.samples[i]
  }
    
  return(dataSet.imputed)
  
}