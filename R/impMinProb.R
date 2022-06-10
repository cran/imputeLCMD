#' @title Imputation by random draws
#' 
#' @description 
#' This function performs missing values imputation by random draws from a gaussian
#   distribution centered in the minimum value observed and with standard deviation
#   equal to the median value of the population of line-wise standard deviations
#' 
#' @param dataSet.mvs expression matrix containing abundances with 
#' MVs (either peptides or proteins)
#' @param q the q-th quantile used to estimate the minimum 
#' value observed for each sample
#' @param tune.sigma coefficient that controls the sd of the MNAR distribution
#' 
#' @return dataset containing complete abundances
#'  
#' @export
#' 
#' @import stats
#' 
impute.MinProb = function(dataSet.mvs,q = 0.01,tune.sigma = 1){
    
  nSamples = dim(dataSet.mvs)[2]
  nFeatures = dim(dataSet.mvs)[1]
    
  dataSet.imputed = dataSet.mvs
  
  # - select the minimum values sample-wise (corresponding to the q-th quantile)
  min.samples = apply(dataSet.imputed,2,quantile,prob = q,na.rm = T)
  
  # - estimate protein-wise standard deviation using only proteins containing 
  #   more than 50% non-NAs
  
  count.NAs = apply(!is.na(dataSet.mvs),1,sum)
  count.NAs = count.NAs/nSamples
  
  dataSet.filtered = dataSet.mvs[which(count.NAs>0.5),]
    
  protSD = apply(dataSet.filtered,1,sd,na.rm = T)
  sd.temp = median(protSD,na.rm = T)*tune.sigma
  
  print(sd.temp)
  
  for (i in 1:(nSamples)){
    
    dataSet.to.impute.temp = rnorm(nFeatures,
                                   mean = min.samples[i],
                                   sd = sd.temp)
    dataSet.imputed[which(is.na(dataSet.mvs[,i])),i] = dataSet.to.impute.temp[which(is.na(dataSet.mvs[,i]))]
    
  }
  
  return(dataSet.imputed)
  
}