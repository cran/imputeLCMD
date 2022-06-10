#' @title Generate expression data
#' 
#' @description this function generates artificial peptide abundance data with DA proteins
#' samples are drawn from a gaussian distribution
#' 
#' @param nSamples1 number of samples in condition 1
#' @param nSamples2 number of samples in condition 2
#' @param meanSamples xxx 
#' @param sdSamples xxx
#' @param nFeatures number of total features
#' @param nFeaturesUp number of features up regulated
#' @param nFeaturesDown number of features down regulated
#' @param meanDynRange mean value of the dynamic range
#' @param sdDynRange sd of the dynamic range
#' @param meanDiffAbund xxx
#' @param sdDiffAbund xxx
#' 
#' @return A list containing the data, the conditions label and the regulation
#' label (up/down/no)
#' 
#' @export
#' 
#' @import stats


generate.ExpressionData = function(nSamples1,
                                   nSamples2,
                                   meanSamples, 
                                   sdSamples,
                                   nFeatures, 
                                   nFeaturesUp, 
                                   nFeaturesDown,
                                   meanDynRange, 
                                   sdDynRange,
                                   meanDiffAbund,
                                   sdDiffAbund){
  
  # generate a matrix of nSamples1 + nSamples2 samples from a Gaussian distribution
  nSamples = nSamples1 + nSamples2
  data = matrix(rnorm(nSamples*nFeatures,meanSamples,sdSamples),nFeatures,nSamples)
  
  # spread the data on the dynamic range ...
  means = rnorm(nFeatures,meanDynRange,sdDynRange)
  data = data + means
  
  # select the groups of samples ...
  conditions = c(rep(1,nSamples1),rep(2,nSamples2))
  
  # define the extra abundance values for up and down expressed features ...
  DE.coef.up = matrix(rnorm(nFeaturesUp*nSamples1,meanDiffAbund,sdDiffAbund),
                      nFeaturesUp,nSamples1)
  DE.coef.down = matrix(rnorm(nFeaturesDown*nSamples2,meanDiffAbund,sdDiffAbund),
                        nFeaturesDown,nSamples2)
  
  # create up and down expressed features
  data[1:nFeaturesUp,conditions==1] = DE.coef.up+data[1:nFeaturesUp,conditions==1]
  data[(nFeaturesUp+1):(nFeaturesUp + nFeaturesDown),conditions==2] =
    DE.coef.down+data[(nFeaturesUp+1):(nFeaturesUp + nFeaturesDown),conditions==2]
  
  # define the labels vector for the features indicating whether they are up/down/no expressed
  labelFeatures = c(rep(1,nFeaturesUp),
            rep(2,nFeaturesDown),
            rep(3,nFeatures - (nFeaturesUp+nFeaturesDown)))
  
  row.names(data) = 1:nFeatures
  
  return(list(data,conditions,labelFeatures))
  
}