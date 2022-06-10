

#' @title Imputation under MCAR and MNAR hypothesis
#' 
#' @description 
#' this function performs missing values imputation under MCAR and MNAR hypothesis
#' 
#' @param dataSet.mvs expression matrix containing abundances with MVs 
#'                         (either peptides or proteins)
#' @param model.selector    - binary vector; "1" indicates MCAR proteins
#' @param method.MAR        - the method to be used for MAR missing values
#'                       - possible values: MLE (default), SVD, KNN
#' @param method.MNAR       - the method to be used for MAR missing values
#                       - possible values: QRILC (default), MinDet, MinProb
#'
#' @return dataset containing complete abundances
#'  
#' @export
#'
impute.MAR.MNAR = function(dataSet.mvs, model.selector, 
                           method.MAR = "KNN", method.MNAR = "QRILC"){
  
  # ___________________________________________________________________________________
  # perform MAR imputation using the specified method
  # -----------------------------------------------------------------------------------
  
    switch(method.MAR,
           MLE = {
             dataSet.MCAR.imputed = impute.MAR(dataSet.mvs,
                                                model.selector,
                                                method = "MLE") 
           },
           
           SVD = {
             dataSet.MCAR.imputed = impute.MAR(dataSet.mvs,
                                                model.selector,
                                                method = "SVD") 
           },
           
           KNN = {
             dataSet.MCAR.imputed = impute.MAR(dataSet.mvs, 
                                                model.selector,
                                                method = "KNN") 
           }     
    )
  
  # ___________________________________________________________________________________
  # perform MAR imputation using the specified method
  # -----------------------------------------------------------------------------------
  
    switch(method.MNAR,
           QRILC = {
             dataSet.complete.obj = impute.QRILC(dataSet.MCAR.imputed,tune.sigma = 0.3) 
             dataSet.complete = dataSet.complete.obj[[1]]
           },
           
           MinDet = {
             dataSet.complete = impute.MinDet(dataSet.MCAR.imputed) 
           },
           
           MinProb = {
             dataSet.complete = impute.MinProb(dataSet.MCAR.imputed) 
           }     
    )
  
  return(dataSet.complete)
  
}