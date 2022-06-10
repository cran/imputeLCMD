

#' @title peptide to protein roll-up
#' 
#' @description 
#' this function performs peptide to protein roll-up
#' 
#' @param pep.Expr.Data matrix of peptide expression data
#' @param rollup.map the map to peptide to protein mapping
#'
#' @return matrix of peptide expression data
#'  
#' @export
#' 
#' @import stats
#'
pep2prot = function(pep.Expr.Data, rollup.map){
  
  protIDs = unique(rollup.map)
  prot.Expr.Data = matrix( ,
                           nrow = length(protIDs),
                           ncol = dim(pep.Expr.Data)[2]
                           )
  
  for (i in 1:length(protIDs)){
    
    temp = protIDs[i]
    pepSet = pep.Expr.Data[which(rollup.map==temp),]
    
    if (length(which(rollup.map==temp))==1){
      prot.Expr.Data[i,] = pepSet
    }
    else{
      prot.Expr.Data[i,] = apply(pepSet,2,median,na.rm = T)
    }
    
  }
    
  return(prot.Expr.Data)
}