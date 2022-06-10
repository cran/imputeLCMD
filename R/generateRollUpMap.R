#' @title Generate roll up map
#' 
#' @description Tthis function generates a map for peptide to protein roll-up
#' 
#' @param nProt number of proteins to map to the peptide expression data
#' @param pep.Expr.Data matrix of peptide expression data
#' 
#' @return the peptide to protein map (for each row in 
#' pep.prot.Map the corresponding value corresponds 
#'  to the index of the protein that peptide is mapped to)
#'  
#' @export

  generate.RollUpMap = function(nProt, pep.Expr.Data){
    
    n = dim(pep.Expr.Data)[1]
    temp = 1:nProt
    pep.prot.Map = rep(0,n)
    pep.prot.Map[sample(temp)] = sample(temp)
    pep.prot.Map[which(pep.prot.Map==0)] = sample.int(nProt, 
                                                      size = (n - nProt), 
                                                      replace = T)
    
    return(pep.prot.Map)
  }