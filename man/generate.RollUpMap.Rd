\name{generate.RollUpMap}
\alias{generate.RollUpMap}
\title{ Generates peptide to protein map. }
\description{ Generates peptide to protein map. For a given peptide expression matrix with \code{nPep} peptides, this functions creates a random map to be used for the aggregation of the peptides into \code{nProt} proteins. }
\usage{ generate.RollUpMap(nProt, pep.Expr.Data) }
\arguments{
  \item{nProt}{ Number of proteins }
  \item{pep.Expr.Data}{ Peptide expression data matrix }
}
\value{ A vector of the length \code{nPep} containing integer random samples in the range \code{[1, nProt]}}
\author{ Cosmin Lazar }
\seealso{\code{\link{generate.ExpressionData}},
\code{\link{pep2prot}}
}
\examples{
  
exprsDataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                          meanSamples = 0, sdSamples = 0.2,
                          nFeatures = 1000, nFeaturesUp = 50, nFeaturesDown = 50,
                          meanDynRange = 20, sdDynRange = 1,
                          meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = exprsDataObj[[1]]
rollUpMap = generate.RollUpMap(round(dim(exprsData)[1]/2),exprsData)

## The function is currently defined as
function (nProt, pep.Expr.Data) 
{
    n = dim(pep.Expr.Data)[1]
    temp = 1:nProt
    pep.prot.Map = rep(0, n)
    pep.prot.Map[sample(temp)] = sample(temp)
    pep.prot.Map[which(pep.prot.Map == 0)] = sample.int(nProt, 
        size = (n - nProt), replace = T)
    return(pep.prot.Map)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
