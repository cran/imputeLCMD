\name{pep2prot}
\alias{pep2prot}
\title{ Peptide to protein mapping. }
\description{ This function performs peptide to protein mapping given a peptide expression matrix and a peptide to protein map. }
\usage{ pep2prot(pep.Expr.Data, rollup.map) }
\arguments{
  \item{pep.Expr.Data}{ Peptide expression data matrix. }
  \item{rollup.map}{ Peptide to protein map. }
}
\value{ A protein expression matrix. Each line is obtained by aggregating the lines in the peptide expression matrix which have the same value in the rollup.map. }
\author{ Cosmin Lazar }
\seealso{ \code{\link{pep2prot}}, \code{\link{generate.RollUpMap}}}
\examples{

# generate expression data matrix
dataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                          meanSamples = 0, sdSamples = 0.2,
                          nFeatures = 2000, nFeaturesUp = 100, nFeaturesDown = 100,
                          meanDynRange = 20, sdDynRange = 1,
                          meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = dataObj[[1]]

# insert 15\% missing data with 100\% missing not at random
m.THR = quantile(exprsData, probs = 0.15)
sd.THR = 0.1
MNAR.rate = 100
exprsData.MD.obj = insertMVs(exprsData,m.THR,sd.THR,MNAR.rate)
exprsData.MD = exprsData.MD.obj[[2]]

# generate rollup.map
rollUpMap = generate.RollUpMap(round(dim(exprsData.MD)[1]/2),exprsData.MD)

# peptide to protein mapping
protExprsData = pep2prot(exprsData.MD,rollUpMap)

## The function is currently defined as
function (pep.Expr.Data, rollup.map) 
{
    protIDs = unique(rollup.map)
    prot.Expr.Data = matrix(, nrow = length(protIDs), ncol = dim(pep.Expr.Data)[2])
    for (i in 1:length(protIDs)) {
        temp = protIDs[i]
        pepSet = pep.Expr.Data[which(rollup.map == temp), ]
        if (length(which(rollup.map == temp)) == 1) {
            prot.Expr.Data[i, ] = pepSet
        }
        else {
            prot.Expr.Data[i, ] = apply(pepSet, 2, median, na.rm = T)
        }
    }
    return(prot.Expr.Data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
