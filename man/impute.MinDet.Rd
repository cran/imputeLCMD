\name{impute.MinDet}
\alias{impute.MinDet}
\title{ Imputation of left-censored missing data using a deterministic minimal value approach. }
\description{ Performs the imputation of left-censored missing data using a deterministic minimal value approach.Considering a peptide/protein expression data matrix with \code{n} columns corresponding to biological samples and \code{p} lines corresponding to peptides/proteins, for each sample (column), the missing entries are replaced with a minimal value observed in that sample. The minimal value observed is estimated as being the \code{q-th} quantile (e.g. \code{q = 0.01}) of the observed values in that sample. }
\usage{ impute.MinDet(dataSet.mvs, q = 0.01) }
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing data. }
  \item{q}{ A scalar used to determine a low expression value to be used for missing data imputation. \code{0 < q < 1}, in this case \code{q} should be set to a low value. The default value is \code{q = 0.01}. }
}
\value{ A complete expression data matrix with missing values imputed. }
\author{ Cosmin Lazar }
\seealso{\code{\link{impute.QRILC}}, \code{\link{impute.MinProb}}, \code{\link{impute.ZERO}}}
\examples{
# generate expression data matrix
exprsDataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                          meanSamples = 0, sdSamples = 0.2,
                          nFeatures = 1000, nFeaturesUp = 50, nFeaturesDown = 50,
                          meanDynRange = 20, sdDynRange = 1,
                          meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = exprsDataObj[[1]]
  
# insert 15\% missing data with 100\% missing not at random

m.THR = quantile(exprsData, probs = 0.15)
sd.THR = 0.1
MNAR.rate = 100
exprsData.MD.obj = insertMVs(exprsData,m.THR,sd.THR,MNAR.rate)
exprsData.MD = exprsData.MD.obj[[2]]

# perform missing data imputation
  
exprsData.imputed = impute.MinDet(exprsData.MD)

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (dataSet.mvs, q = 0.01) 
{
    nSamples = dim(dataSet.mvs)[2]
    dataSet.imputed = dataSet.mvs
    lowQuantile.samples = apply(dataSet.imputed, 2, quantile, 
        prob = q, na.rm = T)
    for (i in 1:(nSamples)) {
        dataSet.imputed[which(is.na(dataSet.mvs[, i])), i] = lowQuantile.samples[i]
    }
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
