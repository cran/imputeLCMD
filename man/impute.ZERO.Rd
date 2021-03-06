\name{impute.ZERO}
\alias{impute.ZERO}
\title{ Imputation of missing entries by \code{0}. }
\description{ This function performs the trivial imputation of missing values by \code{0}. Is is only used for comparison purposes. }
\usage{ impute.ZERO(dataSet.mvs) }
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing data. }
}
\value{ A complete expression data matrix with missing values imputed. }
\author{ Cosmin Lazar }
\seealso{\code{\link{impute.QRILC}},
\code{\link{impute.MinDet}},
\code{\link{impute.MinProb}}}
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
exprsData.imputed = impute.ZERO(exprsData.MD)

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (dataSet.mvs) 
{
    dataSet.imputed = dataSet.mvs
    dataSet.imputed[which(is.na(dataSet.mvs))] = 0
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
