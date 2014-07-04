\name{imputeLCMD-package}
\alias{imputeLCMD-package}
\alias{imputeLCMD}
\docType{package}
\title{ A colection of methods for left-censored missing data imputation }
\description{ The package contains a collection of functions for left-censored missing data imputation. Left-censoring is a special case of missing not at random (MNAR)  mechanism that generates non-responses in proteomics experiments. The package also contains functions to artificially generate peptide/protein expression data (log-transformed) as random draws from a multivariate Gaussian distribution as well as a function to generate missing data (both randomly and non-randomly). For comparison reasons, the package also contains several wrapper functions for the imputation of non-responses that are missing at random. }
\details{
\tabular{ll}{
Package: \tab imputeLCMD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-04\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Cosmin Lazar

Maintainer: Cosmin Lazar <vcosminlazar@gmail.com>
}
\keyword{ package }
\seealso{ \code{\link{impute.QRILC}}, \code{\link{impute.MinDet}}, \code{\link{impute.MinProb}}}
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
obj.QRILC = impute.QRILC(exprsData.MD)
exprsData.imputed = obj.QRILC[[1]]
  
\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}
}
