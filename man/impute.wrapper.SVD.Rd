\name{impute.wrapper.SVD}
\alias{impute.wrapper.SVD}
\title{ SVD-based imputation. }
\description{ This is a wrapper function that performs SVD-based imputation of missing data. The wrapper is built around the \code{\link{pca}} function from the \href{http://master.bioconductor.org/packages/devel/bioc/vignettes/pcaMethods/inst/doc/missingValues.pdf}{pcaMethods} package. }
\usage{ impute.wrapper.SVD(dataSet.mvs, K) }
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing data. }
  \item{K}{ The number of \code{PCs} used. }
}
\value{ A complete expression data matrix with missing values imputed. }
\references{ See package \href{http://master.bioconductor.org/packages/devel/bioc/vignettes/pcaMethods/inst/doc/missingValues.pdf}{pcaMethods} }
\author{ Cosmin Lazar }
\seealso{\code{\link{pca}}, \code{\link{impute.wrapper.KNN}}, \code{\link{impute.wrapper.MLE}}}
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
exprsData.imputed = impute.wrapper.SVD(exprsData.MD,2)

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (dataSet.mvs, K) 
{
    resultSVD = pca(dataSet.mvs, method = "svdImpute", nPcs = K)
    dataSet.imputed = resultSVD@completeObs
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
