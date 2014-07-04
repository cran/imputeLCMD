\name{impute.wrapper.MLE}
\alias{impute.wrapper.MLE}
\title{ MLE-based imputation of missing data. }
\description{ Is is a wrapper function that performs the imputation of missing data using EM algorithm. The wrapper is built around the \code{\link{imp.norm}} function from the \href{http://cran.r-project.org/web/packages/norm/norm.pdf}{norm} package. }
\usage{ impute.wrapper.MLE(dataSet.mvs) }
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing data. }
}
\value{ A complete expression data matrix with missing values imputed. }
\references{ See package \href{http://cran.r-project.org/web/packages/norm/norm.pdf}{norm}. }
\author{ Cosmin Lazar }
\seealso{\code{\link{impute.wrapper.KNN}},
\code{\link{impute.wrapper.SVD}} }
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
  
exprsData.imputed = impute.wrapper.MLE(exprsData.MD)

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (dataSet.mvs) 
{
    s <- prelim.norm(dataSet.mvs)
    thetahat <- em.norm(s, showits = FALSE)
    rngseed(1234567)
    dataSet.imputed <- imp.norm(s, thetahat, dataSet.mvs)
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
