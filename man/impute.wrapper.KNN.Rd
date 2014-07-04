\name{impute.wrapper.KNN}
\alias{impute.wrapper.KNN}
\title{ Wrapper function for KNN imputation. }
\description{ This function is a wrapper around the \code{\link{impute.knn}} function from \href{http://www.bioconductor.org/packages/release/bioc/manuals/impute/man/impute.pdf}{impute} package that performs KNN imputation for a data matrix containing missing entries. }
\usage{ impute.wrapper.KNN(dataSet.mvs, K) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing entries. }
  \item{K}{T he number of neighbors used to infer the missing data. }
}
\details{ Requires \href{http://www.bioconductor.org/packages/release/bioc/manuals/impute/man/impute.pdf}{impute} package. }
\value{ A complete expression data matrix with missing values imputed. }
\references{Hastie, T., Tibshirani, R., Sherlock, G., Eisen, M., Brown, P. and Botstein, D., Imputing Missing Data for Gene Expression Arrays, Stanford University Statistics Department Technical report (1999), http://www-stat.stanford.edu/~hastie/Papers/missing.pdf

Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein and Russ B. Altman, Missing value estimation methods for DNA microarrays BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525
}
\author{ Cosmin Lazar }
\seealso{\code{\link{impute.knn}}, \code{\link{impute.wrapper.SVD}}, \code{\link{impute.wrapper.MLE}}}
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
  
exprsData.imputed = impute.wrapper.KNN(exprsData.MD,15)

## The function is currently defined as
function (dataSet.mvs, K) 
{
    resultKNN = impute.knn(dataSet.mvs, k = K, rowmax = 0.99, 
        colmax = 0.99, maxp = 1500, rng.seed = sample(1:1000, 
            1))
    dataSet.imputed = resultKNN[[1]]
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
