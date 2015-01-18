\name{impute.MAR.MNAR}
\alias{impute.MAR.MNAR}
\title{ Hybrid imputation method }
\description{
Performs the imputation of missing values assuming that the missing values are both MAR and MNAR. It is assumed that the MNAR missing values are left-censored, more precisely, only rows (proteins/peptides) with a mean value below a censoring threshold are considered to contain left-censored missing data. The method relies on a estimation of the left-censoring threshold that is further used to distinguish rows (proteins/peptides) that contain left-censored missing data from those lines who contain random missing data.}
\usage{
impute.MAR.MNAR(dataSet.mvs, model.selector, method.MAR = "KNN", method.MNAR = "QRILC")
}
\arguments{
  \item{dataSet.mvs}{ A data matrix containing missing values. }
  \item{model.selector}{ Vector containing the flag indicators allowing the selection of the appropriate method for the imputation of missing data: the rows (peptides/proteins) corresponding to "1" in the model.selector vector are treated using a MAR specific method while the ones corresponding to "0" are treated using a MNAR specific method.}
  \item{method.MAR}{ The method employed for the imputation of MAR/MCAR missing data. Possible values:  }
  \item{method.MNAR}{ The method employed for the imputation of the left-censored missing data. }
}
\value{ Data matrix containing complete abundances }
\author{ Cosmin Lazar }
\note{ The method makes the assumption that rows (peptides/proteins) with a mean value lower than the censoring threshold are not affected my a random missingness mechanism. In reality, a random missingness mechanism could affect also rows (peptides/proteins) which are below the censoring threshold. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{model.Selector}}, \code{\link{impute.MAR}} }
\examples{
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

# run model.Selector
m.s = model.Selector(exprsData.MD)

# perform MAR-MNAR imputation
exprsData.MAR.imputed = impute.MAR.MNAR(exprsData.MD, m.s, 
                           method.MAR = "KNN", method.MNAR = "QRILC")

## The function is currently defined as
function (dataSet.mvs, model.selector, method.MAR = "KNN", method.MNAR = "QRILC") 
{
    switch(method.MAR, MLE = {
        dataSet.MCAR.imputed = impute.MAR(dataSet.mvs, model.selector, 
            method = "MLE")
    }, SVD = {
        dataSet.MCAR.imputed = impute.MAR(dataSet.mvs, model.selector, 
            method = "SVD")
    }, KNN = {
        dataSet.MCAR.imputed = impute.MAR(dataSet.mvs, model.selector, 
            method = "KNN")
    })
    switch(method.MNAR, QRILC = {
        dataSet.complete.obj = impute.QRILC(dataSet.MCAR.imputed, 
            tune.sigma = 0.3)
        dataSet.complete = dataSet.complete.obj[[1]]
    }, MinDet = {
        dataSet.complete = impute.MinDet(dataSet.MCAR.imputed)
    }, MinProb = {
        dataSet.complete = impute.MinProb(dataSet.MCAR.imputed)
    })
    return(dataSet.complete)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
