\name{impute.MinProb}
\alias{impute.MinProb}
\title{ Imputation of left-censored missing data using stochastic minimal value approach. }
\description{ Performs the imputation of left-censored missing data by random draws from a Gaussian distribution centered in a minimal value. Considering a peptide/protein expression data matrix with \code{n} columns corresponding to biological samples and \code{p} lines corresponding to peptides/proteins, for each sample (column), the mean value of the Gaussian distribution is set to a minimal value observed in that sample. The minimal value observed is estimated as being the \code{q-th} quantile (e.g. \code{q = 0.01}) of the observed values in that sample. The standard deviation is estimated as the median of the peptide/protein-wise standard deviations. Note that when estimating the standard deviation of the Gaussian distribution, only the peptides/proteins which present more than \code{50\%} recorded values are considered. }
\usage{ impute.MinProb(dataSet.mvs, q = 0.01, tune.sigma = 1) }
\arguments{
  \item{dataSet.mvs}{ A data matrix containing left-censored missing data. }
  \item{q}{ A scalar used to determine a low expression value to be used for missing data imputation. \code{0 < q < 1}, in this case q should be set to a low value. The default value is \code{q = 0.01}. }
  \item{tune.sigma}{ A scalar used to control the standard deviation of the Gaussian distribution used for random draws. If the sd is overestimated, than \code{0 < sigma.coef < 1}. The default value is \code{tune.sigma = 1}. }
}
\value{ A complete expression data matrix with missing values imputed. }
\author{ Cosmin Lazar }
\seealso{\code{\link{impute.QRILC}}, \code{\link{impute.MinDet}}, \code{\link{impute.ZERO}}}
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
MNAR.rate = 50
exprsData.MD.obj = insertMVs(exprsData,m.THR,sd.THR,MNAR.rate)
exprsData.MD = exprsData.MD.obj[[2]]

# perform missing data imputation
exprsData.imputed = impute.MinProb(exprsData.MD,0.01,1)

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (dataSet.mvs, q = 0.01, tune.sigma = 1) 
{
    nSamples = dim(dataSet.mvs)[2]
    nFeatures = dim(dataSet.mvs)[1]
    dataSet.imputed = dataSet.mvs
    min.samples = apply(dataSet.imputed, 2, quantile, prob = q, 
        na.rm = T)
    count.NAs = apply(!is.na(dataSet.mvs), 1, sum)
    count.NAs = count.NAs/nSamples
    dataSet.filtered = dataSet.mvs[which(count.NAs > 0.5), ]
    protSD = apply(dataSet.filtered, 1, sd)
    sd.temp = median(protSD, na.rm = T) * tune.sigma
    print(sd.temp)
    for (i in 1:(nSamples)) {
        dataSet.to.impute.temp = rnorm(nFeatures, 
                                        mean = min.samples[i], 
                                        sd = sd.temp)
        dataSet.imputed[which(is.na(dataSet.mvs[, i])), i] = 
        dataSet.to.impute.temp[which(is.na(dataSet.mvs[,i]))]
    }
    return(dataSet.imputed)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
