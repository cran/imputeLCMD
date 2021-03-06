\name{insertMVs}
\alias{insertMVs}
\title{ Generates missing data in a complete data matrix. }
\description{ This function generates missing data in a complete data matrix. Both random and left-censored missing data can be generated. The percentage of all missing data is controlled by \code{mean.THR}. The percentage of missing data which are left-censored is controlled by \code{MNAR.rate}. }
\usage{ insertMVs(original, mean.THR, sd.THR, MNAR.rate) }
\arguments{
  \item{original}{ Original complete data matrix of peptide/protein expression. }
  \item{mean.THR}{ Mean value of the threshold distribution which controls the total missing data rate (\code{mean.THR} should be initially set such that the result of the initial thresholding, in terms of no. of \code{NAs}, equals the total missing values rate). Example: if one wants to generate \code{30\%} missing data, \code{mean.THR} can be set as follows: \code{mean.THR = quantile(exprsData, probs = 0.3)}.}
  \item{sd.THR}{ Standard deviation of the threshold distribution which controls the total missing data rate. \code{sd.THR} is usually set to a small value (e.g. \code{0.1}). }
  \item{MNAR.rate}{ Percentage of MVs which are missing not at random. Among the total number of missing data (\code{NA_rate_TOTAL}) generated by the initial threshold, a percentage of missing data that equals \code{MNAR.rate} is preserved as such; the remaining missing data are replaced with the original values. Next, a proportion of (\code{NA_rate_TOTAL - MNAR.rate}) missing data are generated randomly.
}
}
\value{ A list including elements:
  \item{original}{Original complete data matrix}
  \item{original.mvs}{Data matrix derived from the original by generating missing data }
  \item{pNaNs}{The percetage of missing data generated in the original complete dataset }
}
\author{ Cosmin Lazar }
\seealso{\code{\link{generate.ExpressionData}}}
\examples{

# generate expression data matrix
exprsDataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                          meanSamples = 0, sdSamples = 0.2,
                          nFeatures = 1000, nFeaturesUp = 50, nFeaturesDown = 50,
                          meanDynRange = 20, sdDynRange = 1,
                          meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = exprsDataObj[[1]]
  
# insert 15\% missing data with 50\% missing not at random

m.THR = quantile(exprsData, probs = 0.15)
sd.THR = 0.1
MNAR.rate = 50
exprsData.MD.obj = insertMVs(exprsData,m.THR,sd.THR,MNAR.rate)
exprsData.MD = exprsData.MD.obj[[2]]

\dontrun{
hist(exprsData[,1])
hist(exprsData.MD[,1])
hist(exprsData.imputed[,1])
}

## The function is currently defined as
function (original, mean.THR, sd.THR, MNAR.rate) 
{
    originalNaNs = original
    nProt = nrow(original)
    nSamples = ncol(original)
    thr = matrix(rnorm(nSamples * nProt, mean.THR, sd.THR), nProt, 
        nSamples)
    indices.MNAR = which(original < thr)
    no.MNAR = round(MNAR.rate/100 * length(indices.MNAR))
    temp = matrix(original, 1, nSamples * nProt)
    temp[sample(indices.MNAR, no.MNAR)] = NaN
    indices.MCAR = which(!is.na(temp))
    no.MCAR = floor((100 - MNAR.rate)/100 * length(indices.MNAR))
    print(no.MCAR + no.MNAR)
    temp[sample(indices.MCAR, no.MCAR)] = NaN
    originalNaNs = matrix(temp, nProt, nSamples)
    originalNaNs_adjusted = originalNaNs
    noNaNs_Var = rowSums(is.na(originalNaNs))
    allNaNs_Vars = which(noNaNs_Var == nSamples)
    sampleIndexToReplace = sample(1:nSamples, length(allNaNs_Vars), 
        replace = T)
    for (i in 0:length(sampleIndexToReplace)) {
        originalNaNs_adjusted[allNaNs_Vars[i], sampleIndexToReplace[i]] = original[allNaNs_Vars[i], 
            sampleIndexToReplace[i]]
    }
    pNaNs = length(which(is.na(originalNaNs_adjusted)))/(nSamples * 
        nProt)
    print(pNaNs)
    return(list(original, originalNaNs_adjusted, pNaNs))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
