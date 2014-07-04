\name{generate.ExpressionData}
\alias{generate.ExpressionData}
\title{ Generate Peptide/Protein Expression Data }
\description{ Generates a matrix of peptide/protein expression data. It is assumed that the expression data is log - transformed. Therefore, for each sample the peptides/proteins intensities are randomly drawn from a Gaussian distribution.}
\usage{
generate.ExpressionData(nSamples1, nSamples2, meanSamples, sdSamples, 
                        nFeatures, nFeaturesUp, nFeaturesDown, 
                        meanDynRange, sdDynRange, 
                        meanDiffAbund, sdDiffAbund)
}
\arguments{
  \item{nSamples1}{ Number of samples in condition 1. }
  \item{nSamples2}{ Number of samples in condition 2. }
  \item{meanSamples}{ Mean value of the background noise. }
  \item{sdSamples}{ Standard deviation of the background noise. }
  \item{nFeatures}{ Number of peptides/proteins. }
  \item{nFeaturesUp}{ Number of peptides/proteins up-regulated. }
  \item{nFeaturesDown}{ Number of /peptidesproteins down-regulated. }
  \item{meanDynRange}{ Mean value of the dynamic range of peptide/protein expressions. }
  \item{sdDynRange}{ Standard deviation of the dynamic range of peptide/protein expressions. }
  \item{meanDiffAbund}{ Mean value of the up/down-regulation. }
  \item{sdDiffAbund}{ Standard deviation of the up/down-regulation. }
}
\value{A list including elements:

\item{data }{Peptide/protein expression matrix}
\item{conditions }{Vector indicating the samples in each condition}
\item{labelFeatures }{Vector indicating features which are down/up/no regulated}

}
\author{ Cosmin Lazar }
\seealso{\code{\link{pep2prot}}, \code{\link{generate.RollUpMap}}, \code{\link{insertMVs}}}
\examples{

dataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                          meanSamples = 0, sdSamples = 0.2,
                          nFeatures = 2000, nFeaturesUp = 100, nFeaturesDown = 100,
                          meanDynRange = 20, sdDynRange = 1,
                          meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = dataObj[[1]]
  
\dontrun{
hist(exprsData[,1])
}

## The function is currently defined as
function (nSamples1, nSamples2, meanSamples, sdSamples, nFeatures, 
    nFeaturesUp, nFeaturesDown, meanDynRange, sdDynRange, meanDiffAbund, 
    sdDiffAbund) 
{
    nSamples = nSamples1 + nSamples2
    data = matrix(rnorm(nSamples * nFeatures, meanSamples, sdSamples), 
        nFeatures, nSamples)
    means = rnorm(nFeatures, meanDynRange, sdDynRange)
    data = data + means
    conditions = c(rep(1, nSamples1), rep(2, nSamples2))
    DE.coef.up = matrix(rnorm(nFeaturesUp * nSamples1, meanDiffAbund, 
        sdDiffAbund), nFeaturesUp, nSamples1)
    DE.coef.down = matrix(rnorm(nFeaturesDown * nSamples2, meanDiffAbund, 
        sdDiffAbund), nFeaturesDown, nSamples2)
    data[1:nFeaturesUp, conditions == 1] = DE.coef.up + data[1:nFeaturesUp, 
        conditions == 1]
    data[(nFeaturesUp + 1):(nFeaturesUp + nFeaturesDown), conditions == 
        2] = DE.coef.down + data[(nFeaturesUp + 1):(nFeaturesUp + 
        nFeaturesDown), conditions == 2]
    labelFeatures = c(rep(1, nFeaturesUp), rep(2, nFeaturesDown), 
        rep(3, nFeatures - (nFeaturesUp + nFeaturesDown)))
    row.names(data) = 1:nFeatures
    return(list(data, conditions, labelFeatures))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
