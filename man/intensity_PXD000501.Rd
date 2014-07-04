\name{intensity_PXD000501}
\alias{intensity_PXD000501}
\docType{data}
\title{ Dataset PXD000501 from ProteomeXchange. }
\description{ This dataset contains three biological replicates with three technical replicates each for the conditiones media (CM) and the whole cell lysates (WCL) of C8-D1A cell lines. The dataset within this package contains proteins iBAQ intensity processed using MaxQuant. More information can be found on ProteomeExchange public repository (http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000501) or in the original paper (see reference). }
\usage{data(intensity_PXD000501)}
\format{
  A data frame with 7363 observations on the following 19 variables.
  \describe{
    \item{\code{Protein.IDs}}{Peptides/Proteins names}
    \item{\code{iBAQ.secretome_set1_tech1}}{a numeric vector}
    \item{\code{iBAQ.secretome_set1_tech2}}{a numeric vector}
    \item{\code{iBAQ.secretome_set1_tech3}}{a numeric vector}
    \item{\code{iBAQ.secretome_set2_tech1}}{a numeric vector}
    \item{\code{iBAQ.secretome_set2_tech2}}{a numeric vector}
    \item{\code{iBAQ.secretome_set2_tech3}}{a numeric vector}
    \item{\code{iBAQ.secretome_set3_tech1}}{a numeric vector}
    \item{\code{iBAQ.secretome_set3_tech2}}{a numeric vector}
    \item{\code{iBAQ.secretome_set3_tech3}}{a numeric vector}
    \item{\code{iBAQ.whole_set1_tech1}}{a numeric vector}
    \item{\code{iBAQ.whole_set1_tech2}}{a numeric vector}
    \item{\code{iBAQ.whole_set1_tech3}}{a numeric vector}
    \item{\code{iBAQ.whole_set2_tech1}}{a numeric vector}
    \item{\code{iBAQ.whole_set2_tech2}}{a numeric vector}
    \item{\code{iBAQ.whole_set2_tech3}}{a numeric vector}
    \item{\code{iBAQ.whole_set3_tech1}}{a numeric vector}
    \item{\code{iBAQ.whole_set3_tech2}}{a numeric vector}
    \item{\code{iBAQ.whole_set3_tech3}}{a numeric vector}
  }
}
\source{ Original MaxQuant data: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000501 }
\references{ Han D, Jin J, Woo J, Min H, Kim Y, Proteomic analysis of mouse astrocytes and their secretome by a combination of FASP and StageTip-based, high pH, reversed-phase fractionation. Proteomics, ():(2014). }
\examples{
data(intensity_PXD000501)
}
\keyword{datasets}
