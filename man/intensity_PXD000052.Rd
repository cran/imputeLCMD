\name{intensity_PXD000052}
\alias{intensity_PXD000052}
\docType{data}
\title{ Dataset PXD000052 from ProteomeXchange. }
\description{ This dataset has been collected during a study designed to perform the proteomic analysis of the SLP76 interactome in resting and activated primary mast cells. Four SLP76 replicates (with two analytical replicates each) have been affinity-purified from both resting and activated primary mast cells. The dataset within this package contains proteins intensity processed using MaxQuant. More information can be found on ProteomeExchange public repository (http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000052) or in the original paper (see reference). }
\usage{data(intensity_PXD000052)}
\format{
  A data frame with 1991 observations on the following 17 variables.
  \describe{
    \item{\code{Protein.IDs}}{Peptides/Proteins names}
    \item{\code{iBAQ.stSLP_activ1}}{a numeric vector}
    \item{\code{iBAQ.stSLP_activ2}}{a numeric vector}
    \item{\code{iBAQ.stSLP_activ3}}{a numeric vector}
    \item{\code{iBAQ.stSLP_activ4}}{a numeric vector}
    \item{\code{iBAQ.stSLP_rest1}}{a numeric vector}
    \item{\code{iBAQ.stSLP_rest2}}{a numeric vector}
    \item{\code{iBAQ.stSLP_rest3}}{a numeric vector}
    \item{\code{iBAQ.stSLP_rest4}}{a numeric vector}
    \item{\code{iBAQ.WT_activ1}}{a numeric vector}
    \item{\code{iBAQ.WT_activ2}}{a numeric vector}
    \item{\code{iBAQ.WT_activ3}}{a numeric vector}
    \item{\code{iBAQ.WT_activ4}}{a numeric vector}
    \item{\code{iBAQ.WT_rest1}}{a numeric vector}
    \item{\code{iBAQ.WT_rest2}}{a numeric vector}
    \item{\code{iBAQ.WT_rest3}}{a numeric vector}
    \item{\code{iBAQ.WT_rest4}}{a numeric vector}
  }
}
\source{ Original MaxQuant data: http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000052 }
\references{ Bounab Y, Hesse AM, Iannascoli B, Grieco L, Coute Y, Niarakis A, Roncagalli R, Lie E, Lam KP, Demangel C, Thieffry D, Garin J, Malissen B, Da?ron M, Proteomic analysis of the SH2 domain-containing leukocyte protein of 76 kDa (SLP76) interactome in resting and activated primary mast cells [corrected]. Mol Cell Proteomics, 12(10):2874-89(2013). }
\examples{
data(intensity_PXD000052)
}
\keyword{datasets}
