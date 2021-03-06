% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modifyLayerByBindingFactor.BSgenome.R
\name{modifyLayerByBindingFactor.BSgenome}
\alias{modifyLayerByBindingFactor.BSgenome}
\title{modifyLayerByBindingFactor.BSgenome}
\usage{
modifyLayerByBindingFactor.BSgenome(layerSet, hits, bindingFactor,
  verbose = FALSE)
}
\arguments{
\item{layerSet}{A layerset object}

\item{hits}{The set of hits (GRanges) determined by \code{\link{matchBindingFactor}} or \code{\link{matchBindingFactor.BSgenome}}}

\item{bindingFactor}{An individual bindingFactor object}

\item{verbose}{Give more output}
}
\value{
\code{"LayerList"}
}
\description{
Modify all hits genome-wide on a layerSet built upon a BSgenome object.
}
\examples{
x <- 1   # great!

}
\seealso{
\code{\link{runLayerBinding.BSgenome}}  \code{\link{matchBindingFactor}}  \code{\link{matchBindingFactor.BSgenome}}
}
