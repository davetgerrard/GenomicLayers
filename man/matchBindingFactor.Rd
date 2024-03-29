% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchBindingFactor.R
\name{matchBindingFactor}
\alias{matchBindingFactor}
\title{Find matches for a binding factor on a layer set}
\usage{
matchBindingFactor(
  layerSet,
  bindingFactor,
  clusterGap = 10,
  max.window = 1e+07,
  verbose = FALSE
)
}
\arguments{
\item{layerSet}{the \code{"layerSet"} target}

\item{bindingFactor}{the \code{"bindingFactor"} to match}

\item{clusterGap}{=10 NOT IMPLEMENTED}

\item{max.window}{=10000000 on less powerful computers, break up the search into windows of this size.}

\item{verbose}{output more information to the screen}
}
\value{
\code{"hits"}
}
\description{
Generate a list of matches for a binding factor against a layerSet object.
}
\examples{
x <- 1   # great!

}
