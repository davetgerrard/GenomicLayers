% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createBindingFactor.DNA_motif.R
\name{createBindingFactor.DNA_motif}
\alias{createBindingFactor.DNA_motif}
\title{Create a binding factor object to match a given DNA motif}
\usage{
createBindingFactor.DNA_motif(
  name,
  type = "DNA_motif",
  pwm,
  min.score = "80\%",
  with.score = FALSE,
  patternLength = ncol(pwm),
  stateWidth = patternLength,
  profile.layers = NULL,
  profile.marks = NULL,
  mod.layers = NULL,
  mod.marks = NULL,
  offset = 0,
  offset.method = NULL,
  offset.params = NULL,
  test.layer0.binding = FALSE,
  test.mismatch.rate = 0.1,
  max.pattern.tries = 1000,
  verbose = FALSE
)
}
\arguments{
\item{name}{give the binding factor a name}

\item{type}{"DNA_motif"  to differentiate from other types}

\item{pwm}{NULL put motif here as matrix (see biostrings \code{\link{pwm}})}

\item{min.score}{"80\%"   passed to  \code{\link{matchPWM}}}

\item{with.score}{FALSE   passed to  \code{\link{matchPWM}}}

\item{patternLength}{length of pattern to be matched [ncol(pwm)]}

\item{stateWidth}{the width of pattern to recognise on other layers}

\item{profile.layers}{a vector of named layers to set as a match}

\item{profile.marks}{a vector of 0/1 to match the layers in profile.layers}

\item{mod.layers}{a vector of named layers to alter on a match}

\item{mod.marks}{a vector of 0/1 to set on the mod.layers}

\item{offset}{0  integer value to indicate relative distance from pattern to apply modifications. Very simple.}

\item{offset.method}{NULL   a \code{function} to apply to apply offset. MUST have parameter "n" that is used internally to represent the number of hits.}

\item{offset.params}{NULL  a \code{list} of named parameters to pass to offset.method function}

\item{test.layer0.binding}{when creating, test if the DNA sequence has a match.}

\item{test.mismatch.rate}{proportion of mismatches to tolerate when testing [.1]}

\item{verbose}{set to TRUE for more output}
}
\value{
\code{"bindingFactor"}
}
\description{
Create a new binding factor based on a DNA motif that \emph{may}  also require 
marks on others layers and \emph{may} (when used) set marks on other layers.
Makes use of biostrings function \code{\link{vmatchPattern}}
Can use IUPAC codes and allow mismatches
}
\examples{
require(Biostrings)
data(HNF4alpha)    # from Biostrings package
pwm.HNF4A <- PWM(HNF4alpha)  
bf_motif.1 <- createBindingFactor.DNA_motif(name="HNFA_match",pwm = pwm.HNF4A )

bf_motif.2 <- createBindingFactor.DNA_motif(name="HNFA_alter",pwm = pwm.HNF4A ,
                          min.score="80\%", with.score=FALSE,
                          profile.layers = NULL,
                          profile.marks=NULL,
                          mod.layers="LAYER.1", mod.marks = 1)

}
\seealso{
\code{\link{runLayerBinding}} \code{\link{createBindingFactor.DNA_regexp}}
}
