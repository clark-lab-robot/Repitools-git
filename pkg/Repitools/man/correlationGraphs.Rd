\name{correlationGraphs}
\alias{correlationGraphs}
\title{Creates graphs of the correlation between ChIP-chip data and expression data}
\description{Creates graphs of the correlation between ChIP-chip data and expression data}
\usage{
correlationGraphs(intensities, lookup, compareExpression, ...)
}
\arguments{
  \item{intensities}{A matrix of array data, one column per sample, one row per probe}
  \item{lookup}{A mapping between tiling array probes and expression}
  \item{compareExpression}{A \code{vector} of expression values to correlate the array intensities with}
  \item{...}{Parameters passed into \code{plot}}
}
\value{Graphs are plotted}
\author{Aaron Statham}
\seealso{\code{\link{significancePlots}}}
\examples{
#See examples in the manual
}
