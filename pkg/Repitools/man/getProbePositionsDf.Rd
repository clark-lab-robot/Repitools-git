\name{getProbePositionsDf}
\alias{getProbePositionsDf}
\alias{getProbePositionsDf.AffymetrixCdfFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Translate probe information in a table}
\description{Translates the probe information in the AromaCellPositionFile to a data.frame object}

\usage{
\method{getProbePositionsDf}{AffymetrixCdfFile}(cdf, ..., verbose = -20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{an AffymetrixCdfFile object}
  \item{\dots}{further arguments to send to \code{getCellIndices}}
  \item{verbose}{logical, whether or not to print out progress statements to the screen}
}
\details{
This assumes that the AromaCellPositionFile exist.
}
\value{
A data.frame with 3 columns: chr, position, index
}
\author{Mark Robinson}
\examples{
## not run
# probePositions <- getProbePositionsDf(cdfU)
}
