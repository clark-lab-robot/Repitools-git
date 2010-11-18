\name{enrichmentPlot}
\alias{enrichmentPlot}
\alias{enrichmentPlot,GenomeDataList,BSgenome-method}
\alias{enrichmentPlot,GRangesList,BSgenome-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the distribution of sequencing enrichment.}
\description{Function to generate a plot of the distribution of sequencing reads enrichments.}
\usage{
  \S4method{enrichmentPlot}{GenomeDataList,BSgenome}(rs, organism, seqLen, cols=rainbow(length(rs)), xlim=c(0,20), main="Enrichment Plot", total.lib.size=TRUE, verbose=FALSE, ...)
  \S4method{enrichmentPlot}{GRangesList,BSgenome}(rs, organism, seqLen, cols=rainbow(length(rs)), xlim=c(0,20), main="Enrichment Plot", total.lib.size=TRUE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GRangesList} or \code{GenomeDataList} object of reads to plot enrichment of.}
  \item{organism}{The \code{BSgenome} object to calculate enrichment upon - only used to extract chromosome lengths.}
  \item{seqLen}{The fragment size to be used for extending the sequencing reads.}
  \item{cols}{The line colour for each element of \code{rs}}
  \item{xlim}{\code{xlim} parameter passed to \code{plot}, the default is appropriate for \code{linear} \code{cpgDensityCalc} weighting.}
  \item{main}{\code{main} parameter passed to \code{plot}}
  \item{total.lib.size}{Whether to normalise enrichment values to the total number of reads per lane.}
  \item{verbose}{Print details of processing.}
  \item{...}{Additional graphical parameters to pass to \code{plot}.}
}
\details{
  See \code{enrichmentCalc} for details of how the results are determined.
}
\value{
  A plot is created. The data processed by \code{enrichmentCalc} is invisibly returned. 
}
\author{Aaron Statham}
\examples{
#See examples in the manual.
}
