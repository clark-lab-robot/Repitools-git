\name{cpgDensityPlot}
\alias{cpgDensityPlot}
\alias{cpgDensityPlot.GenomeDataList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the distribution of sequencing reads CpG densities.}
\description{Function to generate a plot of the distribution of sequencing reads CpG densities.}
\usage{
\method{cpgDensityPlot}{GenomeDataList}(rs, cols=rainbow(length(rs)), xlim=c(0,20), main="CpG Density Plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeDataList} object of reads to plot CpG density of}
  \item{cols}{The line colour for each element of \code{rs}}
  \item{xlim}{\code{xlim} parameter passed to \code{plot}, the default is appropriate for \code{linear} \code{cpgDensityCalc} weighting.}
  \item{main}{\code{main} parameter passed to \code{plot}}
  \item{...}{Arguments passed into \code{cpgDensityCalc}, \code{seqLen} and \code{organism} are required}
}
\details{
  See \code{cpgDensityCalc} for details of options for calculating the CpG density.
}
\value{
  A plot is created. The data processed by \code{cpgDensityCalc} is invisibly returned. 
}
\author{Aaron Statham}
\examples{
#See the manual.
}
