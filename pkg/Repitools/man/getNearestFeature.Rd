\name{getNearestFeature}
\alias{getNearestFeature}
\title{Find the closest features of regions}
\description{
	Given a table of regions of interest, and a table of features (e.g. a table of genes), the regions of interest table is returned, with columns that describe the closest feature appended.
}
\usage{
getNearestFeature(peaksTable, featuresTable)
}
\arguments{
  \item{peaksTable}{A data frame with (at least) \code{chr}, \code{start}, and \code{end} columns.}
  \item{featuresTable}{A data frame with (at least) \code{chr}, \code{start}, and \code{end} columns.}
}
\details{
  
}
\value{
A \code{data.frame} with columns \code{feature}, \code{featureStart}, \code{featureEnd} \code{peakLocation}, and \code{distance}. The data frame is sorted to be in the same order as the \code{peaksTable} 
}
\author{Dario Strbenac}
\examples{
	require(ChIPpeakAnno)
	peaksT <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"), start = c(2000010, 19000000, 30000000, 300, 5500, 100000), end = c(2000310, 19000300, 30000300, 600, 5800, 100300))
	featuresT <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"), start = c(1000000, 10000000, 15000000, 1000, 6000, 10000), end = c(2000000, 20000000, 22000000, 5000, 7000, 15000), strand = c('+', '-', '+', '+', '-', '+'))
	getNearestFeature(peaksT, featuresT)
}
