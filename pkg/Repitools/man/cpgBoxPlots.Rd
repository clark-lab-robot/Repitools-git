\name{cpgBoxplots}
\alias{cpgBoxplots}
\alias{cpgBoxplots.AffymetrixCelSet}
\alias{cpgBoxplots.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplots of intensity, binned by Cpg Density }
\description{Either makes a side by side boxplot of two designs, or plots a single boxplot for the difference between the two designs.
}
\usage{
\method{cpgBoxplots}{AffymetrixCelSet}(this, ..., samples=c(1,2), subsetChrs="chr[1-5]", gcContent=7:18, 
                                                     calcDiff=FALSE, verbose=FALSE, nBins=40, pdfFile=NULL,
													 ylim=if (calcDiff) c(-5,6) else c(4,15), 
													 col=if (calcDiff) "salmon" else c("lightgreen","lightblue"),
													 mfrow=if (!is.null(pdfFile)) c(2,2) else c(1,1))
\method{cpgBoxplots}{matrix}(this, ndfTable, ..., samples=c(1,2), subsetChrs="chr[1-5]", gcContent=7:18, 
                                                     calcDiff=FALSE, verbose=FALSE, nBins=40, pdfFile=NULL,
													 ylim=if (calcDiff) c(-5,6) else c(4,15), 
													 col=if (calcDiff) "salmon" else c("lightgreen","lightblue"),
													 mfrow=if (!is.null(pdfFile)) c(2,2) else c(1,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{Either an AffymetrixCelSet or a matrix of intensity data.}
  \item{ndfTable}{In the case of Nimblegen data, a \code{data.frame} with at least columns \code{chr} and \code{sequence}. Must be in the same order of rows as the intensity data.}
  \item{...}{Other arguments for graphing}
  \item{samples}{Which 2 columns from the data matrix to use.}
  \item{subsetChrs}{Which chromosomes to limit the analysis to.}
  \item{gcContent}{A range of GC content, which only probes that have GC content in the range are used for the graphing.}
  \item{calcDiff}{Boolean. Plot the difference between the two samples ?}
  \item{verbose}{Boolean. Print processing output.}
  \item{nBins}{Bins to bin the intensities into.}
  \item{pdfFile}{Name of file to output plots to.}
  \item{ylim}{Y limit of graphs}
  \item{col}{Colour of boxes.}
  \item{mfrow}{Not specified by the user. Rows and columns to draw the plots in.}
}
\details{
  CpG content of probes is calculated in a 600 base window surrounding the probe, with a linearly decresasing weighting further away from the probe. 
}
\value{
 Invisibly returns a list of the plots.
}
\author{Mark Robinson, Dario Strbenac}