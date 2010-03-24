\name{getRMAdata}
\alias{getRMAdata}
\title{Get RMA Data}
\description{A convenient and quick way to get RMA processed data.}
\usage{
getRMAdata(name, chipType, units = NULL, doCore = FALSE, doNorm = FALSE, doExon = FALSE, doMerge = FALSE, doResiduals = FALSE, doWeights = FALSE, force = FALSE)
}
\arguments{
  \item{name}{The name of the directory within the rawData directory that the CEL files are in.}
  \item{chipType}{A character string giving the type of the chip.}
  \item{units}{Which units to fit a model to, and extract.}
  \item{doCore}{Whether to do the core set of probesets, or all of them.}
  \item{doNorm}{Whether to do background correction and quantile normalisation on the CEL set or not.}
  \item{doExon}{Whether to use the exon - based implementation of the RMA algorithm or not.}
  \item{doMerge}{Whether to merge the exons into gene level summaries (if doExon was also TRUE) or not.}
  \item{doResiduals}{Whether to calculate residuals or not.}
  \item{doWeights}{Whether to calculate wieghts or not.}
  \item{force}{Write over the previous cached analysis, if one already exists.}
  
}
\details{
 	This is a convenience function that does all the steps in a traditional RMA pipeline.
 }
\value{
	A list with elements
	\item{name}{Name of the CEL set.}
	\item{chipType}{Chip type of the CEL set.}
	\item{exprs}{A \code{data.frame} with RMA expression data.}
	\item{cdf}{CDF file used.}
	\item{plm}{Probe linear model.}
	\item{qam}{Quality assessment model.}
	\item{residuals}{Residuals of the PLM.}
	\item{weights}{Weights of the PLM.}
}
\author{Mark Robinson, Dario Strbenac} 


