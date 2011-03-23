\name{regionStats} 
\Rdversion{1.1}
\alias{regionStats}
\alias{regionStats,AffymetrixCelSet-method}
\alias{show,RegionStats-method}
\alias{RegionStats-class}
\alias{class:RegionStats}
\alias{regionStats,matrix-method}
\title{Find Regions of significance in microarray data}
\description{
The function finds the highest smoothed score cutoff for a pre-specified FDR. Smoothing is performed over a specified number of basepairs, and regions must have a minimum number of qualifying probes to be considered significant. The FDR is calculated as the ratio of the number of significant regions found in a permutation-based test, to the number found in the actual experimental microarray data. 
}
\usage{
  \S4method{regionStats}{matrix}(x, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ndf, returnTrimmedMeans = FALSE, ...)
  \S4method{regionStats}{AffymetrixCelSet}(x, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ind=NULL, returnTrimmedMeans = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{AffymetrixCelSet} or \code{matrix} of array data to use.}
  \item{design}{A design matrix of how to manipulate }
  \item{fdrLevel}{Cutoff of the maximum acceptable FDR}
  \item{nPermutations}{Number of permutations to use }
  \item{probeWindow}{Size of window, in base pairs, to check for }
  \item{meanTrim}{A number representing the top and bottom fraction of ordered values in a window to be removed, before the window mean is calculated.}
  \item{nProbes}{Minimum number of probes in a window, for the region to qualify as a region of significance.}
  \item{maxGap}{Maximum gap between significant probes allowable.} 
  \item{twoSides}{Look for both significant positive and negative regions.}
  \item{verbose}{Print debugging comments while executing.} 
  \item{ind}{A vector of the positions of the probes on the array}
  \item{ndf}{The Nimblegen Definition File for Nimblegen array data.}
  \item{returnTrimmedMeans}{If TRUE, the values of the trimmed means of the intensities and permuted intensities are also retuned from the function.}
  \item{...}{Arguments passed to internal function for common processing.}
}
\value{
	A \code{RegionStats} object (list) with elements
	\item{regions}{A list of \code{data.frame}. Each \code{data.frame} has columns \code{chr}, \code{start}, \code{end}, \code{score}.}
	\item{tMeanReal}{Matrix of smoothed scores of intensity data. Each column is an experimental design.}
	\item{tMeanPerms}{Matrix of smoothed scores of permuted intensity data. Each column is an experimental design.}
	\item{fdrTables}{List of table of FDR at different score cutoffs. Each list element is for a different experimental design.}
}
\author{Mark Robinson}
\examples{
\dontrun{
library(Repitools)
library(aroma.affymetrix)

# assumes appropriate files are at annotationData/chipTypes/Hs_PromPR_v02/
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02",verbose=-20)
cdfU <- getUniqueCdf(cdf,verbose=-20)

# assumes appropriate files are at rawData/experiment/Hs_PromPR_v02/
cs <- AffymetrixCelSet$byName("experiment",cdf=cdf,verbose=-20)
mn <- MatNormalization(cs)
csMN <- process(mn,verbose=-50)
csMNU <- convertToUnique(csMN,verbose=-20)

#> getNames(cs)
# [1] "samp1"  "samp2"  "samp3"  "samp4"

design <- matrix( c(1,-1,rep(0,length(cs)-2)), ncol=1, dimnames=list(getNames(cs),"elut5_L-P") )

# just get indices of chr7 here
ind <- getCellIndices(cdfU, unit = indexOf(cdfU, "chr7F"), unlist = TRUE, useNames = FALSE)

regs <- regionStats(csMNU, design, ind = ind, probeWindow = 500, verbose = TRUE)
}
}
