\name{regionStats} 
\Rdversion{1.1}
\alias{regionStats}
\alias{regionStats.AffymetrixCelSet}
\alias{regionStats.default}
\title{Find Regions of significance in microarray data}
\description{
The function finds the highest smoothed score cutoff for a pre-specified FDR. Smoothing is performed over a specified number of basepairs, and regions must have a minimum number of qualifying probes to be considered significant. The FDR is calculated as the ratio of the number of significant regions found in a permutation-based test, to the number found in the actual experimental microarray data. 
}
\usage{
\method{regionStats}{default}(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ndf, returnTrimmedMeans = FALSE, ...)
\method{regionStats}{AffymetrixCelSet}(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ind=NULL, returnTrimmedMeans = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{A set of raw microarray intensity data files}
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
	A list with elements
	\item{regions}{A list of \code{data.frame}. Each \code{data.frame} has columns \code{chr}, \code{start}, \code{end}, \code{score}.}
	\item{tMeanReal}{Matrix of smoothed scores of intensity data. Each column is an experimental design.}
	\item{tMeanPerms}{Matrix of smoothed scores of permuted intensity data. Each column is an experimental design.}
	\item{fdrTables}{List of table of FDR at different score cutoffs. Each list element is for a different experimental design.}
}
\author{Mark Robinson}
\examples{
function (cs, design, fdrLevel = 0.05, nPermutations = 5, probeWindow = 600, meanTrim = 0.1, nProbes = 10, verbose = TRUE, ...) 
UseMethod("regionStats")
}
