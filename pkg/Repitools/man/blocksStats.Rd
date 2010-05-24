\name{blocksStats}
\alias{blocksStats}
\alias{blocksStats.AffymetrixCelSet}
\alias{blocksStats.GenomeDataList}
\alias{blocksStats.matrix}
\title{A routine to calculate statistics for regions in the genome}
\description{
For each region of interest or TSS, this routine interrogates probes or sequence data for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.  Regions can be surroundings of TSSs, or can be a table of user - specified regions. The function determines if the \code{start} and \code{end} columns of the \code{coordinatesTable} should be used as regions or as TSSs, based on the \code{useAsRegions} flag.
}
\usage{
\method{blocksStats}{AffymetrixCelSet}(cs, coordinatesTable, annot = NULL, probePositions = NULL, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE, useAsRegions=FALSE, ...)
\method{blocksStats}{GenomeDataList}(cs, coordinatesTable, design, upStream = 0, downStream = 2000, verbose = TRUE, useAsRegions = FALSE, seqLen = NULL, libSize = "lane", Acutoff = NULL, ...)
\method{blocksStats}{matrix}(cs, ndf, coordinatesTable, annot = NULL, probePositions = NULL, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE, useAsRegions=FALSE, ...)
} 
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, or a data.frame with columns called, presumably from normalized data}
  \item{ndf}{A \code{data.frame} with columns \code{chr}, \code{position}, and \code{index}}
  \item{coordinatesTable}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{name}, \code{start}, \code{end}, \code{strand}.}
  \item{annot}{If a mapping with \code{annotationLookup} or \code{annotationBlocksLookup} has already been done, it can be passed in, and avoids re-conmputing of the mapping list within \code{blocksStats}.}
   \item{probePositions}{Only include if specifying the \code{annot} argument. A \code{data.frame} with columns called \code{chr} and \code{position}.}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{upStream}{number of bases upstream to consider in calculation of statistics}
  \item{downStream}{number of bases upstream to consider in calculation of statistics}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{useAsRegions}{If true, the \code{start} and \code{end} columns of \code{coordinatesTable} will define the region to calculate statistics in.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{libSize}{A string that indicates whether to use the total lane count, total count within regions specified by \code{coordinatesTable}, or normalisation to a reference lane by the negative binomial quantile - to - quantile method, as the library size for each lane. For total lane count use \code{"lane"}, for region sums use \code{"inRegions"}, and for the normalisation use \code{"ref"}.}
  \item{Acutoff}{If \code{libSize} is \code{"ref"}, this argument must be provided. Otherwise, it must not. Cutoff on the "A" values to take, before calculating trimmed mean.}
  \item{...}{Arguments passed to internal function for statistics calculations}
}
\details{
The approach simply calculates the contrast for each region and tests whether the collection of these summaries are significantly non-zero.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{coordinatesTable}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}
