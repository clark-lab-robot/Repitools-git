\name{loadPairFile}
\alias{loadPairFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to read Nimblegen tiling array intensities}
\description{
Reads a file in Nimblegen pair format, returning log2 intensities of probes referenced by the supplied ndf data frame.
}
\usage{
loadPairFile(filename, ndf, ncols = 768)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{the name of the pair file which intensities are to be read from.}
  \item{ndf}{a data frame produced by \code{\link{processNDF}}.}
  \item{ncols}{the number of columns of probes on the array - must be the same value as used in \code{processNDF}. The default works for 385K format arrays.}
}
\details{
Reads in intensities from the specified pair file, then matches probes against those specified in the supplied ndf.
}
\value{
  a \code{vector} of log2 intensities, the number of rows of the supplied \code{ndf} in length.
}
\author{Aaron Statham}
\seealso{\code{\link{loadSampleDirectory}} for reading multiple pair files with the same ndf. \code{\link{processNDF}}}
\examples{
# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]
#
## Read in a pair file using the chromosome only NDF
# arrayIntensity <- loadPairFile("Pairs/Array1_532.pair", ndf)
#
}
