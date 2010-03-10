\name{processNDF}
\alias{processNDF}
\title{Reads in a Nimblegen microarray design file (NDF)}
\description{
Reads a Nimblegen microarray design file (NDF file) which describes positions and sequences of probes on a Nimblegen microarray.
}
\usage{
processNDF(filename, ncols = 768)
}
\arguments{
  \item{filename}{the name of the Nimblegen microarray design file}
  \item{ncols}{the number of columns of probes on the array - must be the same value as will be passed to \code{loadPairFile} or \code{loadSampleDirectory}. The default works for 385K format arrays.}
}
\details{
Reads in a Nimblegen microarray design file. This enables the reading in and annotation of Nimblegen microarray data files (pair files).
}
\value{
a data frame containing 
\item{chr}{the chromosome the probe was designed against}
\item{position}{the position of the sequence the probe was designed against (probe centre)}
\item{strand}{the strand the probe was designed against}
\item{index}{the index (x y position) the probe occupies on the array}
\item{sequence}{the actual DNA sequence synthesised onto the array}
\item{GC}{the percent GC content of the probe sequence}
}
\author{Aaron Statham}
\seealso{\code{\link{loadSampleDirectory}}, \code{\link{loadPairFile}}}
\examples{
# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]
}

