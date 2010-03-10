\name{loadSampleDirectory}
\alias{loadSampleDirectory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to read Nimblegen tiling array intensities}
\description{
Reads all files in Nimblegen pair format within the specified directory, returning log2 intensities of probes referenced by the supplied ndf data frame.
}
\usage{
loadSampleDirectory(path, ndf, what="Cy3", ncols = 768)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{the directory containing the pair files to be read.}
  \item{ndf}{a data frame produced by \code{\link{processNDF}}.}
  \item{what}{specifies the channel(s) to be read in - either \code{Cy3}, \code{Cy5}, \code{Cy3/Cy5}, \code{Cy5/Cy3}, \code{Cy3andCy5}, \code{Cy5andCy3}.}
  \item{ncols}{the number of columns of probes on the array - must be the same value as used in \code{processNDF}. The default works for 385K format arrays.}
}
\details{
Reads in intensities of all arrays contained within \code{path}. The parameter \code{what} determines which fluorescent channels are read, and how the are returned.
\code{Cy3} and \code{Cy5} return the log2 intensity of the specified single channel.
\code{Cy3/Cy5} and \code{Cy5/Cy3} return the log2 ratio of the two channels.
\code{Cy3andCy5} and \code{Cy5andCy3} return the log2 intensity of both channels in separate columns of the matrix.
}
\value{
  a \code{matrix} of log2 intensites, with the same number of rows as the supplied \code{ndf} and depending on the value of \code{what} either one or two columns per array.
}
\author{Aaron Statham}
\seealso{\code{\link{loadPairFile}} for reading a single pair files. \code{\link{processNDF}}}
\examples{
# Not run
#
## Read in the NDF file 
# ndfAll <- processNDF("080310_HG18_chr7RSFS_AS_ChIP.ndf")
#
## Subset the NDF to only probes against chromosomes
# ndf <- ndfAll[grep("^chr", ndfAll$chr),]
#
## Read in a directory of pair files, returning both the Cy3 and Cy5 fluorescence in separate columns
# arrayIntensities <- loadSampleDirectory("Arrays", ndf, what="Cy3andCy5")
#
}
