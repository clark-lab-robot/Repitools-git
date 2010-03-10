\name{gcContentCalc}
\alias{gcContentCalc}
\title{Calculate the GC content of regions}
\description{Function to calculate GC content of windows surrounding the given positions.
}
\usage{
gcContentCalc(locations, window = 500, organism, verbose=FALSE, chunkSize=10000000)
}
\arguments{
  \item{locations}{A \code{data.frame}, with columns \code{chr} and \code{position}.}
  \item{window}{Bases around the locations that are in the window. Calculation will consider \code{windowSize / 2 - 1} bases upstream, and \code{windowSize / 2} bases downstream.}
  \item{organism}{The \code{BSgenome} object to calculate mappability upon.}
  \item{verbose}{Print details of processing.}
  \item{chunkSize}{The number of basepairs to process in a chunk, lower amounts use less memory, higher amounts may be faster}
}
\details{
The windows considered will be \code{windowSize / 2 - 1} bases upstream and \code{windowSize / 2} bases downstream of the given position, for each position. The value returned for each region is a percentage of bases in that region that are G or C.
}
\value{
  A vector of GC content percentages, one for each region.
}
\author{Dario Strbenac}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
gcContentCalc(TSSTable, 200, organism=Hsapiens)
}

