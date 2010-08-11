\name{cpgDensityCalc}
\alias{cpgDensityCalc}
\alias{cpgDensityCalc.GenomeDataList}
\alias{cpgDensityCalc.GenomeData}
\alias{cpgDensityCalc.GRanges}
\alias{cpgDensityCalc.data.frame}
\title{Calculate CpG Density in a Window}
\description{Function to calculate CpG density around a position.
}
\usage{
\method{cpgDensityCalc}{GRanges}(rs, seqLen, ...)
\method{cpgDensityCalc}{data.frame}(locations, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, chunkSize=10000000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locations}{A \code{data.frame}, with columns \code{chr} and \code{position}, or a GRanges object.}
  \item{window}{Bases around the locations that are in the window. Calculation will consider \code{windowSize / 2 - 1} bases upstream, and \code{windowSize / 2} bases downstream.}
  \item{wFunction}{Weighting function to use. Can be \code{none}, \code{linear}, \code{log}, or \code{exp}}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{rs}{A \code{GenomeData} or \code{GenomeDataList} object.}
  \item{seqLen}{The fragment size of the sequence reads in \code{rs}, doubles as the \code{window} to be used.}
  \item{verbose}{Print details of processing.}
  \item{chunkSize}{The number of basepairs to process in a chunk, lower amounts use less memory, higher amounts may be faster}
  \item{...}{Arguments passed into the \code{data.frame} method}
}
\details{
  If the version of the data frame with the start, end, and strand columns is given, the window will be created around the TSS.
  
  For weighting scheme \code{none}, this is equivalent to the number of CG matches in the region. For \code{linear} weighting, each match is given a score \code{1/x} where \code{x} is the number of bases from the postition that the match occurred, and the scores are summed. For exponential weighting and logarithmic weighting, the idea is similar, but the scores decay exponentially (\code{exp^-5x/windowSize}) and logarithmically (\code{log2(2 - (distancesForRegion / windowSize)}).
}
\value{
  A vector of CpG densities for each region.
}
\author{Dario Strbenac}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
cpgDensityCalc(TSSTable, 600, organism=Hsapiens)
}
