\name{cpgDensityCalc}
\alias{cpgDensityCalc}
\alias{cpgDensityCalc,GenomeDataList-method}
\alias{cpgDensityCalc,GRangesList-method}
\alias{cpgDensityCalc,GRanges-method}
\alias{cpgDensityCalc,data.frame-method}
\title{Calculate CpG Density in a Window}
\description{Function to calculate CpG density around a position.
}
\usage{
  \S4method{cpgDensityCalc}{GenomeDataList}(x, ...)
  \S4method{cpgDensityCalc}{data.frame}(x, organism, ...)
  \S4method{cpgDensityCalc}{GRangesList}(x, window=500, wFunction=c("none","linear","exp","log"), organism, verbose=TRUE, seqLen=NULL)
  \S4method{cpgDensityCalc}{GRanges}(x, window=500, wFunction=c("none","linear","exp","log"), organism, verbose=TRUE, seqLen=NULL)
}
\arguments{
  \item{x}{A \code{data.frame}, with columns \code{chr} and \code{position}, a \code{GRangesList} object, \code{GenomeDataList} or \code{GRanges}.}
  \item{window}{Bases around the locations that are in the window. Calculation will consider \code{windowSize/2-1} bases upstream, and \code{windowSize/2} bases downstream.}
  \item{wFunction}{Weighting function to use. Can be \code{none}, \code{linear}, \code{log}, or \code{exp}}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{seqLen}{The fragment size of the sequence reads in \code{x}}
  \item{verbose}{Print details of processing.}
  \item{...}{Arguments passed into the \code{GRanges} method}
}
\details{
  If the version of the data frame with the start, end, and strand columns is given, the window will be created around the TSS.
  
  For weighting scheme \code{none}, this is equivalent to the number of CG matches in the region. For \code{linear} weighting, each match is given a score \code{1/x} where \code{x} is the number of bases from the postition that the match occurred, and the scores are summed. For exponential weighting and logarithmic weighting, the idea is similar, but the scores decay exponentially (\code{exp^-5x/windowSize}) and logarithmically (\code{log2(2 - (distancesForRegion / windowSize)}).
}
\value{
  A \code{numeric} vector of CpG densities for each region.
}
\author{Dario Strbenac}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
cpgDensityCalc(TSSTable, organism=Hsapiens, window=600)
}
