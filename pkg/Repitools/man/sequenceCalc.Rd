\name{sequenceCalc}
\alias{sequenceCalc}
\alias{sequenceCalc,GRanges-method}
\alias{sequenceCalc,data.frame-method}
\title{Find occurences of a DNA pattern}
\description{Function to find all occurrences of a DNA pattern in given locations.}
\usage{
  \S4method{sequenceCalc}{GRanges}(locations, organism, pattern, fixed=TRUE, Nmask=FALSE, positions=FALSE, verbose=FALSE, chunkSize=10000000, ...)
  \S4method{sequenceCalc}{data.frame}(locations, window=500, organism, pattern, fixed=TRUE, Nmask=FALSE, positions=FALSE, verbose=FALSE, chunkSize=10000000, ...)
}
\arguments{
  \item{locations}{A \code{data.frame}, with columns \code{chr} and \code{position}, or instead of the column \code{position} there can be columns \code{start}, \code{end}, and \code{strand}, or a \code{GRanges} object of the regions.}
  \item{window}{Bases around the locations that are in the window. Calculation will consider \code{windowSize / 2 - 1} bases upstream, and \code{windowSize / 2} bases downstream.}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{pattern}{The \code{DNAString} to search for.}
  \item{fixed}{Whether to allow degenerate matches.}
  \item{Nmask}{Whether to mask out Ns in the \code{organism} sequence before searching - recommended when \code{fixed=FALSE}.} 
  \item{positions}{If \code{TRUE}  \code{FALSE}}
  \item{verbose}{Print details of processing.}
  \item{chunkSize}{The number of basepairs to process in a chunk, lower amounts use less memory, higher amounts may be faster.}
  \item{...}{Not used.}
}
\details{If the version of the data frame with the start, end, and strand columns is given, the window will be created around the TSS.}
\value{
  If \code{positions} is \code{TRUE}, a list of vectors of positions of matches in relation to the \code{locations}, otherwise a \code{vector} of the number of matches for each element of \code{locations}.
}
\seealso{\code{\link{cpgDensityCalc}}, \code{\link{mappabilityCalc}}, \code{\link{gcContentCalc}}}
\author{Aaron Statham}
\examples{
require(BSgenome.Hsapiens.UCSC.hg18)
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
sequenceCalc(TSSTable, 600, organism=Hsapiens, pattern=DNAString("CG"))
}
