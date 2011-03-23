\name{mappabilityCalc}
\alias{mappabilityCalc}
\alias{mappabilityCalc,GRanges-method}
\alias{mappabilityCalc,data.frame-method}
\title{Calculate The Mappability of a Region}
\description{Function to calculate mappability of windows}
\usage{
  \S4method{mappabilityCalc}{GRanges}(x, organism)
  \S4method{mappabilityCalc}{data.frame}(x, window=500, organism)
}
\arguments{
  \item{x}{A \code{GRanges} object or a \code{data.frame}, with columns \code{chr} and either \code{position} or \code{start}, \code{end} and \code{strand}.}
  \item{window}{Bases around the locations that are in the window. Calculation will consider \code{windowSize/2} bases upstream, and \code{windowSize/2-1} bases downstream.}
  \item{organism}{The \code{BSgenome} object to calculate mappability upon.}
}
\details{
The windows considered will be \code{windowSize/2} bases upstream and \code{windowSize/2-1} bases downstream of the given position, for each position. The value returned for each region is a percentage of bases in that region that are not N (any base in IUPAC nomenclature).
}
\value{
  A vector of mappability percentages, one for each region.
}
\author{Aaron Statham}
\examples{
# require(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
# TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
# mappabilityCalc(TSSTable, 200, organism=Hsapiens36bp)
}

