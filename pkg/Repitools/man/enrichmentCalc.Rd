\name{enrichmentCalc}
\alias{enrichmentCalc}
\alias{enrichmentCalc,GenomeDataList,BSgenome-method}
\alias{enrichmentCalc,GenomeData,BSgenome-method}
\alias{enrichmentCalc,GRanges,BSgenome-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate sequencing enrichment}
\description{Function to calculate enrichment over the whole genome of sequencing reads.
}
\usage{
  \S4method{enrichmentCalc}{GenomeDataList,BSgenome}(rs, organism, seqLen=NULL, ...)
  \S4method{enrichmentCalc}{GenomeData,BSgenome}(rs, organism, seqLen=NULL, do.warn=FALSE)
  \S4method{enrichmentCalc}{GRanges,BSgenome}(rs, organism, seqLen=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeData}, \code{GenomeDataList}, or \code{GRanges} object.}
  \item{organism}{The \code{BSgenome} object to calculate enrichment upon - only used to extract chromosome lengths.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{do.warn}{Warn about coverage objects going over the edges of the chromosome lengths specified by the \code{organism} object}
  \item{...}{If wanting to pass the \code{do.warn} argument from the \code{GenomeDataList} version of the function to the \code{GenomeData} version, this allows it.}
}
\details{
  If \code{rs} has not been extended, it is extended using the supplied \code{seqLen} and then turned into a coverage object. The number of extended reads covering each base pair of the genome is then tabulated, and returned as a \code{data.frame}.
}
\value{
  A \code{data.frame} containing columns \code{coverage} and \code{bases}.
}
\author{Aaron Statham}
\examples{
#See examples in the manual
}
