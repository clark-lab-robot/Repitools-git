\name{enrichmentCalc}
\alias{enrichmentCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate sequencing enrichment}
\description{Function to calculate enrichment over the whole genome of sequencing reads.
}
\usage{
enrichmentCalc(rs, organism, chrs, seqLen=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GRanges} object.}
  \item{organism}{The \code{BSgenome} object to calculate enrichment upon - only used to extract chromosome lengths.}
  \item{chrs}{The names of chromosomes that make up the genome. i.e. Allows exclusion of rand and hap chromosomes from the \code{organism} object.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
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
