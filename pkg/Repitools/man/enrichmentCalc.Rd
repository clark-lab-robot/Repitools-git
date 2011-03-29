\name{enrichmentCalc}
\alias{enrichmentCalc}
\alias{enrichmentCalc,GenomeDataList-method}
\alias{enrichmentCalc,GRangesList-method}
\alias{enrichmentCalc,GRanges-method}
\title{Calculate sequencing enrichment}
\description{Function to calculate enrichment over the whole genome of sequencing reads.}

\section{Usage}{
  \describe{
    \code{enrichmentCalc(x, organism, seqLen=NULL)}
  }
}
\section{Arguments}{
  \describe{
      \item{x:}{A \code{GRangesList}, \code{GRanges}, or \code{GenomeDataList} object.}
      \item{organism:}{The \code{BSgenome} object to calculate enrichment upon - only used to extract chromosome lengths.}
      \item{seqLen:}{If sequencing reads need to be extended, the fragment size to be used}
  }
}

\section{Details}{
  \describe{
    If \code{x} has not been extended, it is extended using the supplied \code{seqLen} and then turned into a coverage object. The number of extended reads covering each base pair of the genome is then tabulated, and returned as a \code{data.frame}.
  }
}
\section{Value}{
  \describe{
    For the \code{GRanges} method, \code{data.frame} containing columns \code{coverage} and \code{bases}.
    For the \code{GRangesList} and \code{GenomeDataList} methods, a list of such \code{data.frame}s.
  }
}
\author{Aaron Statham}
\examples{
#See examples in the manual
}
