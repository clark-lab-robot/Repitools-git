\name{ChromaResults-class}
\docType{class}
\alias{ChromaResults}
\alias{class:ChromaResults}
\alias{ChromaResults-class}

\alias{blocks}
\alias{regions}
\alias{FDRTable}
\alias{cutoff}

\alias{show,ChromaResults-method}
\alias{blocks,ChromaResults-method}
\alias{regions,ChromaResults-method}
\alias{FDRTable,ChromaResults-method}
\alias{cutoff,ChromaResults-method}

\title{ChromaResults class}

\description{
	The \code{ChromaResults} class stores the results of a \code{\link{ChromaBlocks}} run.  
}

\section{Slots of a ChromaResults object}{
    \code{blocks:}{\code{GRanges} of the blocks used across the genome, with their calculated RPKM}
    \code{regions:}{\code{RangesList} of regions determined to be enriched}
    \code{FDRTable:}{\code{data.frame} showing the FDR at each cutoff tested}
    \code{cutoff:}{The cutoff used to determine enrichment}
}

\author{Aaron Statham}
\seealso{\code{\link{ChromaBlocks}}}
