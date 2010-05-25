\name{annotationBlocksCounts}
\alias{annotationBlocksCounts}
\title{Counts the number of sequencing reads within supplied genomic blocks}
\description{
Counts the number of sequencing reads within supplied genomic blocks
}
\usage{
annotationBlocksCounts(rs, annotation, seqLen = NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeData} or \code{GenomeDataList} object}
  \item{annotation}{a RangedData object (preferable) or a \code{data.frame} with required elements \code{chr} (character),  \code{start} (integer), \code{end} (integer) and optional elements \code{strand} ("+" or "-") and \code{name} (character)}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{verbose}{logical, whether to print progress}
}
\value{A \code{matrix} of counts is returned, one column per \code{GenomeData} object supplied in \code{rs}, one row per \code{annotation} supplied}
\author{Aaron Statham}
\seealso{\code{\link{annotationCounts}}, \code{\link{genomeBlocks}}}
\examples{#See the manual}
