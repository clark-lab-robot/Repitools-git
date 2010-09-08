\name{annotationCounts}
\alias{annotationCounts}
\alias{annotationCounts,GenomeDataList-method}
\alias{annotationCounts,GRangesList-method}
\title{Counts the number of sequencing reads surrounding supplied annotations}
\description{
Counts the number of sequencing reads surrounding supplied annotations
}
\usage{
  \S4method{annotationCounts}{GenomeDataList}(rs, annotation, bpUp, bpDown, seqLen = NULL, verbose=TRUE)
  \S4method{annotationCounts}{GRangesList}(rs, annotation, bpUp, bpDown, seqLen = NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeData} or \code{GenomeDataList} object}
  \item{annotation}{a \code{data.frame} with required elements \code{chr} (character),  \code{start} (integer), \code{end} (integer) and optional elements \code{strand} ("+" or "-") and \code{name} (character)}
  \item{bpUp}{The number of bases upstream to look}
  \item{bpDown}{The number of bases downstream to look}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{verbose}{logical, whether to print progress}
}
\value{A \code{matrix} of counts is returned, one column per \code{GenomeData} object supplied in \code{rs}, one row per \code{annotation} supplied}
\author{Aaron Statham}
\seealso{\code{\link{annotationBlocksCounts}}, \code{\link{genomeBlocks}}}
\examples{#See the manual}
