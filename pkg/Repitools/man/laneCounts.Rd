\name{laneCounts}
\alias{laneCounts}
\alias{laneCounts,GenomeDataList-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to count the total number of ChIP-seq reads in each lane.
}
\description{
A function to count the total number of ChIP-seq reads in each lane.}
\usage{
  \S4method{laneCounts}{GenomeDataList}(cs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{A \code{GenomeDataList}, with each element corresponding to the reads for a lane of sequencing.
}
}
\value{
	An \code{integer vector}.
}
\author{Aaron Statham}
\examples{
#See the manual
}
