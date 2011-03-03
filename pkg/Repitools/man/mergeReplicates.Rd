\name{mergeReplicates}
\alias{mergeReplicates}
\alias{mergeReplicates,GenomeDataList-method}
\alias{mergeReplicates,GRangesList-method}
\title{Merge GRanges that are of replicate experiments.}
\description{
  A lane of next generation sequencing data can be stored as a \code{GRanges} object. Sometimes, a \code{GRangesList} of various lanes can have experimental replicates. This function allows the merging of such elements.
}
\usage{
  \S4method{mergeReplicates}{GenomeDataList}(reads, types, verbose = TRUE)
  \S4method{mergeReplicates}{GRangesList}(reads, types, verbose = TRUE)
}
\arguments{
  \item{reads}{A \code{\linkS4class{GenomeDataList}} or \code{\linkS4class{GRangesList}}.}
  \item{types}{A vector the same length as \code{reads}, that gives what type of experiment each element is of.}
  \item{verbose}{Whether to print the progess of processing.}
}

\value{A \code{\linkS4class{GRangesList}} with one element per experiment type.}
\details{The experiment type that each element of the merged list is of, is stored in the first element of the metadata list.}
\author{Dario Strbenac}
