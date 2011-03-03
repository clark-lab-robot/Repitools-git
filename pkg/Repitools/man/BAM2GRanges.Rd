\name{BAM2GenomicRanges}

\alias{BAM2GRanges}
\alias{BAM2GRangesList}
\alias{BAM2GRanges,character-method}
\alias{BAM2GRangesList,character-method}
\title{Read in a (list of) BAM file(s) into a GRanges(List) object.}
\description{
  A wrapper script for coverting the contents of BAM files for use with
  \code{GenomicRanges} classes.
}
\usage{
  \S4method{BAM2GRanges}{character}(path, what = c("rname", "strand", "pos", "qwidth"),
      flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
      verbose = TRUE)
  \S4method{BAM2GRangesList}{character}(paths, what = c("rname", "strand", "pos", "qwidth"),
      flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
      verbose = TRUE)
}
\arguments{
  \item{path}{A character vector of length 1. The path of the BAM file.}
  \item{paths}{A character vector of possibly any length. The paths of the
               BAM files.}
  \item{what}{What attributes of a read to retain. See
              \code{\link[Rsamtools]{scanBam}} and the \code{value} section.}
  \item{flag}{What kinds of reads to retain.  See
              \code{\link[Rsamtools]{ScanBamParam}} and the \code{flag} argument.}
  \item{verbose}{Whether to print the progess of processing.}
}
\value{
  For the single pathname method; a GRanges object. For the multiple pathnames method;
  a GRangesList object.
}
\author{Dario Strbenac}
