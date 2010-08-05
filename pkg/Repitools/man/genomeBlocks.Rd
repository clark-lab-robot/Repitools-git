\name{genomeBlocks}
\alias{genomeBlocks}
\title{
Creates bins across a genome.
}
\description{Creates a compact \code{RangedData} representation of bins across specified chromosomes of a given genome}
\usage{
genomeBlocks(organism, chrs, width, spacing = width)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{A \code{BSgenome} object to create bins across}
  \item{chrs}{A \code{vector} containing which chromosomes to create bins across, may either be numeric indicies or chromosome names}
  \item{width}{The width in base pairs of each bin}
  \item{spacing}{The space between the centres of each adjacent bin. By default is equal to the \code{spacing} parameter, 
  which gives non-overlapping bins. Values larger than \code{spacing} will give overlapping bins, and values smaller than 
  \code{spacing} will give gaps between each bin.}
}
\value{Returns a \code{GRanges} object, compatible with direct usage in \code{annotationBlocksCounts}}
\author{Aaron Statham}
\seealso{\code{\link{annotationBlocksCounts}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
genomeBlocks <- function(organism, chrs, width, spacing=width) {
    chr.length <- seqlengths(organism)[chrs]
    windows <- do.call(c, lapply(names(chr.length), function(x) GRanges(seqnames = x, ranges = IRanges(start=seq.int(spacing/2, chr.length[x], spacing)-width/2+1, end=seq.int(spacing/2, chr.length[x], spacing)+width/2))))
    return(windows)
}

}
