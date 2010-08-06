\name{writeWig}
\alias{writeWig}
\alias{writeWig.AffymetrixCelSet}
\alias{writeWig.GRangesList}
\title{Writes sequencing data out into a wiggle files}
\description{Writes sequencing data out into a wiggle files}
\usage{
\method{writeWig}{AffymetrixCelSet}(rs, design=NULL, log2adjust=TRUE, verbose=TRUE, ...)
\method{writeWig}{GRangesList}(rs, seqLen, design=NULL, sample=20, dropZero=TRUE, normalise=TRUE, verbose=TRUE, ...)
}
\arguments{
  \item{rs}{The sequencing or array data.}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{verbose}{Whether to write progress to screen}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{sample}{At what basepair resolution to sample the genome at}
  \item{dropZero}{Whether to write zero values to the wiggle file - TRUE saves diskspace}
  \item{normalise}{Whether to normalise each lane to its total number of reads, TRUE is suggested}
  \item{...}{Keeps S3 methods happy}
}
\details{
A wiggle file is created for each column in the design matrix (if design is left as NULL, then a file 
is created for each array/lane of sequencing). The filenames are given by the column names of the design
matrix, and if ending in "gz" will be written out as a gzfile. 
}
\value{Wiggle file(s) are created}
\author{
Aaron Statham
}
\examples{
#See examples in the manual
}
