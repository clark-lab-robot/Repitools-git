\name{ChromaBlocks}
\alias{ChromaBlocks}
\title{A function to find areas of enrichment in sequencing data}
\description{
For each region of interest or TSS, this routine interrogates probes or sequence data for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.  Regions can be surroundings of TSSs, or can be a table of user - specified regions. The function determines if the \code{start} and \code{end} columns of the \code{coordinatesTable} should be used as regions or as TSSs, based on the \code{useAsRegions} flag.
}
\usage{
ChromaBlocks(rs.ip, rs.input, organism, chrs, ipWidth=100, inputWidth=500, preset=NULL, blockWidth=NULL, cutoff=NULL, minBlocks=NULL, extend=NULL, verbose=TRUE, seqLen=NULL)
}
\arguments{
  \item{rs.ip}{A \code{GenomeDataList} object containing reads from the Immunoprecipited sample}
  \item{rs.input}{A \code{GenomeDataList} object containing reads from the Input (unenriched) sample}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{chrs}{An \code{integer vector} with the indicies of the chromosomes of the \code{organism} object to analyse}
  \item{ipWidth}{Size in basepairs of the windows to use for the IP samples}
  \item{inputWidth}{Size in basepairs of the windows to use for the Input samples}
  \item{preset}{Either "small", "large" to use cutoffs described in Hawkins et al or \code{NULL} (where \code{blockWidth}, \code{cutoff} and \code{minBlocks} must be specified)}
  \item{blockWidth}{Number of adjacent blocks to consider at once}
  \item{cutoff}{The minimum average RPKM of \code{blockWidth} adjacent blocks}
  \item{minBlocks}{The minimum number of blocks required above \code{cutoff}}
  \item{extend}{Optional: whether to extend significant blocks until adjacent blocks are less than this value}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
}
\value{
Havent worked out what to return yet
}
\author{Aaron Statham}
