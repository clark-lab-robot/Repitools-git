\name{ChromaBlocks}
\alias{ChromaBlocks}
\alias{ChromaBlocks,GenomeDataList,GenomeDataList-method}
\alias{ChromaBlocks,GRangesList,GRangesList-method}
\title{A function to find areas of enrichment in sequencing data}
\description{
This function discovers regions of enrichment in ChIP-seq data, using the method described in Hawkins RD. et al 2010 Cell Stem Cell.
}
\usage{
\S4method{ChromaBlocks}{GenomeDataList,GenomeDataList}(rs.ip, rs.input, ...)
\S4method{ChromaBlocks}{GRangesList,GRangesList}(rs.ip, rs.input, organism, chrs, ipWidth=100, inputWidth=500, preset=NULL, blockWidth=NULL, minBlocks=NULL, extend=NULL, cutoff=NULL, FDR=0.01, nPermutations=5, nCutoffs=20, cutoffQuantile=0.98, verbose=TRUE, seqLen=NULL)
}
\arguments{
  \item{rs.ip}{A \code{GenomeDataList} or \code{GRangesList} object containing reads from the Immunoprecipited sample}
  \item{rs.input}{A \code{GenomeDataList} or \code{GRangesList} object containing reads from the Input (unenriched) sample}
  \item{organism}{The \code{BSgenome} object to calculate CpG density upon.}
  \item{chrs}{An \code{integer vector} with the indicies of the chromosomes of the \code{organism} object to analyse}
  \item{ipWidth}{Size in basepairs of the windows to use for the IP samples}
  \item{inputWidth}{Size in basepairs of the windows to use for the Input samples}
  \item{preset}{Either "small", "large" to use cutoffs described in Hawkins et al or \code{NULL} (where \code{blockWidth}, \code{cutoff} and \code{minBlocks} must be specified)}
  \item{blockWidth}{Number of adjacent blocks to consider at once}
  \item{minBlocks}{The minimum number of blocks required above \code{cutoff}}
  \item{extend}{Optional: whether to extend significant blocks until adjacent blocks are less than this value}
  \item{cutoff}{Optional: the cutoff to use to call regions. If left as \code{NULL} a cutoff will be chosen which satisfied the specified FDR}
  \item{FDR}{The target False Discovery Rate; If \code{cutoff} is not supplied, one will be chosen to satisfy this value}
  \item{nPermutations}{The number of permutations of the data to determine the \code{cutoff} at the supplied \code{FDR}}
  \item{nCutoffs}{The number of different cutoffs to try to satisfy the \code{FDR}, a higher value will give finer resolution but longer processing time}
  \item{cutoffQuantile}{The quantile of the RPKM to use as the maximum cutoff tried; a higher value will give lower resolution but may be needed if a \code{cutoff} satisfying the \code{FDR} cannot be determined with the default value}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{...}{Arguments passed from the \code{GenomeDataList} method to the \code{GRangesList} method.}
}
\value{
  A \code{\link{ChromaResults}} object.
}
\author{Aaron Statham}
\seealso{\code{\link{ChromaResults}}}
