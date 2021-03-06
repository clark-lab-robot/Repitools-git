\name{annotationBlocksLookup}
\alias{annotationBlocksLookup}
\title{Forms a mapping between probe locations and chromosomal blocks (regions)}
\description{
Starting from a table of genome locations for probes (see \code{probes} below) and a table of regions of interest (see \code{annotation} below), this procedure forms a list structure that contains the indices to map from one to the other.
}
\usage{
annotationBlocksLookup(probes, annotation, probeIndex = NULL, verbose=TRUE)
}
\arguments{
  \item{probes}{probe genomic locations, a \code{data.frame} with required elements \code{chr} (character), \code{position} (integer), and optionally \code{strand} ("+" or "-")}
  \item{annotation}{a \code{data.frame} with required elements \code{chr} (character), \code{start} (integer), \code{end} (integer) and optional elements \code{strand} ("+" or "-" or "*") and \code{name} (character). If column \code{strand} is not present, the region is treated as being on both strands (e.g. CpG islands.)}
  \item{probeIndex}{a probe index, default \code{NULL}}
  \item{verbose}{whether to print progress to screen}
}
\details{
Strandedness of probes is ignored, even if it is given. 
}
\value{
  A list with elements
  \item{indexes}{a list for each row in \code{annotation} table, giving a vector of indices to the probe table}
  \item{offsets}{a list for each row in \code{annotation} table, giving a vector (corresponding to \code{indexes}) of offsets (from the beginning of the annotation block) for each probe}
}
\author{Aaron Statham, Mark Robinson}
\seealso{\code{\link{annotationLookup}} which simplifies annotation lookups for constant sized regions}
\examples{
# create example set of probes and gene start sites
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX", strand="+")
genes <- data.frame(chr="chrX", start=c(2100,2200), end=c(2500, 2400), strand=c("+","-"))
rownames(genes) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
annotationBlocksLookup(probeTab, genes)
}
