\name{annotationLookup}
\alias{annotationLookup}
\title{Forms a mapping between probes on a tiling array and regions of a fixed size around chromosomal co-ordinates}
\description{
Starting from a table of genome locations for probes (see \code{probes} below) and a table of TSSs for genes (see \code{annotation} below), this procedure forms a list structure that contains the indices to map from one to the other.
}
\usage{
annotationLookup(probes, annotation, bpUp, bpDown, probeIndex = NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{probe genomic locations, a \code{data.frame} with required elements \code{chr} (character), \code{position} (integer)}
  \item{annotation}{a \code{data.frame} with required elements \code{chr} (character), \code{start} (integer), \code{end} (integer), and optional elements \code{strand} ("+" or "-" or "*") and \code{name} (character). If column \code{strand} is not present, the region is treated as being on both strands (e.g. CpG islands.)}
  \item{bpUp}{The number of bases upstream to look}
  \item{bpDown}{The number of bases downstream to look}
  \item{probeIndex}{a probe index, default \code{NULL}}
  \item{verbose}{whether to print progress to screen}
}
\details{
This function is a wrapper for the generic function \code{annotationBlocksLookup} which can handle annotations of varying sizes. \code{annotationLookup} is appropriate where you wish to map probes that are within a fixed distance of points of annotation e.g gene transcription start sites. Even if strand information is given for probes, it is ignored.
}
\value{
  A list with elements
  \item{indexes}{a list for each gene in \code{annotation} table, giving a vector of indices to the probe table}
  \item{offsets}{a list for each gene in \code{annotation} table, giving a vector (corresponding to \code{indexes}) of offsets (to the TSS) for each probe}
}
\author{Aaron Statham, Mark Robinson}
\seealso{\code{\link{annotationBlocksLookup}}, \code{\link{makeWindowLookupTable}}}
\examples{

# create example set of probes and gene start sites
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX", strand = '-')
genes <- data.frame(chr="chrX", start=c(2100, 1000), end = c(3000, 2200), strand=c("+","-"))
rownames(genes) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
annotationLookup(probeTab, genes, 500, 500)
}
