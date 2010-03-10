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
  \item{probes}{probe genomic locations, a \code{data.frame} with required elements \code{chr} (character) and \code{position} (integer)}
  \item{annotation}{a \code{data.frame} with required elements \code{chr} (character) and \code{position} (integer) and optional elements \code{strand} ("+" or "-") and \code{name} (character)}
  \item{bpUp}{The number of bases upstream to look}
  \item{bpDown}{The number of bases downstream to look}
  \item{probeIndex}{a probe index, default \code{NULL}}
  \item{verbose}{whether to print progress to screen}
}
\details{
This function is a wrapper for the generic function \code{annotationBlocksLookup} which can handle annotations of varying sizes. \code{annotationLookup} is appropriate where you wish to map probes that are within a fixed distance of points of annotation e.g gene transcription start sites.
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
probeTab <- data.frame(position=seq(1000,3000,by=200), chr="chrX")
geneTSS <- data.frame(chr="chrX", position=c(2100,2200), strand=c("+","-"))
rownames(geneTSS) <- paste("gene",1:2,sep="")

# Call annotationLookup() and look at output
aL <- annotationLookup(probeTab, geneTSS, 500, 500)
print(aL)

# Store the results of annotationLookup() in a convenient tabular format
lookupTab <- makeWindowLookupTable(aL$indexes, aL$offsets, starts=seq(-400,200,by=200), ends=seq(-200,400,by=200))
print(lookupTab)

}
