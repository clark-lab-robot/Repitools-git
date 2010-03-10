\name{makeWindowLookupTable}
\alias{makeWindowLookupTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Using the output of 'annotationLookup', create a tabular storage of the indices}
\description{
To allow easy access to the probe-level data for either a gene, or an area of the promoter (over all genes), this routine takes the output of \code{annotationLookup} and organizes the indices into a table, one row for each gene and one column for each region of the promoter.
}
\usage{
makeWindowLookupTable(indexes, offsets, starts, ends)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexes}{a list of indices, e.g. \code{indexes} element from \code{annotationLookup} output}
  \item{offsets}{a list of offsets, e.g. \code{offsets} element from \code{annotationLookup} output}
  \item{starts}{a vector of starts}
  \item{ends}{a vector of ends}
}
\details{
  The vectors \code{starts} and \code{ends} (which should be the same length) determine the number of columns in the output matrix.
}
\value{
  A \code{matrix} with rows for each gene and columns for each bin of the promoter.  \code{NA} signifies that there is no probe in the given distance from a TSS.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}}
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
