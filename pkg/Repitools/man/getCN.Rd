\name{getCN}
\alias{getCN}
\title{Calculate Copy Number and Map To Enriched Regions}
\description{Based on the regions defined in the \code{inputs} element of \code{regionsInfo}, the function first calculates the log2 fold change between the treatment and control columns. Then, based on the regions defined in the \code{enriched} element, it finds which copy number region each enrichment region belongs to, and returns a matrix that holds the copy number state of every enrichment region.
}
\usage{
getCN(regionsInfo, readsInputs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regionsInfo}{A \code{list}, with elements \code{inputs} and \code{enriched}. \code{inputs} is a \code{data.frame} of with columns \code{chr}, \code{start}, \code{end}, that describes what regions from the input counts the copy number should be caculated from.  \code{enriched} is a \code{data.frame}, with the column names \code{chr}, \code{start}, \code{end}, that describes the regions is the enrichment count data that statistics are sought for.}
  \item{readsInputs}{A 2 column matrix of counts. Column 1 must have the counts for the control sample, and column 2 must be the treatment sample.}

}
\details{
 Regions of the input samples that have insufficient counts to calculate a log ratio (i.e. a zero count in either column) must not be present in the \code{readsInputs} matrix. The algorithm used to decide on the copy number regions is Circular Binary Segmentation (Olshen et al. 2004). If an enrichment region maps to a region that was not found to be in a copy number region from the input samples, it is assigned the copy number of 1.
}
\value{
	A 2-column matrix with the same number of rows as \code{regionsInfo$enriched}, that has the copy number of each region. The first column is for the control sample, so it will be all 1. The second column is for the treatment sample.
}
\author{Dario Strbenac}
\references{Olshen, A. B., Venkatraman, E. S., Lucito, R., and Wigler, M. (2004). Circular binary segmentation for the analysis of array-based DNA copy number data. \emph{Biostatistics} 5: 557-572}
\examples{
	inputs <- data.frame(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), start = c(1, 50001, 100001, 1, 10001), end = c(50000, 100000, 150000, 10000, 20000))
	enriched <- data.frame(chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2"), start = c(1, 1001, 2001, 10000001, 10005001, 9001, 9501), end = c(1000, 2000, 3000, 10001000, 10006000, 9500, 10000))
	counts <- matrix(c(25, 39, 3, 10, 22, 29, 38, 5, 19, 31), nrow = 5)
	colnames(counts) <- c("Control", "Treatment")
	getCN(list(inputs = inputs, enriched = enriched), counts)
}

