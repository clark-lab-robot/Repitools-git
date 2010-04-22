\name{findEADs}
\alias{findEADs}
\title{Find Epigenetically Activated Domains}
\description{Given a table of gene positions that has a t - statistic column, the function will find sort the genes into positional order and find consecutive windows of high expression.
}
\usage{
findEADs(statsTable, windowSize = 5, cutoff = 0.05, getFDRs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statsTable}{A \code{data.frame} with (at least) columns \code{chr}, \code{start}, \code{end}, and \code{t}.}
  \item{windowSize}{The number of consecutive genes to consider windows over.}
  \item{cutoff}{The FDR level.}
  \item{getFDRs}{If TRUE, will also return the table of FDRs at a variety of t - statistic cutoffs, from which the t - statistic cutoff for calling EADs is chosen.}
}
\details{
  First, the median over a window of size \code{wnidowSize} is calculated and then associated with the middle gene of the window. The order of the \code{statsTable} rows is randomised, and the median is again calculated.

Next, the procedure for calling EADs is done at a range of t - statistic cutoffs, starting at 1, incrementing by 0.05, and ending at 10. If the median associated with a particular gene is above the t - statistic cutoff, and the actual t - statistic of that gene is positive, then the EAD is progressively expanded outwards in both directions, for as long as the t - statistic remains positive. This contiguous region is an EAD. The number of EADs in the real \code{statsTable} and the permuted table are calculated. The first t - statistic cutoff to give an FDR below \code{cutoff} is chosen as the cutoff to use, and EADs are then called based on this cutoff.
}
\value{
	If \code{getFDRs} is FALSE, then only the \code{statsTable} table, with an additional column, \code{EAD}. If \code{getFDRs} is TRUE, then a list with elements
  \item{table}{The table \code{statsTable} with the additional column \code{EAD}} 
  \item{FDR}{The table of t - statistic cutoffs tried, and their FDRs.}
}
\author{Dario Strbenac}
\references{Saul Bert, in preparation}
\examples{
	chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
	starts <- sample(1:10000000, 25000, replace = TRUE)
	ends <- starts + 10000
	genes <- data.frame(chr = sample(chrs, 25000, replace = TRUE), start = starts, end = ends, t = rnorm(100, 0, 2))	
	findEADs(genes)
}
