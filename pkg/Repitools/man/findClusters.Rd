\name{findClusters}
\alias{findClusters}
\title{Find Clusters Epigenetically Modified Genes}
\description{Given a table of gene positions that has a score column, the function will find sort the genes into positional order and find consecutive windows of high or low scores.
}
\usage{
findClusters(statsTable, posCol, scoreCol, windowSize = 5, cutoff = 0.05, trend = c("down", "up"), nPermutations = 10, getFDRs = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statsTable}{A \code{data.frame} with (at least) columns \code{chr}, \code{start}, \code{end}, and a column of scores.}
  \item{posCol}{A number that gives the column in \code{statsTable} which contains the representative positions of the genes.}
  \item{scoreCol}{A number that gives the column in \code{statsTable} which contains the scores.}
  \item{windowSize}{The number of consecutive genes to consider windows over.}
  \item{cutoff}{The FDR level.}
  \item{trend}{Whether the clusters must have all positive scores (enrichment), or all negative scores (depletion).}
  \item{nPermutations}{How many random tables to generate to use in the FDR calculations.}
  \item{getFDRs}{If TRUE, will also return the table of FDRs at a variety of score cutoffs, from which the score cutoff for calling clusters is chosen.}
  \item{verbose}{Whether to print progress of computations.}
}
\details{
  First, the median over a window of size \code{wnidowSize} is calculated and then associated with the middle gene of the window. The order of the \code{statsTable} rows is randomised, and the median is again calculated, for every randomisation.

Next, the procedure for calling clusters is done at a range ofscore cutoffs.  For trend == "down", the cutoff starts at -1, increments by -0.05, and ends at -10. If the median associated with a particular gene is below the score cutoff, and the actual score of that gene is negative, then the cluster is progressively expanded outwards in both directions, for as long as the score remains negative. 

For trend == "up", the cutoff starts at 1, increments by 0.05, and ends at 10. If the median associated with a particular gene is above the score cutoff, and the actual score of that gene is positive, then the cluster is progressively expanded outwards in both directions, for as long as the score remains positive.

The first score cutoff to give an FDR below \code{cutoff} is chosen as the cutoff to use, and clusters are then called based on this cutoff.
}
\value{
	If \code{getFDRs} is FALSE, then only the \code{statsTable} table, with an additional column, \code{cluster}. If \code{getFDRs} is TRUE, then a list with elements
  \item{table}{The table \code{statsTable} with the additional column \code{cluster}} 
  \item{FDR}{The table of score cutoffs tried, and their FDRs.}
}
\author{Dario Strbenac}
\references{Saul Bert, in preparation}
\examples{
	chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
	starts <- sample(1:10000000, 25000, replace = TRUE)
	ends <- starts + 10000
	strands <- rep('+', 25000)
	genes <- data.frame(chr = sample(chrs, 25000, replace = TRUE), start = starts, end = ends, strand = '+', t = rnorm(100, 0, 2))
	findClusters(genes, 2, 5, trend = "up", nPermutations = 2)
}
