\name{plotClusters}
\alias{plotClusters}
\title{Plot Scores of Cluster Regions}
\description{Given a table of gene positions that has a score column, the function will make a series of bar chart plots, one for each cluster.
}
\usage{
plotClusters(clustersTable, posCol, scoreCol, scoreType, ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clustersTable}{A \code{data.frame} with (at least) columns \code{chr}, \code{start}, \code{end}, and \code{cluster}. Also a column of scores.}
  \item{scoreCol}{A number that gives the column in \code{statsTable} which contains the scores.}
  \item{posCol}{A number that gives the column in \code{statsTable} which contains the representative positions of the genes.}
  \item{scoreType}{A character string of the name of the score. To be used as the y - axis label in the plot.}
  \item{ylim}{The y - axis limits of the plots to be used.}
}
\value{
	No return value.
}
\author{Dario Strbenac}
\examples{
	chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
	starts <- sample(1:10000000, 25000, replace = TRUE)
	ends <- starts + 10000
	strands <- rep('+', 25000)
	genes <- data.frame(chr = sample(chrs, 25000, replace = TRUE), start = starts, end = ends, strand = '+', t = rnorm(100, 0, 2))
	genes <- genes[order(genes$chr, genes$start), ]
	clusterIDs <- rep(0, 25000)
	clusterIDs[96:100] <- 1
	clusterIDs[496:500] <- 2
	clusterIDs[996:1000] <- 3
	genes$cluster = clusterIDs
	plotClusters(genes, 2, 5, "t - statistic", ylim = c(0, 10))
}
