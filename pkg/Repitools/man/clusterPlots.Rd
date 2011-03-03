\name{clusterPlots}
\alias{clusterPlots}
\alias{clusterPlots,ClusteredCoverageList-method}
\alias{clusterPlots,CoverageList-method}

\title{Visualisation of tables of feature coverages.}
\description{
  Takes the output of \link{featureCoverage}, or a modified version of it, and plots a heatmaps lineplots of clustered coverages.
}
\usage{
  \S4method{clusterPlots}{ClusteredCoverageList}(c.list, plot.ord = 1:length(c.list), plot.type = c("line", "heatmap"), cols = NULL, t.name = "Clustered Enrichment")
  \S4method{clusterPlots}{CoverageList}(c.list, scale = function(x) x, capQ = 0.95, cap.type = c("sep", "all"), n.clusters = 5, plot.ord = 1:length(c.list), expr, sort.data = NULL, sort.name = NULL, plot.type = c("line", "heatmap"), cols = NULL, t.name = "Clustered Enrichment")
}
\arguments{
  \item{c.list}{A \linkS4class{CoverageList} or \linkS4class{ClusteredCoverageList} object.}
  \item{scale}{A function to scale all the coverages by. Default : No scaling.}
  \item{capQ}{The quantile of coverages above which to make any bigger coverages equal to the quantile.}
  \item{cap.type}{If \code{"sep"}, then the cap quantile is calculated and applied to each coverage matrix separately. If \code{"all"}, then one cap quantile is calculated based on all of the matrices combined.}
  \item{n.clusters}{Number of clusters to find in the coverage data.}
  \item{plot.ord}{Order of the experiment types to plot.}
  \item{expr}{A vector of expression values.}
  \item{sort.data}{A vector of values to sort the features within a cluster on.}
  \item{sort.name}{Label to place under the \code{sort.data} plot.}
  \item{plot.type}{Style of plot to draw.}
  \item{cols}{The colours to use for the lines in the lineplot or intensities in the heatmap.}
  \item{t.name}{Title to use above all the heatmaps. Ignored when lineplots are drawn.}
}
\value{
  If called with a \code{CoverageList}, then a \code{ClusteredCoverageList} is returned. If called with a \code{ClusteredCoverageList}, then nothing is returned.
}
\details{
  A \code{ClusteredCoverageList} should be created by the user, if they wish to do some custom clustering on the coverage matrices. Otherwise, if the user is happy with k - means clustering, then the \code{CoverageList} object as output by \code{featureCoverage()} can be directly used. If called with a \code{CoverageList}, then the matrices for each coverage type are joined, scaling and capping is done, and k-means clustering is performed to get the cluster membership of each feature. Then a \code{ClusteredCoverageList} object is created, and used. If called with a \code{ClusteredCoverageList}, no scaling or capping is done, so it is the user's responsibility to give the scaled and capped coverage matrices when creating the \code{ClusteredCoverageList} object.

  For the plotting, the median expression level of each cluster is calculated. If a heatmap plot is being drawn, then a heatmap is drawn for every coverage matrix, side-by-side, and a plot of each feature's expression is put alongside the heatmaps. If additional sort vector was given, the data within clusters are sorted on this vector, then a plot of this data is made next to the expression plot.

  For a lineplot, a graph is drawn for each cluster, with the colours being the different coverage types. Because it makes sense that there will be more clusters than there are types of coverage (typically double to triple the number), the plots are not drawn side-by-side, as is the layout for the heatmaps. For this reason, sending the output to a PDF device is necessary.
}
\seealso{\code{\link{featureCoverage}} for generating coverage matrices.}
\author{Dario Strbenac}
