\name{clusterPlots}
\alias{clusterPlots}
\alias{clusterPlots,list-method}
\title{Create graphical summaries of averaged read coverage across a promoter associated with expression.}
\description{
Clusters the coverage around gene TSSs for a set of ChIP-seq experiments and plots either the heatmaps of these alongside a vector of gene expressions, and gene widths or per - cluster lineplots.
}
\usage{
  \S4method{clusterPlots}{list}(coverageList, capQuantile = 0.95, nClusters = 5, expr, sortData = NULL, sortName = NULL, plotType = c("line", "heatmap"), cols = NULL, tName = "Clustered Enrichment")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coverageList}{A \code{list} object of matrices of coverage. Each matrix has rows for each gene, and columns for each sampling position. The row names of the matrix should be the gene names.}
  \item{capQuantile}{Intensities above this quantile will be set to this quantile.}
  \item{nClusters}{The number of clusters to find. k-means method is used.}
  \item{expr}{A vector of expression values for the genes. Must have same number of rows as \code{geneAnno}. The relationship to the gene annotation table is by the order of expression values corresponding to the order of rows of the annotation.}
  \item{sortData}{A vector as long as the number of rows of coverage list. Gives some measurement of the genes involved, like width, for example.}
  \item{sortName}{An intuitive name for the sort data. Will be used as an axis label for the plot of the sort data.}
  \item{plotType}{Type of binned plot to draw.}
  \item{cols}{A vector of colours to use. The length of this also determines how many levels of intensities there are.}
  \item{tName}{Title to put above the heatmaps.}
}
\details{
  The \code{coverageList} list should be the list that is output from the function \code{getTSScoverage}. The code has been divided into two functions because the coverage and smoothing calculations take time, and by having the heatmap plotting separate, the user can change the number of clusters or the colour scheme, and see the results quickly, compared to the smoothed coverage having to be executed each time.

  If \code{plotType} is 'line', then the average coverage at each position for each mark is drawn as a line. Each plot is for one cluster. If a heatmap is requested, then the heatmap image contains a row for every gene. Each heatmap panel is for one mark.
}
\value{
  A list the same length at the number of clusters (invisibly). Each element contains gene names that belong to that cluster. The name of each element is the median expression level of that cluster.
}
\seealso{\code{\link{getTSScoverage}}}
\author{Dario Strbenac}
\examples{
  # Not run.
  #
  # geneAnno <- read.csv("geneAnno.csv")
  # load("exampleGRangesList.Rdata")
  # load("exprVector.Rdata")
  # tablesList <- getTSScoverage(exampleGRangesList, exptTypes = c("MBD2", "MBD2", "H3K27me3", "H3K27me3"), geneAnno = geneAnno)
  # clusteredPlots(tablesList, geneAnno, nClusters = 3, expr = exprVector, tName = "A Sample Plot")
}

