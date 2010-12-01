\name{clusteredHeatmap}
\alias{clusteredHeatmap}
\alias{clusteredHeatmap,GRangesList-method}
\title{Create heatmaps of averaged read coverage across a promoter and and plot with expression.}
\description{
Clusters the coverage around gene TSSs for a set of ChIP-seq experiments and plots the heatmaps of these alongside a vector of gene expressions, and gene widths.
}
\usage{
  \S4method{clusteredHeatmap}{GRangesList}(readsIPs, exptTypes, geneAnno, upPct = 50, downPct = 100, cvgResolution = 10, smoothingWidth = 2001, capQuantile = 0.95, nClusters = 5, expr, cols = colorpanel(64, "blue", "white", "red"), tName = "Clustered Enrichment")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readsIPs}{A \code{GRangesList} object of lanes of sequencing data.}
  \item{exptTypes}{A character vector of the same length as \code{readsIPs}, that gives what epigenetic mark each sample is of.}
  \item{geneAnno}{A \code{data.frame} with at least columns \code{chr}, \code{start}, \code{end}, and \code{strand}.}
  \item{upPct}{How many percent upstream of the TSS to plot.}
  \item{downPct}{How many percent upstream of the TSS to plot.}
  \item{cvgResolution}{How often to take a measurement.}
  \item{smoothingWidth}{Must be an odd number. The sliding window size used for smoothing.}
  \item{capQuantile}{Intensities above this quantile will be set to this quantile.}
  \item{nClusters}{The number of clusters to find. k-means method is used.}
  \item{expr}{A vetor of expression values for the gemes. Must have same number of rows as \code{geneAnno}. The relationship to the gene annotation table is by the order of expression values corresponding to the order of rows of the annotation.}
  \item{cols}{A vector of colours to use. The length of this also determines how many levels of intensities there are.}
  \item{tName}{Title to put above the heatmaps.}
}
\details{
	The region of coverage of \code{upPct} percent of the width of each gene upstream to \code{downPct} percent of the width of each gene downstream is extracted. Flanking regions of half the smoothing width are also extracted at each end, so that even the first or last base's coverage of the region can be smoothed. A sequence of sampling points is determined by \code{upPct}, \code{downPct}, and \code{cvgResolution}. Starting at \code{upPct} of the width of each gene upstream, a sample of the smoothed coverage is made every \code{cvgResolution} percent, until the last sample is made at \code{downPct} percent of the gene width downstream. Genes that are close to the edge of chrosomomes and don't have coverage data at all of the sampling points are excluded from the heatmap.
}
\value{
  A plot to the current device.
}
\author{Dario Strbenac}
\examples{
  # Not run.
  # 
  # geneAnno <- read.csv("geneAnno.csv")
  # load("exprVector.Rdata")
  # load("exampleGRangesList.Rdata")
  # clusteredHeatmap(exampleGRangesList, exptTypes = c("MBD2", "MBD2", "H3K27me3", "H3K27me3") geneAnno, expr = exprVector, tName = "A Sample Plot")
}

