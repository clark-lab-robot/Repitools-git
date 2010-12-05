\name{getTSScoverage}
\alias{getTSScoverage}
\alias{getTSScoverage,GRangesList-method}
\title{Find smoothed coverage around TSS.}
\description{
Creates a list of tables of coverages, one table for each sample type.
}
\usage{
  \S4method{getTSScoverage}{GRangesList}(readsIPs, seqLen = 300, exptTypes, geneAnno, up = 50, down = 100, distType = c("percent", "base"), cvgResolution = 10, smoothingWidth = 2001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readsIPs}{A \code{GRangesList} object of lanes of sequencing data.}
  \item{seqLen}{Length to extend the short reads to.}
  \item{exptTypes}{A character vector of the same length as \code{readsIPs}, that gives what epigenetic mark each sample is of.}
  \item{geneAnno}{A \code{data.frame} with at least columns \code{chr}, \code{start}, \code{end}, \code{strand}, and \code{name}.}
  \item{up}{How far upstream of the TSS to find coverage of.}
  \item{down}{How far upstream of the TSS to find coverage of.}
  \item{distType}{Whether the \code{up} and \code{down} values refer to a percentage of the gene width, or the number of bases.}
  \item{cvgResolution}{How often to take a measurement.}
  \item{smoothingWidth}{Must be an odd number. The sliding window size used for smoothing.}
}
\details{
	The region of coverage of surrounding each gene is extracted. If \code{distType} is \code{percent}, then the region is \code{up} \% upstream and \code{down} \% downstream of the wdith of each gene. On the other hand, if \code{distType} is \code{base}, then \code{up} simply refers the number of bases upstream and \code{down} refers to the number of bases downstream to take.  Flanking regions of half the smoothing width are also extracted at each end, so that even the first or last base's coverage of the region can be smoothed. \code{cvgResolution} describes how often to sample from the smoothed coverage. For example, if \code{distType} is \code{percent} and \code{up}, \code{down}, \code{cvgResolution} are all their default values, then the bases sampled will be -50\%, -40\%, -30\%, ..., 90\%, 100\% of the gene width. Genes that are close to the edge of chrosomomes and don't have coverage data at all of the sampling points are excluded from the table.
}
\value{
  A list of tables of coverage. These can then be visualised with the \code{clusteredHeatmap} function.
}

\seealso{\code{\link{clusteredHeatmap}}}
\author{Dario Strbenac}
\examples{
  # Not run.
  # 
  # geneAnno <- read.csv("geneAnno.csv")
  # load("exampleGRangesList.Rdata")
  # tablesList <- getTSScoverage(exampleGRangesList, exptTypes = c("MBD2", "MBD2", "H3K27me3", "H3K27me3"), geneAnno = geneAnno)
}


