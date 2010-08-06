\name{significancePlots}
\alias{significancePlots}
\alias{significancePlots.GenomeDataList}
\alias{significancePlots.AffymetrixCelSet}
\alias{significancePlots.matrix}
\title{Create line plots of averaged signal across a promoter compared to random sampling}
\description{
Creates a plot where the average signal across a promoter of supplied gene lists is compared to random samplings of all genes, with a shaded confidence area.
}
\usage{
\method{significancePlots}{GRangesList}(rs, coordinatesTable, design=NULL, upStream=7500, downStream=2500, by=100, bw=300, total.lib.size=TRUE, seqLen=NULL, verbose=FALSE, ...)
\method{significancePlots}{AffymetrixCelSet}(cs, probeMap=NULL, coordinatesTable=NULL, upStream=7500, downStream=2500, by=100, bw=300, log2adjust=TRUE, verbose=FALSE, ...)
\method{significancePlots}{matrix}(dataMatrix, lookupTable, geneList, titles, nSamples=1000, confidence=0.975, legend.plot="topleft", cols=rainbow(length(geneList)), removeZeros=TRUE, useMean=FALSE, ...)
}
\arguments{
  \item{dataMatrix}{\code{matrix} of data to use}
  \item{rs}{\code{GRangesList} of sequencing data to use}
  \item{probeMap}{A pre-made mapping of probes to a coordinatesTable}
  \item{cs}{\code{AffymetrixCelSet} of array data to use}
  \item{lookupTable}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{design}{Optional for GDL data - design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{geneList}{Names \code{list} of \code{logical} or \code{integer} vectors, specifying the genes to be averaged and plotted}
  \item{titles}{\code{character} vector of titles for the plots}
  \item{nSamples}{number of times to randomly sample from all genes}
  \item{confidence}{percentage confidence interval to be plotted (must be > 0.5 and < 1.0}
  \item{legend.plot}{Where to plot the legend - directly passed to \code{legend}. NA suppresses the legend}
  \item{cols}{the colour for each of the genelists supplied}
  \item{removeZeros}{Remove zero values.}
  \item{useMean}{Whether to use mean instead of median within gene lists (useful for seq data).}  
  \item{coordinatesTable}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{position}, \code{strand}}
  \item{upStream}{number of bases upstream to plot}
  \item{downStream}{number of bases upstream to plot}
  \item{by}{how many bases between data points}
  \item{bw}{how many bases up and downstream to be counted}
  \item{total.lib.size}{whether to use the total number of sequencing reads per lane to normalise.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{verbose}{Print details of processing.}
  \item{...}{arguments passed from the \code{GenomeDataList} method into the \code{matrix} method.}
}
\details{
For each column of data in \code{dataMatrix} a plot is created showing the average signal of the genes specified in each list element of \code{geneList} compared to \code{nSamples} random samplings of all genes, with \code{confidence} % intervals shaded. If an element of \code{geneList} is a \code{logical} vector, its length must be the same as the number of rows in \code{lookupTable}.
}
\value{
  A plot.
}
\author{Aaron Statham}
\examples{
  # See examples in manual.
}
