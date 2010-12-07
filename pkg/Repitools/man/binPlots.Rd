\name{binPlots}
\alias{binPlots}
\alias{binPlots,GenomeDataList-method}
\alias{binPlots,GRangesList-method}
\alias{binPlots,AffymetrixCelSet-method}
\alias{binPlots,matrix-method}
\title{Create line plots of averaged signal across a promoter}
\description{
Using a specified ordering of genes, they are split into multiple bins.  In each bin, the signal across is summarized and displayed visually.  The plots are created in an array, one row for each column of the data matrix and one column for each ordering specified.
}
\usage{
  \S4method{binPlots}{GenomeDataList}(x, coordinatesTable, design=NULL, upStream=7500, downStream=2500, by=100, bw=300, libSize="lane", seqLen=NULL, verbose=FALSE, Acutoff=NULL, ...)
  \S4method{binPlots}{GRangesList}(x, coordinatesTable, design=NULL, upStream=7500, downStream=2500, by=100, bw=300, libSize="lane", seqLen=NULL, verbose=FALSE, Acutoff=NULL, ...)
  \S4method{binPlots}{AffymetrixCelSet}(x, probeMap=NULL, coordinatesTable=NULL, upStream=7500, downStream=2500, by=100, bw=300, log2adjust=TRUE, verbose=FALSE, ...)
  \S4method{binPlots}{matrix}(x, lookupTable, ordering, ordLabel, plotType = c("line", "heatmap", "terrain", "boxplot"), nbins = 10, cols = NULL, lwd = 3, lty = 1, sameScale = TRUE, symmScale = FALSE, verbose = FALSE, removeZeros=TRUE, useMean=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{matrix}, \code{GenomeDataList} or \code{GRangesList} object, or \code{AffymetrixCelSet} object of data to use}
  \item{probeMap}{A pre-made mapping of probes to a coordinatesTable}
  \item{coordinatesTable}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{name}, \code{start}, \code{end}, \code{strand}.}
  \item{lookupTable}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{design}{Optional for GDL data - design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{ordering}{\code{data.frame} of either numeric or factor variables, with the same number of rows as \code{lookupTable}}
  \item{ordLabel}{Character string that describes what type of data the ordering is. e.g. "log2 expression". Used to label relevant plot axis.}
  \item{plotType}{type of plot}
  \item{nbins}{number of bins}
  \item{cols}{vector of colours}
  \item{lwd}{line width (either scalar or vector)}
  \item{lty}{line type (either scalar or vector)}
  \item{sameScale}{should the scale on all plots be the same}
  \item{symmScale}{should the scale on plots be symmetrical ie centred around 0}
  \item{verbose}{Print details of processing.}
  \item{removeZeros}{Remove zero values.}
  \item{useMean}{Whether to use mean instead of median within bins (useful for seq data).}
  \item{upStream}{number of bases upstream to plot}
  \item{downStream}{number of bases upstream to plot}
  \item{by}{how many bases between data points}
  \item{bw}{how many bases up and downstream to be counted}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{libSize}{Either the string \code{"ref"} to normalise to a reference or \code{"lane"} to multiply each library to be the same size.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{Acutoff}{If \code{libSize} is \code{"ref"}, this argument must be provided. Otherwise, it must not. Cutoff on the "A" values to take, before calculating trimmed mean.}
  \item{...}{further arguments to send to \code{bxp}}
}
\details{
If \code{plotType="line"}, a line is plotted for each bin across the promoter.
  
If \code{plotType="heatmap"}, a series of bins are plotted as a heatmap.  This can be useful to display a larger number of bins.

If \code{plotType="terrain"}, a series of bins are plotted as a 3d-terrain map.  This can be useful to display a larger number of bins.

The option \code{plotType="boxplot"} is really only suitable for 1- or 2-group plots.  In this case, instead of plotting the summarized data, it plots the whole distribution of signal across the whole promoter.
}
\value{
  Either a single- or multiple-panel figure.
}
\author{Mark Robinson}
\examples{
  # See examples in manual.
}
