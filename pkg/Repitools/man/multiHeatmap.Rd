\name{multiHeatmap}
\alias{multiHeatmap}
\title{Superfigure plots}
\description{
This function takes a list of matrices and plots heatmaps for each one.  There are several features for the spacing (X and Y), colour scales, titles and label sizes.  If a matrix has row and/or column names, these are added to the plot.
}
\usage{
multiHeatmap(dataList, colourList, titles = NULL, main = "", showColour = TRUE, xspace = 1, cwidth = 0.5, ystarts = c(0.05, 0.9, 0.925, 0.95, 0.98), rlabelcex = 1, clabelcex = 1, titlecex = 1.2, maincex = 1.5, scalecex = 0.7, offset=.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataList}{A \code{list} of matrices to be plotted as different panels}
  \item{colourList}{A \code{list} of colourscales (if length 1, it is copied for all panels of the plot)}
  \item{titles}{A vector of panel titles}
  \item{main}{A main title}
  \item{showColour}{logical or logical vector, whether to plot the colour scale}
  \item{xspace}{The space between the panels (relative to number of columns).  This can be either a scalar or a vector of \code{length(dataList)+1} }
  \item{cwidth}{widths of the colour scales relative to the width of the panels}
  \item{ystarts}{A vector of length 5 of numbers between 0 and 1 giving the relative Y positions of where the heatmaps, colourscale labels, colour scales, panel titles and main title (respectively) start}
  \item{rlabelcex}{character expansion factor for row labels}
  \item{clabelcex}{character expansion factor for column labels}
  \item{titlecex}{character expansion factor for panel titles}
  \item{maincex}{character expansion factor for main title}
  \item{scalecex}{character expansion factor for colour scale labels}
  \item{offset}{small offset to adjust scales for point beyond the colour scale boundaries}
}
\value{
This function is called for its output, a plot in the current device.
}
\author{Mark Robinson}
\examples{
library(gplots)

cL <- NULL
br <- seq(-3,3,length=101)
col <- colorpanel(low="blue",mid="grey",high="red",n=100)
cL[[1]] <- list(breaks=br,colors=col)
br <- seq(-2,2,length=101)
col <- colorpanel(low="green",mid="black",high="red",n=100)
cL[[2]] <- list(breaks=br,colors=col)
br <- seq(0,20,length=101)
col <- colorpanel(low="black",mid="grey",high="white",n=100)
cL[[3]] <- list(breaks=br,colors=col)

testD <- list(matrix(runif(400),nrow=20),matrix(rnorm(100),nrow=20),matrix(rpois(100,lambda=10),nrow=20))
colnames(testD[[1]]) <- letters[1:20]
rownames(testD[[1]]) <- paste("row",1:20,sep="")

multiHeatmap(testD,cL,xspace=1)
}
