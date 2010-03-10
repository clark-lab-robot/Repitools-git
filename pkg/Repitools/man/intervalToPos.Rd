\name{intervalToPos}
\alias{intervalToPos}
\title{Add a Position Column to a Data Frame}
\description{
A convenience method to take a \code{data.frame} that has (at least) \code{start} and \code{end} columns, and adds a \code{position} column.
}
\usage{
intervalToPos(regions)
}
\arguments{
  \item{regions}{A data frame with (at least) \code{start} and \code{end} columns.}
  }
\details{
	For each genomic region, the position is calculated the midpoint of the start and end, and rounded (to the nearest even number) where this does not give an integer.
}
\value{
The same \code{data.frame} as was passed in, with the additional column \code{position}.
}
\author{Dario Strbenac}
\examples{
	rdf <- data.frame(chr=c(1, 2, 20), start=c(100, 200, 300), end=c(400,500,600), strand=c('+', '-', '+'))
	intervalToPos(rdf)
}


