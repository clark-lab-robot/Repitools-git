\name{correlationGraphs}
\alias{correlationGraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correlationGraphs(intensities, lookup, compareExpression, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intensities}{ ~~Describe \code{intensities} here~~ }
  \item{lookup}{ ~~Describe \code{lookup} here~~ }
  \item{compareExpression}{ ~~Describe \code{compareExpression} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (intensities, lookup, compareExpression, ...) 
{
    corScores = matrix(NA, nrow = ncol(intensities), ncol = length(lookup), 
        dimnames = list(colnames(intensities), names(lookup)))
    for (i in 1:ncol(intensities)) corScores[i, ] = sapply(lookup, 
        scoreCorrelation, intensities[, i], compareExpression, 
        minProbes = 2)
    cols = rainbow(nrow(corScores))
    plot(0, type = "n", xlim = c(-7500, 2500), ylim = c(-1, 1), 
        xlab = "Distance from TSS", ylab = "Correlation with Expression", 
        ...)
    for (i in 1:nrow(corScores)) lines(as.integer(colnames(corScores)), 
        corScores[i, ], col = cols[i])
    legend("topleft", legend = rownames(corScores), col = cols, 
        lty = 1)
  }
}
