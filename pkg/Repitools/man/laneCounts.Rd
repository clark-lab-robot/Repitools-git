\name{laneCounts}
\alias{laneCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to count the total number of ChIP-seq reads in each lane.
}
\description{
A function to count the total number of ChIP-seq reads in each lane.}
\usage{
laneCounts(cs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{A \code{GenomeDataList}.
}
}
\details{
}
\value{
	An \code{integer vector}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cs) 
{
    stopifnot("GenomeDataList" \%in\% class(cs))
    require(IRanges)
    require(BSgenome)
    lib.sizes <- IRanges::as.list(BSgenome::gdapply(cs, function(x) length(unlist(x))))
    lib.sizes <- lapply(lib.sizes, IRanges::as.list)
    lib.sizes <- sapply(lib.sizes, function(x) sum(unlist(x)))
    return(lib.sizes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
