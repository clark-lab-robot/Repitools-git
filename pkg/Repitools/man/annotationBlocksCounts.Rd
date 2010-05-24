\name{annotationBlocksCounts}
\alias{annotationBlocksCounts}
\title{Counts the number of sequencing reads within supplied genomic blocks}
\description{
Counts the number of sequencing reads within supplied genomic blocks
}
\usage{
annotationBlocksCounts(rs, annotation, seqLen = NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeData} or \code{GenomeDataList} object}
  \item{annotation}{a RangedData object (preferable) or a \code{data.frame} with required elements \code{chr} (character),  \code{start} (integer), \code{end} (integer) and optional elements \code{strand} ("+" or "-") and \code{name} (character)}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{verbose}{logical, whether to print progress}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (rs, annotation, seqLen = NULL) 
{
    if (class(rs) == "GenomeData") 
        rs <- GenomeDataList(list(rs))
    anno.counts <- matrix(as.integer(NA), nrow = nrow(annotation), 
        ncol = length(rs), dimnames = list(annotation$name, names(rs)))
    anno.ranges <- IRanges(start = annotation$start, end = annotation$end)
    for (i in 1:length(rs)) {
        if (!class(rs[[i]][[1]]) == "IRanges") {
            if (is.null(seqLen)) 
                stop("If rs has not been extended, seqLen must be supplied")
            rs[[i]] <- extendReads(rs[[i]], seqLen = seqLen)
        }
        for (chr in unique(annotation$chr)) {
            which.anno <- annotation$chr == chr
            if (is.null(rs[[i]][[chr]])) 
                anno.counts[which.anno] <- 0
            else anno.counts[which.anno, i] <- as.table(findOverlaps(anno.ranges[which.anno], 
                rs[[i]][[chr]]))
        }
    }
    anno.counts
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
