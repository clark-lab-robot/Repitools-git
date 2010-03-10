\name{annotationCounts}
\alias{annotationCounts}
\title{Counts the number of sequencing reads surrounding supplied annotations}
\description{
Counts the number of sequencing reads surrounding supplied annotations
}
\usage{
annotationCounts(rs, annotation, bpUp, bpDown, seqLen = NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{A \code{GenomeData} or \code{GenomeDataList} object}
  \item{annotation}{a \code{data.frame} with required elements \code{chr} (character),  \code{start} (integer), \code{end} (integer) and optional elements \code{strand} ("+" or "-") and \code{name} (character)}
  \item{bpUp}{The number of bases upstream to look}
  \item{bpDown}{The number of bases downstream to look}
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
function (rs, annotation, bpUp, bpDown, seqLen = NULL) 
{
    if (class(rs) == "GenomeData") 
        rs <- GenomeDataList(list(rs))
    anno.ranges <- IRanges(start = ifelse(annotation$strand == 
        "+", annotation$position - bpUp, annotation$position - 
        bpDown), end = ifelse(annotation$strand == "+", annotation$position + 
        bpDown, annotation$position + bpUp))
    anno <- data.frame(chr = annotation$chr, start = ifelse(annotation$strand == 
        "+", annotation$position - bpUp, annotation$position - 
        bpDown), end = ifelse(annotation$strand == "+", annotation$position + 
        bpDown, annotation$position + bpUp), name = annotation$name)
    annotationBlocksCounts(rs, anno, seqLen)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
