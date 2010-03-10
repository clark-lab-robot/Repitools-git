\name{writeWig}
\alias{writeWig}
\title{Writes sequencing data out into a wiggle file}
\description{Writes sequencing data out into a wiggle file}
\usage{
writeWig(rs.cov, filename, sample = 20, dropZero = TRUE, normalise=1, verbose = TRUE)
}
\arguments{
  \item{rs.cov}{The data to be written, must be a \code{GenomeData} object which \code{coverage} has been called upon}
  \item{filename}{File to write out to - if ending in "gz", a gzipped is written out}
  \item{sample}{At what basepair resolution to sample the genome at}
  \item{dropZero}{Whether to write zero values to the wiggle file - TRUE saves diskspace}
  \item{normalise}{Divide the coverage amount by this number - simple way to normalise to number of reads}
  \item{verbose}{Whether to write progress to screen}
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
function(rs.cov, filename, sample=20, dropZero=TRUE, verbose=TRUE) {
if (class(rs.cov)!="GenomeData" | class(rs.cov[[1]])!="Rle") stop("rs.cov must be a coverage GenomeData object!")
f1 <- if (grepl("gz$",filename)) gzfile(filename, open="wt") else file(filename, open="wt")
scipen <- getOption("scipen")
options(scipen=100)
for (i in 1:length(rs.cov)) {
cat("writing out",names(rs.cov)[i],"\n")
writeLines(paste("variableStep chrom=", names(rs.cov)[i]," span=", sample, sep=""), f1)
bp <- 1:(length(rs.cov[[i]]) \%/\% sample)*sample
bp.score <- as.integer(rs.cov[[i]][bp])
temp <- if (dropZero) paste(bp[bp.score!=0], bp.score[bp.score!=0], sep="\t")
else paste(bp, bp.score, sep="\t")
writeLines(temp, f1)
  }
options(scipen=scipen)
close(f1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
