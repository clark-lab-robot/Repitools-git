\name{checkProbes}
\alias{checkProbes}
\title{Check Probe Specificity for Some Regions}
\description{Given a table of regions, and a table of probe mappings to the genome, a plot is created across every region of how many probe mapped to each position.
}
\usage{
checkProbes(regionsTable, probesTable, outFile = "Probe Specificity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regionsTable}{A \code{data.frame} with (at least) columns \code{chr}, \code{start}, \code{end}, and \code{name}.  The start and end of a region describe the window analysed, and the position describes the location of the TSS.}
  \item{probesTable}{A \code{data.frame} describing where the probes hit with (at least) columns \code{name} (identifier of a probe), \code{chr}, \code{start}, and \code{end}.}
  \item{outFile}{Name of output PDF, without the .pdf extension.}
}
\value{
  A PDF of plots is created in the current directory, one for each region in the table of regions. The dots in the plot show where a probe hits (the x - axis) and how many places in total the probe hits in the genome (y -axis).
}
\author{Dario Strbenac}
\examples{
	pTable <- data.frame(name = c("probe A", "probe B", "probe C", "probe C", "probe C"), strand = c('+', '-', '+', '-', '-'), chr = c("chr1", "chr2", "chr1", "chr2", "chr2"), start = c(20, 276, 101, 101, 151), end = c(44, 300, 125, 125, 175))
	rTable <- data.frame(name = c("gene 1", "gene 2", "gene 3"), chr = c("chr1", "chr2", "chr2"), strand = c('+', '-', '+'), start = c(20, 50, 75), end = c(200, 300, 400))
	rTable$pos <- ifelse(rTable$strand == '+', rTable$start, rTable$end)
	checkProbes(rTable, pTable)
}

