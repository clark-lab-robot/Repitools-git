\name{processAffyAnnotation}
\alias{processAffyAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A simple routine to read and process Affy annotation files (It may not be very general, so far just tested on Gene 1.0 ST transcript.csv files)}
\description{
Reads the CSV file and only returns a subset of the columns and has an option to peel off the RefSeq identifier.
}
\usage{
processAffyAnnotation(csvFile, skip=19, keepColumns=c("probeset_id","seqname","strand","start","stop","total_probes","category"), getRefseq=FALSE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csvFile}{filename of Affymetrix CSV file}
  \item{skip}{number of lines at top of CSV file to skip}
  \item{keepColumns}{a vector of column names to keep from the original CSV file}
  \item{getRefseq}{logical, whether to process the 'gene_assignment' column for }
  \item{verbose}{logical, whether to write progress to the screen}
}
\value{
a data.frame with a row for each Affymetrix probeset id.
}
\author{Mark Robinson}
\examples{

## not run
# chrs <- paste( "chr", c(1:22,"X","Y","M"), sep="" )
# names(chrs) <- 1:25
# paa <- processAffyAnnotation("../hugene/HuGene-1_0-st-v1.na27.hg18.transcript.csv")
# paa <- paa[ paa$category=="main" & paa$type %in% c("RefSeq","GenBank") & paa$seqname %in% chrs, ]

}
