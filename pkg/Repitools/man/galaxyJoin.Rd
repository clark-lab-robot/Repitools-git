\name{galaxyJoin}
\alias{galaxyJoin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper for the python script to merge ranges}
\description{
Takes 2 tables of coordinates (or 2 appropriately formatted text files), calls the python script to match up the coordinates, grabs the output and returns it as an R data.frame.
}
\usage{
galaxyJoin(tabOrFile1 = NULL, tabOrFile2 = NULL, cols1 = 1:3, cols2 = 1:3, mincols = 1, fill = c("none", "right", "left", "both"), delete = TRUE, script = "gops_join.py")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tabOrFile1}{first table of coordinates, or file}
  \item{tabOrFile2}{second table of coordinates, or file}
  \item{cols1}{the 3 columns in first table to use -- chromosome, start, end}
  \item{cols2}{the 3 columns in second table to use -- chromosome, start, end}
  \item{mincols}{amount of overlap required}
  \item{fill}{how to return the output (see below)}
  \item{delete}{logical, whether to delete files afterwards}
  \item{script}{the name of the perl script to call}
}
\details{

Some more details of the command line usage of the python:

usage: %prog bed_file_1 bed_file_2 out_file
   -1, --cols1=N,N,N,N: Columns for start, end, strand in first file
   -2, --cols2=N,N,N,N: Columns for start, end, strand in second file
   -m, --mincols=N: Require this much overlap (default 1bp)
   -f, --fill=N: none, right, left, both

}
\value{
Returns a data.frame of the overlap, according to what is specified in the \code{fill} argument
}
\author{Mark Robinson}
\examples{

## tab1 <- data.frame(chr="chr1",start=seq(0,10000,length=11),end=seq(500,10500,length=11))
## tab2 <- data.frame(chr="chr1",start=seq(250,10250,length=51),end=seq(300,10300,length=51))
## options(galaxyPath="~/software/galaxy_dist/tools/new_operations/",scipen=3)
## gj <- galaxyJoin(tab1, tab2)
## gj <- galaxyJoin(tab1, tab2, fill="both")

}
