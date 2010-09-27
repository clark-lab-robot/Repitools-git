\name{genQC}
\alias{genQC}
\title{Generate Quality Checking Information For Sequencing Data}
\description{The function generates a number of tables of data. See below. In addition to the returned tables, a graphical summary of the data is output as a PDF that is called "QC <exptName> plots.pdf"}
\usage{
genQC(paths, exptName, laneNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{A vector of character strings containing absolute paths to each lane in the experiment.}
  \item{exptName}{The name of the experiment which the lanes are about.}
  \item{laneNames}{Name of the sample in each sequencing lane.}

}
\value{
	A list as long as the number of lanes. Each element of the list contains a list with elements :
	\item{avgQual}{Average quality scores, per base.}
	\item{freqAll}{Frequency of bases A, G, C, T, N over the length of the reads.}
	\item{freqAligned}{Frequency of bases A, G, C, T, N over the length of only the uniquely aligned reads.}
	\item{mmTable}{Frequencies of mismatches at each position.}
}
\author{Dario Strbenac}
\examples{
	# Not given.
}


