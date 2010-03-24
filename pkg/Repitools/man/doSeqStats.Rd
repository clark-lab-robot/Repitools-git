\name{doSeqStats}
\alias{doSeqStats}
\title{Calculate Statistics for Sequencing Data}
\description{
	The fold changes are first calculated over the input regions using the input reads data. These values are then mapped to the regions of the enriched data. Fold changes, p - values, and Z - scores, are then given for each enriched region in a table.
}
\usage{
doSeqStats(reads, whichInputs, whichControl, whichTreat, minCount, blockSize, blocksTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{A \code{GenomeDataList} of the sequencing data.}
  \item{whichInputs}{A vector of length 2, containing indices for \code{reads} that contain the inputs. First element is the index of the control input, and the second index is of the treatment input}
  \item{whichControl}{A vector of indicies for \code{reads} for the enriched controls.}
  \item{whichTreat}{A vector of indices for \code{reads} for the enriched treatments.}
  \item{minCount}{A vector of length 4 that gives the minimum count of an input column for a particular region, minimum total count for both input columns, minumum count for any enriched region in any column, and the minimum total count for any enriched region over all columns.}
  \item{blockSize}{A vector of length 2. The first element is the size of regions along the genome to be used for the input data, covering the whole genome. The second element is the size of such regions to be used for the enriched data.}
    \item{blocksTable}{A data frame with four columns: \code{name}, \code{chr}, \code{start}, \code{end}.}
}
\details{
	All arguments except blocksTable are mandatory. Regions that do not meet any of the cutoffs speicified in minCount are eliminated from further analysis. The input fold changes calculated and grouped together by using the function \code{getCN()}. If \code{blocksTable} is given, only one number needs to be passed as \code{blockSize}. As discribed, the first number is the size of the regions on the input dataset to consider, and it is ensured for short blocks that at least two copy number blocks cover a block, by making the first copy number block start at blockSize elements to the left of the block and the last block can end upto \code{blockSize} bases to the right of the right end of a block.
}
\value{
	A \code{data.frame} of statistics for each of the enriched regions.
}
\author{Dario Strbenac}
\examples{
	# See user guide.
}


