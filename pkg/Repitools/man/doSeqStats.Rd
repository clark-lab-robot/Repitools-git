\name{doSeqStats}
\alias{doSeqStats}
\alias{doSeqStats,GenomeDataList-method}
\alias{doSeqStats,GRangesList-method}
\title{Calculate Statistics for Sequencing Data}
\description{
	The fold changes are first calculated over the input regions using the input reads data. These values are then mapped to the regions of the enriched data. Fold changes, p - values, and Z - scores, are then given for each enriched region in a table.
}
\usage{
  \S4method{doSeqStats}{GenomeDataList}(reads, seqLen, whichInputs, whichControl, whichTreat, minCount, blockSize, CNlevels, blocksTable, bpUp, bpDown)
  \S4method{doSeqStats}{GRangesList}(reads, seqLen, whichInputs, whichControl, whichTreat, minCount, blockSize, CNlevels, blocksTable, bpUp, bpDown)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{A \code{GenomeDataList} or \code{GRangesList} of the sequencing data.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used.}
  \item{whichInputs}{A vector of length 2, containing indices for \code{reads} that contain the inputs. First element is the index of the control input, and the second index is of the treatment input.}
  \item{whichControl}{A vector of indicies for \code{reads} for the enriched controls.}
  \item{whichTreat}{A vector of indices for \code{reads} for the enriched treatments.}
  \item{minCount}{A vector of length 2 that gives the minimum total count for both input columns and the minimum total count for any enriched region over all columns.}
  \item{blockSize}{A vector of length 2. The first element is the size of regions along the genome to be used for the input data, covering the whole genome. The second element is the size of such regions to be used for the enriched data. If \code{blocksTable} is given, then just a single number specifying the size of regions along the input data.}
  \item{CNlevels}{An integer. How many discrete copy number levels to make from the input data.}
  \item{blocksTable}{A data frame with four columns: \code{name}, \code{chr}, \code{start}, \code{end}.}
  \item{bpUp}{Only used if blocksTable provided. This many bases upstream of the TSS. If strand information is not present in \code{blocksTable}, then the windows are created this many bases upstream and downstream of the centre of the blocks.}
  \item{bpDown}{Only used if blocksTable provided and contains strand information. The windows extend this many bases downstream of the TSS.}
}
\details{
	All arguments except \code{blocksTable} are mandatory. Regions that do not meet any of the cutoffs speicified in minCount are eliminated from further analysis. The input fold changes calculated and grouped together by using the function \code{getCN()}. If \code{blocksTable} is given, only one number needs to be passed as \code{blockSize}. As described, the first number is the size of the regions on the input dataset to consider, and it is ensured for short blocks that at least two copy number blocks cover a block, by making the first copy number block start at \code{blockSize} elements to the left of the block and the last block can end upto \code{blockSize} bases to the right of the right end of a block.
}
\value{
	A \code{data.frame} of statistics for each of the enriched regions.
}
\author{Dario Strbenac}
\examples{
	# See user guide.
}


