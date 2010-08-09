enrichmentCalc <- function(rs, organism, seqLen=NULL) {
	require(GenomicRanges)

	if (class(rs) != "GRanges")
		stop("rs must be a GRanges object.")
	chrs <- levels(seqnames(rs))
	seqlengths(rs) <- seqlengths(organism)[chrs]
	coverage <- colSums(table(coverage(resize(rs, seqLen))))
	coverageTable <- cbind(as.numeric(names(coverage)), coverage)
	colnames(coverageTable) <- c("coverage", "bases")

	return(coverageTable)
}
