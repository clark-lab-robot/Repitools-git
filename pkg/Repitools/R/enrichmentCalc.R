enrichmentCalc <- function(rs, organism, chrs, seqLen=NULL) {
	require(GenomicRanges)

	if (class(rs) != "GRanges") {
		stop("rs must be a GRanges object.")

	seqlengths(rs) <- seqlengths(organism)[chrs]
	coverage <- colSums(table(coverage(resize(rs, seqLen))))
	coverageTable <- cbind(as.numeric(names(coverage)), coverage)
	colnames(coverageTable) <- c("coverage", "bases")

	return(coverageTable)
}
