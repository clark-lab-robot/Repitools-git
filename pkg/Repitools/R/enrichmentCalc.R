setGeneric("enrichmentCalc", function(rs, organism, ...){standardGeneric("enrichmentCalc")})

setMethod("enrichmentCalc", c("GenomeDataList", "BSgenome"), function(rs, organism, seqLen=NULL, ...) {
	return(lapply(IRanges::as.list(rs), enrichmentCalc, organism, seqLen, ...))
})

setMethod("enrichmentCalc", c("GenomeData", "BSgenome"), function(rs, organism, seqLen=NULL, do.warn=FALSE) {
	if (class(rs) != "GenomeData") {
	} else if (class(rs[[1]]) != "Rle") {
		if (class(rs[[1]]) != "IRanges") {
			if (is.null(seqLen)) 
				stop("If rs has not been extended, seqLen must be supplied")
			rs <- chipseq::extendReads(rs, seqLen)
		}
		rs <- BSgenome::gdapply(rs, IRanges::coverage)
	}
	rs <- IRanges::as.list(rs)
	max.cov <- max(sapply(rs, max))
	cov.table <- numeric(max.cov+1)
	names(cov.table) <- 0:max.cov
	chr.lengths <- GenomicRanges::seqlengths(organism)
	for (i in 1:length(rs)) {
		temp <- IRanges::table(rs[[i]])
		cov.table[names(temp)] <- cov.table[names(temp)]+temp
		if (length(rs[[i]])>chr.lengths[names(rs)[i]]) cov.table["0"] <- cov.table["0"] + length(rs[[i]])-chr.lengths[names(rs)[i]]
			else if (do.warn) warning("Coverage of ",names(rs)[i]," extends off end")
	}
	cov.table <- data.frame(coverage=0:max.cov, bases=cov.table, row.names=NULL)
	return(cov.table)
})

setMethod("enrichmentCalc", c("GRanges", "BSgenome"), function(rs, organism, seqLen=NULL) {
	require(GenomicRanges)
	
	rs <- resize(rs, seqLen)
	chrs <- levels(seqnames(rs))
	seqlengths(rs) <- seqlengths(organism)[chrs]
	coverage <- colSums(table(coverage(resize(rs, seqLen))))
	coverageTable <- data.frame(as.numeric(names(coverage)), coverage)
	colnames(coverageTable) <- c("coverage", "bases")

	return(coverageTable)
})
