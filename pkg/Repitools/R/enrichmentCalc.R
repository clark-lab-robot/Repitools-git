setMethodS3("enrichmentCalc", "GenomeDataList", function(rs, organism, seqLen, ...) {
	return(lapply(IRanges::as.list(rs), enrichmentCalc, organism, seqLen, ...))
})


setMethodS3("enrichmentCalc", "GenomeData", function(rs, organism, seqLen=NULL, do.warn=FALSE, ...) {
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
	chr.lengths <- BSgenome::seqlengths(organism)
	for (i in 1:length(rs)) {
		temp <- IRanges::table(rs[[i]])
		cov.table[names(temp)] <- cov.table[names(temp)]+temp
		if (length(rs[[i]])>chr.lengths[names(rs)[i]]) cov.table["0"] <- cov.table["0"] + length(rs[[i]])-chr.lengths[names(rs)[i]]
			else if (do.warn) warning("Coverage of ",names(rs)[i]," extends off end")
	}
	cov.table <- data.frame(coverage=0:max.cov, bases=cov.table, row.names=NULL)
	return(cov.table)
})
