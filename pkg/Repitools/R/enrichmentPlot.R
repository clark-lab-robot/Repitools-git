enrichmentPlot <- function(rs, organism, seqLen, cols=rainbow(length(rs)), xlim=c(0,20), main="Enrichment Plot", total.lib.size=TRUE, verbose=FALSE, ...) {
	if (length(cols)!=length(rs)) stop("rs and cols must have the same number of elements.")
	if (verbose) cat("Extending reads\n")
	rs.ext <- endoapply(rs, resize, seqLen)
	if (verbose) cat("Creating coverage object\n")
	rs.cov <- IRanges::lapply(rs.ext, coverage)
	if (verbose) cat("Calculating enrichment\n")
	rs.enrich <- enrichmentCalc(rs.cov, organism, ...)
	if (total.lib.size) {
		if (verbose) cat("Normalising to reads per lane\n")
		rs.counts <- elementLengths(rs)
		for (i in 1:length(rs)) rs.enrich[[i]]$coverage <- rs.enrich[[i]]$coverage/(rs.counts[[i]]/1000000)
	}
	plot(x=rs.enrich[[1]]$coverage, y=rs.enrich[[1]]$bases, type="l", col=cols[1], xlim=xlim, main=main, ylab="Frequency", log="y",
			xlab=if (total.lib.size) "Normalised Enrichment Level of reads" else "Enrichment Level of reads")
	if (length(rs)>1) for (i in 2:length(rs)) {
		lines(x=rs.enrich[[i]]$coverage, y=rs.enrich[[i]]$bases, col=cols[i])
	}
	legend("topright", lty=1, col=cols, legend=names(rs))
	invisible(rs.enrich)
}
