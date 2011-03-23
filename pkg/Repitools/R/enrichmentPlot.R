setGeneric("enrichmentPlot", function(x, ...) standardGeneric("enrichmentPlot"))

setMethod("enrichmentPlot", "GenomeDataList",
    function(x, ...)
{
    enrichmentPlot(.GDL2GRL(x), ...)
})

setMethod("enrichmentPlot", "GRangesList",
    function(x, organism, seq.len, cols=rainbow(length(x)), xlim=c(0,20), main="Enrichment Plot", total.lib.size=TRUE, verbose=TRUE, ...)
{
    if (length(cols)!=length(x)) stop("x and cols must have the same number of elements.")
    if (verbose) message("Calculating enrichment")
    x.enrich <- enrichmentCalc(x, organism, seq.len, verbose)
    if (total.lib.size) {
	    if (verbose) message("Normalising to reads per lane")
	    x.counts <- elementLengths(x)
	    for (i in 1:length(x)) x.enrich[[i]]$coverage <- x.enrich[[i]]$coverage/(x.counts[[i]]/1000000)
    }
    plot(x=x.enrich[[1]]$coverage, y=x.enrich[[1]]$bases, type="l", col=cols[1], xlim=xlim, main=main, ylab="Frequency", log="y",
		    xlab=if (total.lib.size) "Normalised Enrichment Level of reads" else "Enrichment Level of reads", ...)
    if (length(x)>1) for (i in 2:length(x)) {
	    lines(x=x.enrich[[i]]$coverage, y=x.enrich[[i]]$bases, col=cols[i], ...)
    }
    legend("topright", lty=1, col=cols, legend=names(x), ...)
    invisible(x.enrich)
})
