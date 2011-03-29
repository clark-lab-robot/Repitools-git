setGeneric("cpgDensityPlot", function(x, ...){standardGeneric("cpgDensityPlot")})

setMethod("cpgDensityPlot", "GenomeDataList", function(x, ...) {
    cpgDensityPlot(.GDL2GRL(x), ...)
})

setMethod("cpgDensityPlot", "GRangesList", function(x, cols=rainbow(length(x)), xlim=c(0,20), lty = 1, lwd = 1, main="CpG Density Plot", verbose=TRUE, ...) {
	if (length(cols)!=length(x)) stop("x and cols must have the same number of elements.")
    if (verbose) message("Calculating CpG density")
	x.cpg <- cpgDensityCalc(x, verbose=verbose, ...)
	x.den <- lapply(x.cpg, density)
	ymax <- max(sapply(x.den, function(u) max(u$y)))
	plot(x=x.den[[1]]$x, y=x.den[[1]]$y, type="l", col=cols[1], xlim=xlim, ylim=c(0,ymax), main=main, ylab="Frequency", xlab="CpG Density of reads", lty = lty, lwd = lwd)
	if (length(x)>1) for (i in 2:length(x)) {
		lines(x=x.den[[i]]$x, y=x.den[[i]]$y, col=cols[i], lty = lty, lwd = lwd)
	}
	legend("topright", col=cols, legend=names(x), lty = lty, lwd = lwd)
	invisible(x.cpg)
})
