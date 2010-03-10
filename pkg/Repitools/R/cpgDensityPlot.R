setMethodS3("cpgDensityPlot", "GenomeDataList", function(rs, cols=rainbow(length(rs)), xlim=c(0,20), main="CpG Density Plot", ...) {
	if (length(cols)!=length(rs)) stop("rs and cols must have the same number of elements.")
	rs.cpg <- cpgDensityCalc(rs, ...)
	rs.den <- lapply(rs.cpg, density)
	ymax <- max(sapply(rs.den, function(u) max(u$y)))
	plot(x=rs.den[[1]]$x, y=rs.den[[1]]$y, type="l", col=cols[1], xlim=xlim, ylim=c(0,ymax), main=main, ylab="Frequency", xlab="CpG Density of reads")
	if (length(rs)>1) for (i in 2:length(rs)) {
		lines(x=rs.den[[i]]$x, y=rs.den[[i]]$y, col=cols[i])
	}
	legend("topright", lty=1, col=cols, legend=names(rs))
	invisible(rs.cpg)
})
