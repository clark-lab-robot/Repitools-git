setMethodS3("cpgDensityCalc", "GenomeDataList", function(rs, seqLen, ...) {
	return(lapply(IRanges::as.list(rs), cpgDensityCalc, seqLen, ...))
})


setMethodS3("cpgDensityCalc", "GenomeData", function(rs, seqLen, ...) {
	rs.midpt <- vector(mode='list', length=length(rs))
	names(rs.midpt) <- names(rs)
	
	for (chr in names(rs)) if (length(rs[[chr]][["+"]])+length(rs[[chr]][["+"]])>0) {
		rs.midpt[[chr]] <- data.frame(chr=chr, position=c(rs[[chr]][["+"]]+seqLen/2, rs[[chr]][["-"]]-seqLen/2), stringsAsFactors=FALSE)
	}
	rs.midpt <- do.call(rbind, rs.midpt)
	
	return(cpgDensityCalc(rs.midpt, window=seqLen, ...))
			
})

setMethodS3("cpgDensityCalc", "data.frame", function(locations, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, chunkSize=10000000, ...) {
	wFunction <- match.arg(wFunction)

	if(wFunction == "none") {
		cpgDensity <- sequenceCalc(locations, window, organism, DNAString("CG"), verbose=verbose, chunkSize=chunkSize)
	} else {
		CGfinds <- sequenceCalc(locations, window, organism, DNAString("CG"), verbose=verbose, positions=TRUE)
		distances <- lapply(CGfinds, function(positionsInRegion) {if(!is.null(positionsInRegion)) abs(positionsInRegion)})
		if(wFunction == "linear") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(1 - (distancesInRegion / (window / 2))))
		} else if(wFunction == "log") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(log2(2 - (distancesInRegion / (window / 2)))))
		} else { # Exponential decay was requested.
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(exp(-5 * distancesInRegion / (window / 2))))	
		}
		rm(CGfinds)
	}
	gc()
	return(cpgDensity)
})
