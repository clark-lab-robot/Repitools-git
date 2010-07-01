setMethodS3("cpgDensityCalc", "GenomeDataList", function(rs, seqLen, ...) {
            return(lapply(IRanges::as.list(rs), cpgDensityCalc, seqLen, ...))
        })


setMethodS3("cpgDensityCalc", "GenomeData", function(rs, seqLen, ...) {
            require(chipseq)

	    rs <- extendReads(rs, seqLen = seqLen)
	    rs <- as(rs, "RangedData")
            return(cpgDensityCalc(rs, window=seqLen, ...))
        })

setMethodS3("cpgDensityCalc", "RangedData", function(locations, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, chunkSize=10000000, ...) {
            wFunction <- match.arg(wFunction)
            if(wFunction == "none") {
                cpgDensity <- sequenceCalc(locations, organism, DNAString("CG"), verbose=verbose, chunkSize=chunkSize)
            } else {
                CGfinds <- sequenceCalc(locations, organism, DNAString("CG"), verbose=verbose, positions=TRUE)
                distances <- lapply(CGfinds, function(positionsInRegion) {if(!is.null(positionsInRegion)) abs(positionsInRegion)})
                if(wFunction == "linear") {
                    cpgDensity <- mapply(function(dists, aSpan) sum(1 - (dists / (aSpan / 2))), distances, width(locations))
                } else if(wFunction == "log") {
                    cpgDensity <- mapply(function(dists, aSpan) sum(log2(2 - (dists / (aSpan / 2)))), distances, width(locations))
                } else { # Exponential decay was requested.
                    cpgDensity <- mapply(function(dists, aSpan) sum(exp(-5 * dists / (aSpan / 2))), distances, width(locations))	
                }
                rm(CGfinds)
            }
            gc()
            return(cpgDensity)
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
