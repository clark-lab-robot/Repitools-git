setMethodS3("cpgDensityCalc", "GRanges", function(rs, seqLen = width(rs[1]), ...)
{
	    rs <- resize(rs, seqLen)
	    positionsDF <- data.frame(chr = as.character(seqnames(rs)), position = round((start(rs) + end(rs)) / 2))

            return(cpgDensityCalc(positionsDF, window=seqLen, ...))
})

setMethodS3("cpgDensityCalc", "data.frame", function(locations, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, chunkSize=10000000, ...)
{
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
