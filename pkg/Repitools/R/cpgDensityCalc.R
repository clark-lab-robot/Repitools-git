setGeneric("cpgDensityCalc", signature = "x", function(x, ...){standardGeneric("cpgDensityCalc")})

setMethod("cpgDensityCalc", "GenomeDataList", function(x, seqLen, ...) {
            return(lapply(IRanges::as.list(x), cpgDensityCalc, seqLen, ...))
        })


setMethod("cpgDensityCalc", "GenomeData", function(x, seqLen, window = seqLen / 2, ...) {
            rs.midpt <- vector(mode='list', length=length(x))
            names(rs.midpt) <- names(x)
            
            for (chr in names(x)) if (length(x[[chr]][["+"]])+length(x[[chr]][["+"]])>0) {
                    rs.midpt[[chr]] <- data.frame(chr=chr, position=c(x[[chr]][["+"]]+seqLen/2, x[[chr]][["-"]]-seqLen/2), stringsAsFactors=FALSE)
                }
            rs.midpt <- do.call(rbind, rs.midpt)
            
            return(cpgDensityCalc(rs.midpt, ...))
            
        })

setMethod("cpgDensityCalc", "GRanges", function(x, seqLen, window = seqLen / 2, ...)
{
	    x <- resize(x, seqLen)
	    positionsDF <- data.frame(chr = as.character(seqnames(x)), position = round((start(x) + end(x)) / 2))

            return(cpgDensityCalc(positionsDF, ...))
})

setMethod("cpgDensityCalc", "data.frame", function(x, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, chunkSize=10000000, ...)
{
            wFunction <- match.arg(wFunction)
            
            if(wFunction == "none") {
                cpgDensity <- sequenceCalc(x, window, organism, DNAString("CG"), verbose=verbose, chunkSize=chunkSize)
            } else {
                CGfinds <- sequenceCalc(x, window, organism, DNAString("CG"), verbose=verbose, positions=TRUE)
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
