gcContentCalc <- function(locations, window = 500, organism, verbose=FALSE, chunkSize=10000000) {	
	return(sequenceCalc(locations, window, organism, DNAString("S"), fixed=FALSE, Nmask=TRUE, verbose=verbose, chunkSize=chunkSize)/window)
}
