mappabilityCalc <- function(locations, window = 500, organism, verbose=FALSE, chunkSize=10000000) {	
	return(1-sequenceCalc(locations, window, organism, DNAString("N"), verbose=verbose, chunkSize=chunkSize)/window)
}
