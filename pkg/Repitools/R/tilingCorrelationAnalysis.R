.makeWindowLookupList <- function(indexes, offsets, starts, ends) {

	getProbesInWindow <- function(indexes, offsets, start, end) {
	#indexes & offsets are integer vectors at this stage
	#start & end are atomic integers - offset co-ordinates we are interested in
		if (length(indexes)==0) return(integer()) else return(indexes[which(offsets>start&offsets<end)])
	}

	lookupList <- vector(mode='list', length=length(starts))
	names(lookupList) <- (starts+ends)/2
	for (i in 1:length(starts)) {
		cat("Processing",starts[i],"to",ends[i],"\n")
		lookupList[[i]] = vector(mode='list', length=length(indexes))
		if (!is.null(names(indexes))) names(lookupList[[i]]) <- names(indexes)
		for (j in 1:length(indexes)) lookupList[[i]][[j]] = getProbesInWindow(indexes[[j]], offsets[[j]], starts[i], ends[i])

	}
	return(lookupList)
}


makeWindowLookupTable <- function(indexes, offsets, starts, ends) {
	lookupTable <- matrix(NA,nrow=length(indexes),ncol=length(starts))
	colnames(lookupTable) <- (starts+ends)/2
	rownames(lookupTable) <- names(indexes)

	ind <- unlist(indexes)
	off <- unlist(offsets, use.names=FALSE)
	genes <- unlist(mapply(rep, 1:length(indexes), sapply(indexes, length)))

	off.IRanges <- IRanges(start=off, width=1)
	o <- findOverlaps(query=off.IRanges, subject=IRanges(start=starts, end=ends))@matchMatrix
	o <- tapply(o[,1], o[,2], list)

	for (i in 1:length(o)) {
		pos <- as.integer(names(o)[i])
		genes.mid <- tapply(o[[i]], genes[o[[i]]], function(x, midpt) {
			return(ind[x[which.min(abs(off[x]-midpt))]])
		}, (starts[pos]+ends[pos])/2)
		lookupTable[as.integer(names(genes.mid)),pos] <- genes.mid
	}
	return(lookupTable)
}




.scoreCorrelation <- function(lookup, intensities, correlateTo, minProbes=1, cor.method="pearson") {
	windowMeans <- function(indexes, intensities) {
		return(mean(intensities[indexes]))
	}

	useLookup = sapply(lookup, length)>=minProbes
	return(cor(sapply(lookup[useLookup], windowMeans, intensities),correlateTo[useLookup],method=cor.method))
}


.scoreIntensity <- function(lookup, intensities, minProbes=1, removeZeros=FALSE, returnMatrix=FALSE, useMean=FALSE) {
	windowMeans <- function(indexes, intensities) {
	    x <- intensities[indexes]
		if(removeZeros)
		  x <- x[x != 0]
		return(mean(x))
	}
	
	if( class(lookup)=="list") {
	  useLookup = sapply(lookup, length)>=minProbes
	  temp <- sapply(lookup[useLookup], windowMeans, intensities)
	  if (useMean) return(mean(temp, na.rm=TRUE)) else return(median(temp, na.rm=TRUE))
	} else if (class(lookup)=="matrix") {
	  d <- lookup
	  for(i in 1:ncol(lookup)) {
	    d[,i] <- intensities[  lookup[,i] ]
	  }
	  if (removeZeros)
	    d[d==0] <- NA
	  if (returnMatrix)
	    return(d)
	  if (useMean) return( apply(d,2,mean,na.rm=TRUE) ) else return( apply(d,2,median,na.rm=TRUE) )
	} else {
	  stop("lookup is neither a 'list' or a 'matrix'.")
	}
}


#Draw graph of correlation vs distance from TSS
correlationGraphs <- function(intensities, lookup, compareExpression, ...) {
	corScores = matrix(NA, nrow=ncol(intensities), ncol=length(lookup), dimnames=list(colnames(intensities),names(lookup)))
	for (i in 1:ncol(intensities)) corScores[i,] = sapply(lookup, .scoreCorrelation, intensities[,i], compareExpression, minProbes=2)

	cols = rainbow(nrow(corScores))
	plot(0, type="n", xlim=c(-7500,2500), ylim=c(-1,1), xlab="Distance from TSS", ylab="Correlation with Expression", ...)
	for (i in 1:nrow(corScores)) lines(as.integer(colnames(corScores)), corScores[i,], col=cols[i])
	legend("topleft", legend=rownames(corScores), col=cols, lty=1)
}

