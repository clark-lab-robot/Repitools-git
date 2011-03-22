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


