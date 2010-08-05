setMethodS3("sequenceCalc", "GRanges", function(locations, organism, pattern, fixed=TRUE, Nmask=FALSE, positions=FALSE, verbose=FALSE, chunkSize=10000000, ...) {
	chr.length <- seqlengths(organism)
	elementMetadata(locations)$chunk <- paste(seqnames(locations), trunc(start(locations)/chunkSize), sep="/")
	scores <- if (positions) vector(mode='list', length=length(locations)) else numeric(length(locations))
	for(chunk in unique(elementMetadata(locations)$chunk)) {		
		if (verbose) cat(chunk, " ")
		chr <- gsub("/.*", "", chunk)
		thisChunk <- which(elementMetadata(locations)$chunk==chunk)
		#Grab the smallest chunk of the chromosome possible
		chrStart <- max(1, min(start(locations)[thisChunk]))
		chrEnd <- min(chr.length[chr], max(end(locations)[thisChunk]))
		if (verbose) cat(chrStart, chrEnd, "\n")
		chrseq <- subseq(organism[[chr]], chrStart, chrEnd)
		if (Nmask) chrseq <- mask(chrseq, pattern=pattern) 
		matches <- start(matchPattern(pattern, chrseq, fixed=fixed))
		if (length(matches)==0) next
		matches.IRanges <- IRanges(start=matches, width=1)
		loc.IRanges <- IRanges(start = start(locations)[thisChunk] - chrStart + 1, end = end(locations)[thisChunk] - chrStart + 1)
		loc.overlaps <- findOverlaps(query=matches.IRanges, subject=loc.IRanges)
		loc.overlaps <- tapply(loc.overlaps@matchMatrix[,1], loc.overlaps@matchMatrix[,2], list)
		if (length(loc.overlaps)==0) next
		thisChunk2 <- thisChunk[as.integer(names(loc.overlaps))]
		if (positions) {
			elementMetadata(locations)$position <- round((start(locations) + end(locations)) / 2)
			temp <- mapply(function(x,y) matches[x]-y+chrStart-1, loc.overlaps, elementMetadata(locations)$position[thisChunk2], SIMPLIFY=FALSE)
			scores[thisChunk2] <- temp
		} else scores[thisChunk2] <- sapply(loc.overlaps, length)
	}
	if (verbose) cat("\n")
	return(scores)
})

setMethodS3("sequenceCalc", "data.frame", function(locations, window=500, organism, pattern, fixed=TRUE, Nmask=FALSE, positions=FALSE, verbose=FALSE, chunkSize=10000000, ...) {
	chr.length <- seqlengths(organism)
	if(is.null(locations$position)) locations$position <- ifelse(locations$strand == '+', locations$start, locations$end)
	if (any((locations$position<window/2) | (locations$position+window/2>chr.length[locations$chr])))
		warning("Not all locations' windows are within chromosome boundaries.")
	locations$chunk <- paste(locations$chr, trunc(locations$position/chunkSize), sep="/")
	scores <- if (positions) vector(mode='list', length=nrow(locations)) else numeric(nrow(locations))
	for(chunk in unique(locations$chunk)) {		
		if (verbose) cat(chunk, " ")
		chr <- gsub("/.*", "", chunk)
		thisChunk <- which(locations$chunk==chunk)
		#Grab the smallest chunk of the chromosome possible
		chrStart <- max(1, min(locations$position[thisChunk])-window/2+1)
		chrEnd <- min(chr.length[chr], max(locations$position[thisChunk])+window/2)
		if (verbose) cat(chrStart, chrEnd, "\n")
		chrseq <- subseq(organism[[chr]], chrStart, chrEnd)
		if (Nmask) chrseq <- mask(chrseq, pattern=pattern) 
		matches <- start(matchPattern(pattern, chrseq, fixed=fixed))
		if (length(matches)==0) next
		matches.IRanges <- IRanges(start=matches, width=1)
		loc.IRanges <- IRanges(start=locations$position[thisChunk]-window/2+1-chrStart, end=locations$position[thisChunk]+window/2-chrStart)
		loc.overlaps <- findOverlaps(query=matches.IRanges, subject=loc.IRanges)
		loc.overlaps <- tapply(loc.overlaps@matchMatrix[,1], loc.overlaps@matchMatrix[,2], list)
		if (length(loc.overlaps)==0) next
		thisChunk2 <- thisChunk[as.integer(names(loc.overlaps))]
		if (positions) {
			temp <- mapply(function(x,y) matches[x]-y+chrStart-1, loc.overlaps, locations$position[thisChunk2], SIMPLIFY=FALSE)
			scores[thisChunk2] <- temp
		} else scores[thisChunk2] <- sapply(loc.overlaps, length)
	}
	if (verbose) cat("\n")
	return(scores)
})
