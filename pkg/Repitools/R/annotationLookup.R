annotationBlocksLookup <- function(probes, annotation, probeIndex=NULL, verbose=TRUE) {
	require(GenomicRanges)
	probeRanges <- GRanges(seqnames = probes$chr, ranges = IRanges(probes$position, width = 1))
	if("strand" %in% colnames(annotation))
		annoStrands <- annotation$strand
	else
		annoStrands <- Rle('*', nrow(annotation))
	if("name" %in% colnames(annotation))
		annoNames <- annotation$name
	else
		annoNames <- rownames(annotation)
	if(is.null(probeIndex)) probeIndex <- 1:nrow(probes)
	annotationRanges <- GRanges(seqnames = annotation$chr, ranges = IRanges(annotation$start, annotation$end), strand = annoStrands)
	matches <- findOverlaps(annotationRanges, probeRanges)@matchMatrix

	indexes <- vector("list", nrow(annotation))
	offsets <- vector("list", nrow(annotation))
	names(indexes) <- annoNames
	names(offsets) <- annoNames
	for(index in 1:nrow(annotation))
	{
		currAnno <- which(matches[, 1] == index)
		if(length(currAnno) != 0)
		{
			indexes[[index]] <- probeIndex[matches[currAnno, 2]]
			namedOffsets <- start(probeRanges)[indexes[[index]]] - start(annotationRanges)[index]
			names(namedOffsets) <- indexes[[index]]
			offsets[[index]] <- namedOffsets
		}
	}

	if (verbose) cat("Mapping done.\n")
	return(list(indexes = indexes, offsets = offsets))
}

annotationLookup <- function(probes, annotation, bpUp, bpDown, probeIndex=NULL, verbose=TRUE) {
	newStart <- numeric(nrow(annotation))
	newEnd <- numeric(nrow(annotation))
	newStrand <- if (is.null(annotation$strand)) '*' else annotation$strand
	whichPlus <- which(newStrand == '+')
	whichMinus <- which(newStrand == '-')
	whichBoth <- which(newStrand == '*')
	newStart[whichPlus] <- annotation$start[whichPlus] - bpUp
	newStart[whichMinus] <- annotation$end[whichMinus] - bpDown
	newStart[whichBoth] <- round((annotation$start[whichBoth] + annotation$end[whichBoth]) / 2) - bpUp
	newEnd[whichPlus] <- annotation$start[whichPlus] + bpDown
	newEnd[whichMinus] <- annotation$end[whichMinus] + bpUp
	newEnd[whichBoth] <- round((annotation$start[whichBoth] + annotation$end[whichBoth]) / 2) + bpDown
	
	regions <- data.frame(chr = annotation$chr, start = newStart, end = newEnd, strand = newStrand, name = if("name" %in% colnames(annotation)) annotation$name else 1:nrow(annotation))
	if(is.null(probeIndex)) probeIndex <- 1:nrow(probes)
	annot <- annotationBlocksLookup(probes, regions, probeIndex, verbose)

	annot$offsets <- mapply(function(aStrand, aPos, pInds){
		if(aStrand == '+')
			newOffsets <- (aPos - probes$position[match(pInds, probeIndex)]) * -1
		else
			newOffsets <- aPos - probes$position[match(pInds, probeIndex)]
		names(newOffsets) <- pInds
		return(newOffsets)
	}, annotation$strand, ifelse(annotation$strand == '+', annotation$start, annotation$end), annot$indexes, SIMPLIFY = FALSE)
	names(annot$offsets) <- names(annot$indexes)

	if (verbose) cat("Mapping done.\n")
	return(annot)
}

annotationBlocksCounts <- function(rs, annotation, seqLen=NULL, verbose=TRUE) {
    require(GenomicRanges)
    if (!class(rs) == "GRangesList")
    		stop("rs must be a GRangesList")
    if (class(annotation)=="data.frame") {
    	if (is.null(annotation$name)) annotation$name <- 1:nrow(annotation)
        annotation <- GRanges(seqnames=annotation$chr, ranges = IRanges(start = annotation$start, end = annotation$end), name = annotation$name)
    } else {
    	stopifnot(class(annotation)=="GRanges")
    	if(!"name" %in% names(elementMetadata(annotation))) elementMetadata(annotation)$name <- 1:nrow(annotation)
    }

    anno.counts <- matrix(as.integer(NA), nrow=length(annotation), ncol=length(rs), dimnames=list(elementMetadata(annotation)[, "name"], names(rs)))
    rs <- endoapply(rs, resize, seqLen)
    for (i in 1:length(rs)) {
        if(verbose) cat("Counting in", names(rs)[i], "\n")
	anno.counts[, i] <- countOverlaps(annotation, rs[i])
    }
    if(verbose)
    	cat("Counting successful.\n")
    return(anno.counts)
}

annotationCounts <- function(rs, annotation, bpUp, bpDown, seqLen=NULL, verbose=TRUE) {
	require(GenomicRanges)
	if (!class(rs) == "GRangesList")
		stop("rs must be a GRangesList")
	anno = annotation
	if (is.null(anno$strand)) anno$strand <- "*"
	anno$position <- if(anno$strand == '+') anno$start else if(anno$strand == '-') anno$end else round((anno$start + anno$end) / 2)
	if (is.null(anno$name)) anno$name <- 1:nrow(annotation)
	anno$start=ifelse(anno$strand=="+", anno$position-bpUp, anno$position-bpDown)
        anno$end=ifelse(anno$strand=="+", anno$position+bpDown, anno$position+bpUp)
	annotationBlocksCounts(rs, anno, seqLen, verbose)
}
