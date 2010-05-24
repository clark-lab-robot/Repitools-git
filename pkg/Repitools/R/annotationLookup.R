annotationBlocksLookup <- function(probes, annotation, probeIndex=NULL, verbose=TRUE) {
#probes = dataframe of $chr, $position and $strand ("+" or "-")
#annotation = dataframe of $chr, $start, $end, $strand ("+" or "-") and $name or rownames = annotation name


	processChromosome <- function(probePositions, annotation) {
		numAnnot = nrow(annotation)
		annotProbes = list(indexes=vector(mode='list', length=numAnnot), offsets=vector(mode='list', length=numAnnot))
		if (length(probePositions)==0) return(annotProbes) #if no probes on this chromosome return empty annotations	
		require(IRanges)
		probes.IRanges <- IRanges(start=probePositions, width=1)
		annotation.IRanges <- IRanges(start=annotation$start, end=annotation$end)
		anno.overlaps <- findOverlaps(query=probes.IRanges, subject=annotation.IRanges)
		anno.overlaps <- tapply(anno.overlaps@matchMatrix[,1], anno.overlaps@matchMatrix[,2], list)
		annotProbes$indexes[as.integer(names(anno.overlaps))] <- anno.overlaps
		annotProbes$offsets <- mapply(function(x ,y) probePositions[x]-y, annotProbes$indexes, annotation$start, SIMPLIFY = FALSE)
		annotProbes$indexes <- lapply(annotProbes$indexes, function(x) as.integer(names(probePositions[x])))
		return(annotProbes)
  }


	if (is.null(annotation$strand)) { #dont bother with strandedness anywhere
		probesStrandChr <- probes$chr
		annotationStrandChr <- annotation$chr
	} else { #create separate plus & minus chromosomes
		if (is.null(probes$strand)) stop("error: if 'annotation' contains strand information, 'probes' must contain strand information as well")
		probesStrandChr <- paste(probes$chr, probes$strand, sep="")
		annotationStrandChr <- paste(annotation$chr, annotation$strand, sep="")
	}

	#split by strand AND chromosome simultaneously
	annotChr = split(1:nrow(annotation), annotationStrandChr)
	annot = list(indexes=vector(mode='list', length=nrow(annotation)), offsets=vector(mode='list', length=nrow(annotation)))
	if (verbose) cat("Processing mapping between probes and genes.\n")

	for (i in annotChr) {
		thisChr = annotationStrandChr[i[1]]
		
		#Grab the subset of probes on that chromosome
		tempIndex = which(probesStrandChr==thisChr)
		tempProbes = probes$position[tempIndex]
		#Use probeIndex supplied, or assume probes are in order
		if (is.null(probeIndex)) names(tempProbes) <- tempIndex else names(tempProbes) <- probeIndex[tempIndex]

		#Process the chromosome
		tempAnnot = processChromosome(tempProbes, annotation[i,])
		annot$indexes[i] = tempAnnot$indexes
		annot$offsets[i] = tempAnnot$offsets
	}
	if (!is.null(rownames(annotation))) {
		names(annot$indexes) <- annotation$name
		names(annot$offsets) <- annotation$name
	} else {
		names(annot$indexes) <- rownames(annotation)
		names(annot$offsets) <- rownames(annotation)
	}
	if (verbose) cat("Mapping done.\n")
	return(annot)
	#returns $indexes = a list for each annotation entry with the indexes of the probes within the block
	#	 $offsets = a list for each annotation entry with the offsets from the beginning of the block
	
}



annotationLookup <- function(probes, annotation, bpUp, bpDown, probeIndex=NULL, verbose=TRUE) {
#probes = dataframe of $chr and $position
#annotation = dataframe of $chr, $position, $strand ("+" or "-") and $name or rownames = annotation name
#if annotation has no strand, assume are + strand
	if (is.null(annotation$strand)) annotation$strand <- "+"
	annotationTemp <- data.frame(chr=annotation$chr, 
                                     start=annotation$position+ifelse(annotation$strand=="+",-bpUp, -bpDown),
                                     end=annotation$position+ifelse(annotation$strand=="+",+bpDown, +bpUp),
                                     name=rownames(annotation), stringsAsFactors=F)
	annot <- annotationBlocksLookup(probes, annotationTemp, probeIndex, verbose)

	annot$offsets[annotation$strand=="+"] = lapply(annot$offsets[annotation$strand=="+"], function(x,bpOff) {return(x-bpOff)}, bpUp)
	annot$offsets[annotation$strand=="-"] = lapply(annot$offsets[annotation$strand=="-"], function(x,bpOff) {return(rev(bpOff-x))}, bpDown)
	annot$indexes[annotation$strand=="-"] = lapply(annot$indexes[annotation$strand=="-"], rev)
	if (!is.null(rownames(annotation))) {
		names(annot$indexes) <- annotation$name
		names(annot$offsets) <- annotation$name
	} else {
		names(annot$indexes) <- rownames(annotation)
		names(annot$offsets) <- rownames(annotation)
	}

	return(annot)
}


annotationBlocksCounts <- function(rs, annotation, seqLen=NULL, verbose=TRUE) {
    if (class(rs)=="GenomeData") rs <- GenomeDataList(list(rs))
    if (class(annotation)=="data.frame") {
        annotation <- RangedData(Ranges=IRanges(start=annotation$start, end=annotation$end), space=annotation$chr)
    } else stopifnot(class(annotation)=="RangedData") 
    anno.counts <- matrix(as.integer(NA), nrow=nrow(annotation), ncol=length(rs), dimnames=list(annotation$name, names(rs)))
    if (!class(rs[[1]][[1]])=="IRanges") {
        if (is.null(seqLen)) stop("If rs has not been extended, seqLen must be supplied")
        rs <- chipseq::extendReads(rs, seqLen=seqLen)
    }
    #convert rs to RangesList
#    rs <- lapply(rs, function(x) RangesList(as.list(x)))
    anno.chr <- names(annotation)
    for (i in 1:length(rs)) {
        
        if(verbose)
            cat(names(rs)[i], "\n", sep="")
        #Handle missing chromosomes
        rs.chr <- names(rs[[i]])
        nochr <- which(!anno.chr %in% rs.chr)
        nochr.list <- rep(list(IRanges()), length(nochr))
        names(nochr.list) <- anno.chr[nochr]
        rs[[i]] <- RangesList(c(as.list(rs[[i]]), nochr.list))
        rs.chr <- names(rs[[i]])
        both.chr <- intersect(anno.chr, rs.chr)
        anno.counts[,i] <- unlist(countOverlaps(annotation[both.chr], rs[[i]][both.chr], type="any"))
    }
    anno.counts
}

annotationCounts <- function(rs, annotation, bpUp, bpDown, seqLen=NULL, verbose=TRUE) {
	require(chipseq)
	if (class(rs)=="GenomeData") rs <- GenomeDataList(list(rs))
	if (is.null(annotation$strand)) annotation$strand <- "+"
	if (is.null(annotation$name)) annotation$name <- 1:nrow(annotation)
	anno <- data.frame(chr=annotation$chr,
                     start=ifelse(annotation$strand=="+", annotation$position-bpUp, annotation$position-bpDown), 
                       end=ifelse(annotation$strand=="+", annotation$position+bpDown, annotation$position+bpUp),
                      name=annotation$name)
	annotationBlocksCounts(rs, anno, seqLen, verbose)
}

