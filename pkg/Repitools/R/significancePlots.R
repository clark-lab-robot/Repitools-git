setMethodS3("significancePlots", "default", function(rs, coordinatesTable, design=NULL, upStream=7500, downStream=2500, by=100, bw=300, total.lib.size=TRUE, seqLen=NULL, verbose=FALSE, ...) {
	coordinatesTable$position <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
	rownames(coordinatesTable) <- coordinatesTable$name
	blockPos <- seq.int(-upStream, downStream, by)
	if (verbose) cat("made blockPos\n")
	annoBlocks <- data.frame(chr=rep(coordinatesTable$chr, each=length(blockPos)),
		start=rep(coordinatesTable$position-bw, each=length(blockPos)),
		end=rep(coordinatesTable$position+bw, each=length(blockPos)),
		strand=rep(coordinatesTable$strand, each=length(blockPos)))
	annoBlocks$start[annoBlocks$strand=="+"] <- annoBlocks$start[annoBlocks$strand=="+"] + blockPos
	annoBlocks$end[annoBlocks$strand=="+"] <- annoBlocks$end[annoBlocks$strand=="+"] + blockPos
	annoBlocks$start[annoBlocks$strand=="-"] <- annoBlocks$start[annoBlocks$strand=="-"] - blockPos
	annoBlocks$end[annoBlocks$strand=="-"] <- annoBlocks$end[annoBlocks$strand=="-"] - blockPos
	if (verbose) cat("made annoBlocks\n")
	if (!is.null(design)) {
		stopifnot(all(design %in% c(-1,0,1)), nrow(design)==length(rs))
		inUse <- !apply(design==0,1,all)
		design <- design[inUse,, drop=FALSE]
	} else inUse <- rep(TRUE, length(rs))
	annoCounts <- annotationBlocksCounts(rs[inUse], annoBlocks, seqLen, verbose)
	if (total.lib.size) {
		if (verbose) cat("normalising to total library sizes\n")
		totalReads <- if(class(rs) == "GenomeDataList") laneCounts(rs[inUse]) else elementLengths(rs[inUse])
		annoCounts <- t(t(annoCounts)/totalReads)*1000000
	}
	if (verbose) cat("made annoCounts\n")
	if (!is.null(design)) {
		if (verbose) cat("applying design matrix\n")
		design <- apply(design, 2, function(x) {
					x[x==1] <- 1/sum(x==1)
					x[x==-1] <- -1/sum(x==-1)
					return(x)
				})
		annoCounts <- annoCounts %*% design 
	}
	annoTable <- matrix(1:nrow(annoCounts), byrow=TRUE, ncol=length(blockPos), nrow=nrow(coordinatesTable), dimnames=list(NULL, blockPos))
	if (verbose) cat("made annoTable\n")
	significancePlots(annoCounts, annoTable, removeZeros=FALSE, useMean=TRUE, ...)
})

setMethodS3("significancePlots", "AffymetrixCelSet", function(cs, probeMap=NULL, coordinatesTable=NULL, upStream=7500, downStream=2500, by=100, bw=300, log2adjust=TRUE, verbose=FALSE, ...) {			
	if (is.null(probeMap)) {
		if (is.null(coordinatesTable)) stop("Either probeMap or coordinatesTable must be supplied!")
		probePositions <- getProbePositionsDf( getCdf(cs), verbose=verbose )
		coordinatesTable$position <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
		rownames(coordinatesTable) <- coordinatesTable$name
		
		# run lookup twice.  first to get a list of smaller list of probes to use
		annot <- annotationLookup(probePositions, coordinatesTable, upStream+bw, downStream+bw, verbose=verbose)
		pb <- unique(unlist(annot$indexes, use.names=FALSE))
		probePositions <- probePositions[pb,]
		annot <- annotationLookup(probePositions, coordinatesTable, upStream+bw, downStream+bw, verbose=verbose)
		lookupT <- makeWindowLookupTable(annot$indexes, annot$offsets,
				starts = seq(-upStream-bw, downStream-bw, by), ends = seq(-upStream+bw, downStream+bw, by))
	} else {
		if (verbose) cat("Using supplied probeMap\n")
		probePositions <- probeMap$probePositions
		lookupT <- probeMap$lookupT
	}
	
	dmM <- extractMatrix(cs, cells = probePositions$index, verbose = verbose)
	if (log2adjust) dmM <- log2(dmM)
	
	significancePlots(dmM, lookupT, ...)
	invisible(list(lookupT=lookupT, probePositions=probePositions))
})


setMethodS3("significancePlots", "matrix", function(dataMatrix, lookupTable, geneList, titles=colnames(dataMatrix), nSamples=1000, confidence=0.975, legend.plot="topleft", cols=rainbow(length(geneList)), removeZeros=TRUE, useMean=FALSE, ...) {
	#Test geneList for sanity
	for (i in 1:length(geneList)) if (class(geneList[[i]])=="logical") {
		if (length(geneList[[i]])!=nrow(lookupTable)) 
		  stop("boolean geneList element length must equal num of rows in lookupTable")
	} else if (class(geneList[[i]])=="integer") {
		if(max(geneList[[i]])>nrow(lookupTable)) 
		  stop("geneList element value greater than num of rows in lookupTable") 
	} else stop("geneList elements must a be boolean or integer vector")
	stopifnot(confidence>0.5, confidence<1)
	if (length(legend.plot)!=ncol(dataMatrix)) if (length(legend.plot)!=1) stop("legend.plot must be either same length as columns in dataMatrix or 1") else legend.plot <- rep(legend.plot, ncol(dataMatrix))
	x.p <- as.numeric(colnames(lookupTable))
	#grab Intensities for all genes
	geneList.max <- which.max(sapply(geneList, FUN=function(u) if(class(u)=="logical") sum(u) else length(u)))
	for (i in 1:ncol(dataMatrix)) {
		sMat <- .scoreIntensity(lookupTable, dataMatrix[,i], removeZeros=removeZeros, returnMatrix=TRUE, useMean=useMean)
		sMat.geneList <- lapply(geneList, FUN=function(u) sMat[u, ])
		if (useMean) trace.geneList <- sapply(sMat.geneList, FUN=function(u) apply(u, 2, mean, na.rm=TRUE)) else
		trace.geneList <- sapply(sMat.geneList, FUN=function(u) apply(u, 2, median, na.rm=TRUE))

		#choose nSamples random genelists
		inds <- lapply(seq_len(nSamples), FUN=function(u) sample(nrow(sMat), nrow(sMat.geneList[[geneList.max]])))
		if (useMean) meds <- sapply(inds, FUN=function(u) apply(sMat[u,], 2, mean, na.rm=TRUE)) else
		meds <- sapply(inds, FUN=function(u) apply(sMat[u,], 2, median, na.rm=TRUE))
		meds.conf <- apply(meds, 1, quantile, p=c(1-confidence, 0.5, confidence))
	
		#plot tiem
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="n", lty=c(2,1,2,1), lwd=c(1,3,1,3), xlab="Position relative to TSS", ylab="Signal", main=titles[i], ...)
		polygon(x=c(x.p, rev(x.p)), y=c(meds.conf[1,], rev(meds.conf[3,])), col="lightblue")
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="l", lty=c(2,1,2,rep(1, length(geneList))), lwd=c(1,3,1,rep(3, length(geneList))), add=TRUE, col=c("blue","blue","blue",cols))
		if (!is.na(legend.plot[i])) legend(legend.plot[i], legend=names(geneList), col=cols, lwd=3)
	}
})

