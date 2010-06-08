setMethodS3("writeWig", "GenomeDataList", function(rs, seqLen, design=NULL, sample=20, dropZero=TRUE, normalise=TRUE, verbose=TRUE, ...) {
	scipen <- getOption("scipen")
	options(scipen=100)
	if (verbose) cat("Extending the reads by ", seqLen,"bp\n", sep="" )
	rs.ext <- extendReads(rs, seqLen)
	if (verbose) cat("Creating coverage objects\n")
	rs.cov <- gdapply(rs.ext, coverage)
	if (normalise) {
		if (verbose) cat("Calculating normalisation counts\n")
		norm <- laneCounts(rs)/1000000
	}
	if (is.null(design)) {
		if (verbose) cat("Default design matrix\n")
		design <- matrix(0, nrow=length(rs), ncol=length(rs), dimnames=list(NULL, paste(names(rs), ".wig.gz", sep="")))
		for (i in 1:length(rs)) design[i,i] <- 1
	}
	for (i in 1:ncol(design)) {
		filename <- colnames(design)[i]
		if (verbose) cat("Processing column",i,"of design matrix -",filename,"\n")
		f1 <- if (grepl("gz$",filename)) gzfile(filename, open="wt") else file(filename, open="wt")

		thisDes <- which(design[,i]!=0)
		thisMul <- design[thisDes,i]
		rs.this <- rs.cov[thisDes]		
		for (j in 1:length(rs.this[[1]])) {
			if (verbose) cat(names(rs.this[[1]])[j]," ")
			#Write header for this chromosome
			writeLines(paste("variableStep chrom=", names(rs.this[[1]])[j]," span=", sample, sep=""), f1)
			#ensure all lanes being used have the same length
			chrMax <- max(sapply(rs.this, function(x) length(x[[j]])))
			for (k in 1:length(rs.this)) rs.this[[k]][[j]] <- c(rs.this[[k]][[j]], Rle(0, chrMax-length(rs.this[[k]][[j]]))) 
			bp <- 1:(length(rs.this[[1]][[j]]) %/% sample)*sample
			bp.score <- as.integer(rs.this[[1]][[j]][bp])/norm[thisDes[1]]*thisMul[1]
			if (length(thisDes)>1) for (k in 2:length(thisDes)) bp.score <- bp.score + as.integer(rs.this[[k]][[j]][bp])/norm[thisDes[k]]*thisMul[k]
			temp <- if (dropZero) paste(bp[bp.score!=0], bp.score[bp.score!=0], sep="\t")
				else paste(bp, bp.score, sep="\t")
			writeLines(temp, f1)
		}
		if (verbose) cat("\n")
		close(f1)
	}
	options(scipen=scipen)
})

setMethodS3("writeWig", "AffymetrixCelSet", function(rs, design=NULL, log2adjust=TRUE, verbose=TRUE, ...) {
	scipen <- getOption("scipen")
	options(scipen=100)
	if (is.null(design)) {
		if (verbose) cat("Default design matrix\n")
		design <- matrix(0, nrow=length(rs), ncol=length(rs), dimnames=list(NULL, paste(getNames(rs), ".wig.gz", sep="")))
		for (i in 1:length(rs)) design[i,i] <- 1
	}

	w <- which( rowSums(design != 0) > 0 )
	rs <- extract(rs, w, verbose=verbose)
	probePositions <- getProbePositionsDf( getCdf(rs), verbose=verbose )
	dmP <- extractMatrix(rs, cells=probePositions$index, verbose=verbose)
	if(log2adjust == TRUE) diffs <- log2(dmP) %*% design[w,] else diffs <- dmP %*% design[w,]
	
	for (i in 1:ncol(design)) {
		filename <- colnames(design)[i]
		if (verbose) cat("Processing column",i,"of design matrix -",filename,"\n")
		f1 <- if (grepl("gz$",filename)) gzfile(filename, open="wt") else file(filename, open="wt")
		for (chr in unique(probePositions$chr)) {
			if (verbose) cat(chr," ")
			thisChr <- which(probePositions$chr==chr)
			#Write header for this chromosome
			writeLines(paste("variableStep chrom=", chr," span=1",sep=""), f1)
			temp <- paste(probePositions$position[thisChr], diffs[thisChr,i], sep="\t")
			writeLines(temp, f1)
		}
				if (verbose) cat("\n")
		close(f1)
	}
	options(scipen=scipen)
})
