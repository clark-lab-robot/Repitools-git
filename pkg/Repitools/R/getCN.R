getCN <- function(regionsInfo, readsInputs)
{
	require(DNAcopy)
	require(IRanges)	

	Mvalues <- log2((readsInputs[, 2]/ sum(readsInputs[, 2])) / (readsInputs[, 1] / sum(readsInputs[, 1])))
	cnObject <- CNA(chrom = regionsInfo$inputs[, "chr"], maploc = regionsInfo$inputs[, "start"], data.type = "logratio", genomdat = Mvalues, sampleid = paste(colnames(readsInputs)[2], "/", colnames(readsInputs)[1], "Fold Change"))
	readsCounts <- colSums(readsInputs)
	wts <- ((readsCounts[2] - readsInputs[, 2]) / (readsInputs[, 2] * readsCounts[2]) + (readsCounts[1] - readsInputs[, 1]) / (readsInputs[, 1] * readsCounts[1]))^-1
	nonZero <- which(wts > 0)
	cnObject <- cnObject[nonZero, ]
	wts <- wts[nonZero]
	cnObject <- segment(smooth.CNA(cnObject), weights = wts, p.method = "perm")
	cnObject$out[, "loc.end"] <- cnObject$out[, "loc.end"] + regionsInfo$inputs[1, "end"] - regionsInfo$inputs[1, "start"] # Extend CNV region to the end of the interval, since all positions are starts.
	
	enrichRegionsRanges <- RangedData(IRanges(start = regionsInfo$enriched[, "start"], end = regionsInfo$enriched[, "end"]), space = regionsInfo$enriched[, "chr"])
	copyRanges <- RangedData(IRanges(start = cnObject$out[, "loc.start"], end = cnObject$out[, "loc.end"]), space = cnObject$out[, "chrom"])
	map <- findOverlaps(enrichRegionsRanges, copyRanges, select = "first")
							   
	scaleFactors <- numeric()
	cnPerChr <- split(cnObject$out, cnObject$out[, "chrom"])
	for(chrIndex in 1:length(map))	
	{	
		whichCurrChr <- which(regionsInfo$enriched[, "chr"] == names(map)[[chrIndex]])
		if(!is.null(cnPerChr[[names(map)[chrIndex]]]))
		{
			scaleFactors[whichCurrChr] <- 2^(cnPerChr[[names(map)[chrIndex]]][map[[names(map)[chrIndex]]], "seg.mean"])
		} else {
			scaleFactors[whichCurrChr] <- NA
		}
	}
	
	scaleFactors[which(is.na(scaleFactors))] <- 1 # Regions of no change.
	

	factorMatrix <- cbind(1, scaleFactors)
	
	colnames(factorMatrix) <- colnames(readsInputs)
	rownames(factorMatrix) <- rownames(regionsInfo$enriched)
	return(factorMatrix)
}

