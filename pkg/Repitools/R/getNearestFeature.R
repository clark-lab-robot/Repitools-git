getNearestFeature <- function(peaksTable, featuresTable)
{
	require(ChIPpeakAnno)
	
	peaksTable$primaryKey = 1:nrow(peaksTable)
	featuresTable$primaryKey = 1:nrow(featuresTable)
	
	if("strand" %in% colnames(peaksTable))
	{
		peaksRangedData <- RangedData(ID = peaksTable$primaryKey, space = peaksTable$chr, strand = peaksTable$strand, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	} else {
		if("name" %in% colnames(peaksTable)) # Maybe the peaks are CpG islands from a database, and do have a name.
		{
			peaksRangedData <- RangedData(ID = peaksTable$primaryKey, name = peaksTable$name, space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
		} else {
			peaksRangedData <- RangedData(ID = peaksTable$primaryKey, space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
		}
	}

	if("strand" %in% colnames(featuresTable))
	{
		featuresRangedData <- RangedData(ID = featuresTable$primaryKey, name = featuresTable$name, space = featuresTable$chr, strand = featuresTable$strand, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "TSS"
	} else {
		featuresRangedData <- RangedData(ID = featuresTable$primaryKey, name = featuresTable$name, space = featuresTable$chr, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "middle"
	}	
	
	annoPeaks <- as.data.frame(annotatePeakInBatch(peaksRangedData, AnnotationData = featuresRangedData))
	featLoc <- ifelse(as.data.frame(values(featuresRangedData))[, "strand"] == '+', start(featuresRangedData), end(featuresRangedData))
	annoPeaks[, "distancetoFeature"] <- abs(round((annoPeaks[, 2] + annoPeaks[, 3]) / 2) - featLoc[as.numeric(as.character(annoPeaks$feature))])
	annoPeaks$feature <- as.data.frame(unlist(values(featuresRangedData)[, "name"]))[, 1][as.numeric(as.character(annoPeaks$feature))]
	annoPeaks <- annoPeaks[order(as.numeric(as.character(annoPeaks$peak))), ]
	orginalOrder <- order(unlist(values(peaksRangedData)[, "ID"]))
	annoPeaks <- annoPeaks[orginalOrder, ]
	colnames(annoPeaks)[which(colnames(annoPeaks) == "insideFeature")] <- "peakLocation"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "distancetoFeature")] <- "distance"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "start_position")] <- "featureStart"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "end_position")] <- "featureEnd"
	annoPeaks <- annoPeaks[, c(1:3, 8:10, 7, 11, 12)]
	if("name" %in% colnames(peaksTable))
			annoPeaks <- cbind(name = peaksTable$name, annoPeaks)
	return(annoPeaks)
}
