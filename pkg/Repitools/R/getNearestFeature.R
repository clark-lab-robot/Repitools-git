getNearestFeature <- function(peaksTable, featuresTable)
{
	require(ChIPpeakAnno)
	
	peaksTable$primaryKey = 1:nrow(peaksTable)
	featuresTable$primaryKey = 1:nrow(featuresTable)
	
	if("strand" %in% colnames(peaksTable))
	{
		peaksRangedData <- RangedData(ID = peaksTable$primaryKey, space = peaksTable$chr, strand = peaksTable$strand, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	} else {

		peaksRangedData <- RangedData(ID = peaksTable$primaryKey, space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	}

	if("strand" %in% colnames(featuresTable))
	{
		featuresRangedData <- RangedData(ID = featuresTable$primaryKey, space = featuresTable$chr, strand = featuresTable$strand, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "TSS"
	} else {
		featuresRangedData <- RangedData(ID = featuresTable$primaryKey, space = featuresTable$chr, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "middle"
	}	
	
	annoPeaks <- as.data.frame(annotatePeakInBatch(peaksRangedData, AnnotationData = featuresRangedData, PeakLocForDistance = "middle", FeatureLocForDistance = featureLoc))
	annoPeaks <- annoPeaks[order(as.numeric(as.character(annoPeaks$peak))), ]
	orginalOrder <- order(unlist(values(peaksRangedData)[, "ID"]))
	annoPeaks <- annoPeaks[orginalOrder, ]
	colnames(annoPeaks)[which(colnames(annoPeaks) == "insideFeature")] <- "peakLocation"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "distancetoFeature")] <- "distance"
	annoPeaks <- cbind(peaksTable[, 1:(ncol(peaksTable) - 1)], featuresTable[as.numeric(as.character(annoPeaks$feature)), 1:(ncol(featuresTable) - 1)], annoPeaks[, c("peakLocation", "distance")])
	return(annoPeaks)
}
