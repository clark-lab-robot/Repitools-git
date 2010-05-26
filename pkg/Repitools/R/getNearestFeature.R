getNearestFeature <- function(peaksTable, featuresTable)
{
	require(ChIPpeakAnno)
	if(packageDescription("ChIPpeakAnno")$Version < 1.5.3)
		warning("Your version of ChIPpeakAnno is not 1.5.3 or greater. Versions before this had bugs in them.")
	
	if("strand" %in% colnames(peaksTable))
	{
		peaksRangedData <- RangedData(space = peaksTable$chr, strand = peaksTable$strand, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	} else {
		peaksRangedData <- RangedData(space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	}

	if("strand" %in% colnames(featuresTable))
	{
		featuresRangedData <- RangedData(name = featuresTable$name, space = featuresTable$chr, strand = featuresTable$strand, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "TSS"
	} else {
		featuresRangedData <- RangedData(name = featuresTable$name, space = featuresTable$chr, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "middle"
	}	
	
	annoPeaks <- as.data.frame(annotatePeakInBatch(peaksRangedData, AnnotationData = featuresRangedData, PeakLocForDistance = "middle"))
	annoPeaks <- annoPeaks[order(annoPeaks$peak), c("feature", "space", "strand", "start_position", "end_position", "insideFeature", "distancetoFeature")]
	colnames(annoPeaks)[which(colnames(annoPeaks) == "insideFeature")] <- "peakLocation"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "distancetoFeature")] <- "distance"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "start_position")] <- "featureStart"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "end_position")] <- "featureEnd"
	annoPeaks$feature <- featuresTable[annoPeaks$feature, "name"]
	annoPeaks <- annoPeaks[, -2]
	return(annoPeaks)
}
