getNearestFeature <- function(peaksTable, featuresTable)
{
	require(ChIPpeakAnno)
	
	if("strand" %in% colnames(peaksTable))
	{
		peaksRangedData <- RangedData(space = peaksTable$chr, strand = peaksTable$strand, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
	} else {
		if("name" %in% colnames(peaksTable)) # Maybe the peaks are CpG islands from a database, and do have a name.
		{
			peaksRangedData <- RangedData(name = peaksTable$name, space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
		} else {
			peaksRangedData <- RangedData(space = peaksTable$chr, ranges = IRanges(start = peaksTable$start, end = peaksTable$end))
		}
	}

	RDspaces <- space(peaksRangedData)
	RDstarts <- start(peaksRangedData)
	RDends <- end(peaksRangedData)
	peaksRDtoTableMap <- mapply(function(chr, start, end) which(RDspaces == chr & RDstarts == start & RDends == end), 
                                   peaksTable$chr, peaksTable$start, peaksTable$end)

	if("strand" %in% colnames(featuresTable))
	{
		featuresRangedData <- RangedData(name = featuresTable$name, space = featuresTable$chr, strand = featuresTable$strand, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "TSS"
	} else {
		featuresRangedData <- RangedData(name = featuresTable$name, space = featuresTable$chr, ranges = IRanges(start = featuresTable$start, end = featuresTable$end))
		featureLoc <- "middle"
	}	
	
	annoPeaks <- as.data.frame(annotatePeakInBatch(peaksRangedData, AnnotationData = featuresRangedData, PeakLocForDistance = "middle"))
	annoPeaks$feature <- unlist(values(featuresRangedData)[, "name"])[as.numeric(annoPeaks$feature)]
	annoPeaks <- annoPeaks[order(as.numeric(annoPeaks$peak)), ]
	annoPeaks <- annoPeaks[peaksRDtoTableMap, ]
	colnames(annoPeaks)[which(colnames(annoPeaks) == "insideFeature")] <- "peakLocation"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "distancetoFeature")] <- "distance"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "start_position")] <- "featureStart"
	colnames(annoPeaks)[which(colnames(annoPeaks) == "end_position")] <- "featureEnd"
	annoPeaks <- annoPeaks[, -2]
	return(annoPeaks)
}
