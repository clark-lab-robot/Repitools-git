setGeneric("clusteredHeatmap", signature = "readsIP", function(readsIP, ...){standardGeneric("clusteredHeatmap")})

setMethod("clusteredHeatmap", "GRanges", function(readsIP, geneAnno, upPct = 50, downPct = 100, cvgResolution = 10, smoothingWidth = 2001, capQuantile = 0.95, nClusters = 5, expr, cols = colorpanel(64, "blue", "white", "red"), tName = "Clustered Enrichment")
{
	coverageIP <- coverage(readsIP)

	pctPositions <- seq(-upPct, downPct, cvgResolution)
	widths <- geneAnno$end - geneAnno$start + 1
	geneAnno$tss <- ifelse(geneAnno$strand == '+', geneAnno$start, geneAnno$end)
	
	TSScoverageList <- list()
	for(geneIndex in 1:nrow(geneAnno))
	{
		# Get smoothing region. Extra flanks so that we can smooth the edge of the region nicely.

		flank = floor(smoothingWidth / 2)

		if(geneAnno[geneIndex, "strand"] == '+')
		{
			region <- round(geneAnno[geneIndex, "start"] - upPct / 100 * widths[geneIndex] - flank) : round(geneAnno[geneIndex, "start"] + (downPct / 100 * widths[geneIndex]) + flank)
			coverageGene <- tryCatch(coverageIP[[geneAnno[geneIndex, "chr"]]][region], error = function(errMsg) "Error getting all bases' coverage.")
		} else {
			region <- round(geneAnno[geneIndex, "end"] - downPct / 100 * widths[geneIndex] - flank) : round(geneAnno[geneIndex, "end"] + (upPct / 100 * widths[geneIndex]) + flank)
			coverageGene <- tryCatch(rev(coverageIP[[geneAnno[geneIndex, "chr"]]][region]), error = function(errMsg) "Error getting all bases' coverage.")
		}

		# If could get all bases, then smooth.
		if(class(coverageGene) == "Rle") {coverageGene <- runmean(coverageGene, smoothingWidth, endrule = "constant")} else {TSScoverageList[[geneIndex]] <- "Error"; next;}

		# Get positions of precentage distances, then get the coverage at that point.
		samplingPositions <- round(flank + (0:(length(pctPositions)-1) / (length(pctPositions)-1)) * (length(region) - 2 * flank)) + 1
		reprCoverage <- as.numeric(coverageGene[samplingPositions])
		names(reprCoverage) <- pctPositions
		TSScoverageList[[geneIndex]] <- reprCoverage
	}

	# Find which gene didn't have all coverage available. Eliminate from expression vector.
	whichEliminate <- which(sapply(TSScoverageList, function(geneCoverage) length(geneCoverage) == 1))
	if(length(whichEliminate) > 0)
	{
		TSScoverageList <- TSScoverageList[-whichEliminate]
		expr <- expr[-whichEliminate]
	}
	TSScoverageTable <- do.call(rbind, TSScoverageList)
	maxIntensity = quantile(TSScoverageTable, capQuantile)
	TSScoverageTable[TSScoverageTable > maxIntensity] = maxIntensity

	# Do the k-means clustering.
	clusterID <- kmeans(TSScoverageTable, nClusters, iter.max = 100)$cluster

	# Re-arrange the ChIP and expression data.
	clusterOrder <- order(clusterID)
	TSScoverageTable <- TSScoverageTable[clusterOrder, ]
	expr <- expr[clusterOrder]

	# Plot heatmap.
	layout(rbind(c(1,2,3)), widths=c(1,3,1))
	par(mai=c(1.02,0.50,0.82,0.05))
	
	par(oma = c(0, 0, 0, 0))
	heatmapRange <- range(TSScoverageTable)
	nBins = length(cols)
	image(y=seq(1/nBins/2, 1-(1/nBins/2), 1/nBins), z=rbind(1:nBins), col = cols, axes = FALSE, xlab = "Signal Intensity", ylab = NA)
	axis(2, at=(0:nBins)/nBins, labels=format(seq(heatmapRange[1], heatmapRange[2], length.out = nBins + 1), digits = 2))

	par(mai=c(1.02,0.05,0.82,0.05))
	image(pctPositions, 1:nrow(TSScoverageTable), t(TSScoverageTable), xlab="Position (%) relative to TSS", yaxt = "n", ylab = "Gene", col = cols)

	# Add lines delimiting the cluster boundaries.
	boundaries <- cumsum(table(clusterID))[-nClusters]
	abline(h = boundaries, lwd = 2)

	par(mai=c(1.02,0.05,0.82,0.50))
	plot(expr, y = 1:length(expr), yaxs = 'i', pch = 19, xlab = "log2 expression", ylab = NA, yaxt = "n", cex = 0.5)

	par(oma = c(0, 0, 2, 0))
	mtext(tName, line = 0, outer = TRUE)
})
