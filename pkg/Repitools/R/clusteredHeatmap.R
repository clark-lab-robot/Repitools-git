setGeneric("clusteredHeatmap", signature = "readsIPs", function(readsIPs, ...){standardGeneric("clusteredHeatmap")})

setMethod("clusteredHeatmap", "GRangesList", function(readsIPs, exptTypes, geneAnno, upPct = 50, downPct = 100, cvgResolution = 10, smoothingWidth = 2001, capQuantile = 0.95, nClusters = 5, expr, cols = colorpanel(64, "blue", "white", "red"), tName = "Clustered Enrichment")
{
	# Pool any replicates found.
	typeCounts <- table(exptTypes)
	if(any(typeCounts > 1))
	{
		oldOrder <- unique(exptTypes)
		repTypes <- names(typeCounts)[typeCounts > 1]
		pooledIPs <- lapply(repTypes, function(repType)
		{
			whichReps <- which(exptTypes == repType)
			do.call(c, unname(as.list(readsIPs[whichReps])))
		})
		names(pooledIPs) <- repTypes
		uniqueIdxs <- exptTypes %in% names(typeCounts)[typeCounts == 1]
		uniqueIPs <- as.list(readsIPs[uniqueIdxs])
		names(uniqueIPs) <- exptTypes[uniqueIdxs]
		readsIPs <- c(pooledIPs, uniqueIPs)
		readsIPs <- GRangesList(readsIPs[oldOrder])
	}

	# Scale all coverages to be as if there were 10 million reads in each GRanges object.
	coverageIPs <- lapply(readsIPs, function(readsIP) coverage(readsIP) * 10000000 / length(readsIP))

	pctPositions <- seq(-upPct, downPct, cvgResolution)
	widths <- geneAnno$end - geneAnno$start + 1
	geneAnno$tss <- ifelse(geneAnno$strand == '+', geneAnno$start, geneAnno$end)
	
	flank = floor(smoothingWidth / 2)
	# Get matrix for each IP mark.
	coverageMatrices <- lapply(coverageIPs, function(coverageIP)
	{
		TSScoverageList <- list()
		for(geneIndex in 1:nrow(geneAnno))
		{
			# Get smoothing region. Extra flanks so that we can smooth the edge of the region nicely.
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
		names(TSScoverageList) <- 1:nrow(geneAnno)	
	
		# Find which gene didn't have all coverage available. Eliminate from the coverage list.
		whichEliminate <- which(sapply(TSScoverageList, function(geneCoverage) length(geneCoverage) == 1))
		if(length(whichEliminate) > 0)
		{
			TSScoverageList <- TSScoverageList[-whichEliminate]
		}
		TSScoverageTable <- do.call(rbind, TSScoverageList)
		rownames(TSScoverageTable) <- names(TSScoverageList)
		return(TSScoverageTable)
	})

	# Now only keep in the expression and width vectors, the genes with all coverage sample positions.
	expr <- expr[as.integer(rownames(coverageMatrices[[1]]))]		
	widths <- widths[as.integer(rownames(coverageMatrices[[1]]))]
	geneAnno <- geneAnno[as.integer(rownames(coverageMatrices[[1]])), ]

	# Find the maximum score allowable, then make any scores bigger be the maximum score.
	allScores <- do.call(cbind, coverageMatrices)
	maxIntensity = quantile(allScores, capQuantile)
	coverageMatrices <- lapply(coverageMatrices, function(coverageMatrix)
	{
		coverageMatrix[coverageMatrix > maxIntensity] = maxIntensity
		return(coverageMatrix)
	})

	# Do the k-means clustering for all marks together.
	clusterID <- kmeans(allScores, nClusters, iter.max = 100)$cluster
	
	# Get median expression for each cluster. Find ascending order.
	medClusterExpr <- tapply(expr, factor(clusterID), median)
	exprOrder <- order(medClusterExpr)
	
	# Get order by gene width next.
	sortedIndices <- unlist(lapply(exprOrder, function(cluster)
	{
		whichGenes <- which(clusterID == cluster)
		widthOrder <- order(widths[whichGenes])
		orderedGenes <- whichGenes[widthOrder]
		return(orderedGenes)
	}))

	# Re-arrange the ChIP and expression data and vector of gene widths.
	coverageMatrices <- lapply(coverageMatrices, function(coverageMatrix) coverageMatrix[sortedIndices, ])
	expr <- expr[sortedIndices]
	widths <- widths[sortedIndices]

	# Plot heatmap.
	layout(rbind(1:(length(coverageMatrices) + 3)), widths=c(1, rep(3, length(coverageMatrices)), 1, 1))
	par(mai=c(1.02,0.50,0.82,0.05))
  	
	par(oma = c(0, 0, 0, 0))
	heatmapRange <- range(allScores)
	nBins = length(cols)
	image(y=seq(1/nBins/2, 1-(1/nBins/2), 1/nBins), z=rbind(1:nBins), col = cols, axes = FALSE, xlab = "Signal Intensity", ylab = NA)
	axis(2, at=(0:nBins)/nBins, labels=format(seq(heatmapRange[1], heatmapRange[2], length.out = nBins + 1), digits = 2))

	par(mai=c(1.02,0.05,0.82,0.05))
	mapply(function(coverageMatrix, exptName)
	{
		image(pctPositions, 1:nrow(geneAnno), t(coverageMatrix), xlab="Position (%) relative to TSS", yaxt = "n", ylab = "Gene", col = cols, main = exptName)
		
		# Add lines delimiting the cluster boundaries.
		boundaries <- cumsum(table(clusterID)[exprOrder])[-nClusters]
		abline(h = boundaries, lwd = 2)
	}, coverageMatrices, unique(exptTypes))

	par(mai=c(1.02,0.05,0.82,0.50))
	plot(expr, y = 1:length(expr), yaxs = 'i', pch = 19, xlab = "log2 expression", ylab = NA, yaxt = "n", cex = 0.5)
	plot(widths, y = 1:length(widths), yaxs = 'i', pch = 19, xlab = "Gene widths", ylab = NA, yaxt = "n", cex = 0.5)	

	par(oma = c(0, 0, 2, 0))
	mtext(tName, line = 0, outer = TRUE)
})
