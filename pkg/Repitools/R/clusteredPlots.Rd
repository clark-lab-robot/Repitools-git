setGeneric("clusteredBinplots", signature = "coverageList", function(coverageList, ...){standardGeneric("clusteredBinplots")})

setMethod("clusteredBinplots", "list", function(coverageList, geneAnno, capQuantile = 0.95, nClusters = 5, expr, plotType = c("line", "heatmap"), cols = NULL, logWidth = TRUE, tName = "Clustered Enrichment")
{
	plotType <- match.arg(plotType)

	if(is.null(cols) == TRUE)
	{
		require(gplots)
		if(plotType=="line")
		{
			cols <- colorpanel(length(coverageList), "blue", "green", "red")
		} else {
		  	cols <- colorpanel(64, "blue", "white", "red")
		}
	}

	widths <- geneAnno$end - geneAnno$start + 1

	# Now only keep in the expression and width vectors, the genes with all coverage sample positions.
	expr <- expr[as.integer(rownames(coverageList[[1]]))]		
	widths <- widths[as.integer(rownames(coverageList[[1]]))]
	geneAnno <- geneAnno[as.integer(rownames(coverageList[[1]])), ]

	# Find the maximum score allowable, then make any scores bigger be the maximum score.
	allScores <- do.call(cbind, coverageList)
	maxIntensity = quantile(allScores, capQuantile)
	coverageList <- lapply(coverageList, function(coverageMatrix)
	{
		coverageMatrix[coverageMatrix > maxIntensity] = maxIntensity
		return(coverageMatrix)
	})

	# Do the k-means clustering for all marks together.
	set.seed(100)
	clusterID <- kmeans(allScores, nClusters, iter.max = 100)$cluster
	
	# Get median expression for each cluster. Find ascending order.
	medClusterExpr <- tapply(expr, factor(clusterID), median)
	exprOrder <- order(medClusterExpr)

	# Get x-axis positions and labels.
	posLabels <- colnames(coverageList[[1]])
	positions <- as.integer(gsub('%', '', posLabels)) # Get raw position if labels have percentage signs.

	if(plotType == "line")
	{
		# Get the cluster average profile within each mark for each cluster.
				
		profilesList <- lapply(coverageList, function(coverageMatrix)
		{
			lapply(1:nClusters, function(index)
			{
				whichGenes <- which(clusterID == index)
				colMeans(coverageMatrix[whichGenes, ])	
			})
		})

		# Group each cluster from all epigenetic marks.

		profilesByCluster <- lapply(1:nClusters, function(index)
		{
			lapply(profilesList, function(markList)
			{
				markList[[exprOrder[index]]]
			})
		})

		# Plot the lineplots.
		layout(rbind(1:2), widths=c(3, 1))
		invisible(mapply(function(clustProfile, expr)
		{
			par(mai=c(1,2,0.8,0.5))
			plot.new()
			plot.window(xlim = c(positions[1], positions[length(positions)]), ylim = c(0, maxIntensity * 1.1), xaxs = "i", yaxs = "i")
			title(main = paste("Within Cluster Coverage (Median Expression : ", round(expr, 2), ")", sep = ''))
			title(xlab = "Position")
			title(ylab = "Smoothed Coverage")
			axis(1, positions, labels = posLabels)
			axis(2)
			mapply(function(clustMarkProfile, aCol)
			{
				lines(positions, clustMarkProfile, col = aCol)
			}, clustProfile, cols)

			plot.new()
			par(mai=c(1,0.1,0.2,0.1))
			legend("topright", legend = names(profilesList), title = "Mark", col = cols, lwd = 2)
		}, profilesByCluster, medClusterExpr[exprOrder]))

	} else { # Plot a heatmap

		# Get order by gene width next.
		sortedIndices <- unlist(lapply(exprOrder, function(cluster)
		{
			whichGenes <- which(clusterID == cluster)
			widthOrder <- order(widths[whichGenes])
			orderedGenes <- whichGenes[widthOrder]
			return(orderedGenes)
		}))

		# Re-arrange the ChIP and expression data and vector of gene widths.
		coverageList <- lapply(coverageList, function(coverageMatrix) coverageMatrix[sortedIndices, ])
		expr <- expr[sortedIndices]
		widths <- widths[sortedIndices]

		# Plot heatmap.
		layout(rbind(1:(length(coverageList) + 3)), widths=c(1, rep(3, length(coverageList)), 2, 1))
		par(mai=c(1.02,0.50,0.82,0.05))
  	
		par(oma = c(0, 0, 0, 0))
		heatmapRange <- range(allScores)
		nBins = length(cols)
		image(y=seq(1/nBins/2, 1-(1/nBins/2), 1/nBins), z=rbind(1:nBins), col = cols, axes = FALSE, xlab = "Signal Intensity", ylab = NA)
		axis(2, at=(0:nBins)/nBins, labels=format(seq(heatmapRange[1], heatmapRange[2], length.out = nBins + 1), digits = 2))

		par(mai=c(1.02,0.05,0.82,0.05))
		mapply(function(coverageMatrix, exptName)
		{
			image(positions, 1:nrow(geneAnno), t(coverageMatrix), xlab="Position relative to TSS", xaxt = "n", yaxt = "n", ylab = "Gene", col = cols, main = exptName)
			axis(1, positions, labels = posLabels)		

			# Add lines delimiting the cluster boundaries.
			boundaries <- cumsum(table(clusterID)[exprOrder])[-nClusters]
			abline(h = boundaries, lwd = 2)
		}, coverageList, names(coverageList))

		if(logWidth == TRUE) widths <- log(widths)
		par(mai=c(1.02,0.05,0.82,0.50))
		plot(expr, y = 1:length(expr), yaxs = 'i', pch = 19, xlab = "log2 expression", ylab = NA, yaxt = "n", cex = 0.5)
		plot(widths, y = 1:length(widths), yaxs = 'i', pch = 19, xlab = ifelse(logWidth == TRUE, "log(Gene widths)", "Gene widths") , ylab = NA, yaxt = "n", cex = 0.5)	

		par(oma = c(0, 0, 2, 0))
		mtext(tName, line = 0, outer = TRUE)
	}

	returnList <- split(geneAnno$name, factor(clusterID, levels = exprOrder))
	names(returnList) <- medClusterExpr[exprOrder]
	invisible(returnList)
})
