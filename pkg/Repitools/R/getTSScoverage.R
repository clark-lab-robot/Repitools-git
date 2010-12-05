setGeneric("getTSScoverage", signature = "readsIPs", function(readsIPs, ...){standardGeneric("getTSScoverage")})

setMethod("getTSScoverage", "GRangesList", function(readsIPs, seqLen = 300, exptTypes, geneAnno, up = 50, down = 100, distType = c("percent", "base"), cvgResolution = 10, smoothingWidth = 2001)
{
	distType <- match.arg(distType)
	readsIPs <- endoapply(readsIPs, resize, seqLen)

	# Pool any replicates.
	cat("Pooling any replicates.\n")
	typesInds <- split(1:length(exptTypes), exptTypes)

	readsIPs <- mapply(function(typeInds, typeName)
	{
		cat("Processing ", typeName, '\n')
		if(length(typeInds) == 1)
		{
			readsIPs[[typeInds]]
		} else {
			pooledReads <- GRanges()
			for(replicateIndex in 1:length(typeInds))
			{
				pooledReads <- c(pooledReads, readsIPs[[typeInds[replicateIndex]]])
			}
			return(pooledReads)
		}
	}, typesInds, names(typesInds))
  	gc()

	positions <- seq(-up, down, cvgResolution)
	if(distType == "percent")
		positionsLabels <- paste(positions, '%', sep = '')
	else
		positionsLabels <- positions
	widths <- geneAnno$end - geneAnno$start + 1
	geneAnno$tss <- ifelse(geneAnno$strand == '+', geneAnno$start, geneAnno$end)
	flank = floor(smoothingWidth / 2)
  
	cat("Getting positions of regions to extract.\n")
	# Get base positions to extract in 5' to 3' direction for every gene.
	regionBasesList <- mapply(function(start, end, strand, width)
	{
		if(strand == '+')
		{
			if(distType == "percent")
				region <- c(round(start - up / 100 * width - flank), round(start + (down / 100 * width) + flank))
			else
				region <- c(round(start - up - flank), round(start + down + flank))
		} else {
			if(distType == "percent")
				region <- c(round(end + up / 100 * width + flank), round(end - down / 100 * width - flank))
			else
				region <- c(round(end + up + flank), round(end - down - flank))
		}
		return(region)
	}, geneAnno$start, geneAnno$end, geneAnno$strand, widths, SIMPLIFY = FALSE)

	# Get sampling positions in the extracted region for every gene.
	if(distType == "percent")
	{
		samplingPositions <- lapply(widths, function(width)
		{
			round(flank + (0:(length(positions)-1) / (length(positions)-1)) * width) + 1
		})
	} else {
		samplingPositions <- round(flank + (0:(length(positions)-1) * cvgResolution)) + 1
	}

	cat("Getting coverage matrices.\n")
	# Get coverage sampling matrix for each IP mark.
	coverageMatrices <- lapply(readsIPs, function(readsIP)
	{
		# Scale all coverages to be as if there were 10 million reads in each GRanges object.
		# Make the coverage vectors be numeric. Much faster to access than Rle.
		coverageIP <- coverage(readsIP) * 10000000 / length(readsIP)
		
		TSScoverageList <- mapply(function(chr, regionBases, geneIndex)
		{
			# Get smoothing region. Extra flanks so that we can smooth the edge of the region nicely.
			coverageGene <- tryCatch(coverageIP[[chr]][regionBases[1]:regionBases[2]], error = function(errMsg) "Error getting all bases' coverage.")

			# If could get all bases, then smooth.
			if(class(coverageGene) == "Rle") {coverageGene <- runmean(coverageGene, smoothingWidth, endrule = "constant")} else {return("Error")}

			# Get the coverage at the sampling positions.
			if(distType == "percent")
				reprCoverage <- as.numeric(coverageGene[samplingPositions[[geneIndex]]])
			else
				reprCoverage <- as.numeric(coverageGene[samplingPositions])

			names(reprCoverage) <- positionsLabels
			return(reprCoverage)	
		}, geneAnno$chr, regionBasesList, 1:nrow(geneAnno), SIMPLIFY = FALSE)
		names(TSScoverageList) <- 1:nrow(geneAnno)
	
		# Find which genes didn't have all coverage available. Eliminate from the coverage list.
		whichEliminate <- which(sapply(TSScoverageList, function(geneCoverage) length(geneCoverage) == 1))
		if(length(whichEliminate) > 0)
		{
			TSScoverageList <- TSScoverageList[-whichEliminate]
		}
		TSScoverageTable <- do.call(rbind, TSScoverageList)
		rownames(TSScoverageTable) <- names(TSScoverageList)
		return(TSScoverageTable)
	})

	rm(readsIPs)
	gc()
	return(coverageMatrices)
})
