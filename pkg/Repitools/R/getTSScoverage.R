setGeneric("getTSScoverage", signature = "readsIPs", function(readsIPs, ...){standardGeneric("getTSScoverage")})

setMethod("getTSScoverage", "GRangesList", function(readsIPs, exptTypes, geneAnno, up = 50, down = 100, distType = c("percent", "base"), cvgResolution = 10, smoothingWidth = 2001)
{
	distType <- match.arg(distType)

	# Pool any replicates
	readsIPs <- split(unlist(readsIPs), rep(exptTypes, elementLengths(readsIPs)))
  inds <- split(1:length(exptTypes), unique(exptTypes))
  
  readsIPs <- lapply(inds, FUN=function(u) {
    gr <- GRanges()
    cat(":")
    for(i in 1:length(u)) {
      cat(".")
      gr <- c(gr,IPlist[[ u[i] ]])
    }
    return(gr)
  })
  gc()
  cat("\n")

	positions <- seq(-up, down, cvgResolution)
	if(distType == "percent")
		positionsLabels <- paste(positions, '%', sep = '')
	else
		positionsLabels <- positions
	widths <- geneAnno$end - geneAnno$start + 1
	geneAnno$tss <- ifelse(geneAnno$strand == '+', geneAnno$start, geneAnno$end)
	flank = floor(smoothingWidth / 2)

	# Get base positions to extract in 5' to 3' direction for every gene.
	regionBasesList <- lapply(1:nrow(geneAnno), function(geneIndex)
	{
		if(geneAnno[geneIndex, "strand"] == '+')
		{
			if(distType == "percent")
				region <- c(round(geneAnno[geneIndex, "start"] - up / 100 * widths[geneIndex] - flank), round(geneAnno[geneIndex, "start"] + (down / 100 * widths[geneIndex]) + flank))
			else
				region <- c(round(geneAnno[geneIndex, "start"] - up - flank), round(geneAnno[geneIndex, "start"] + down + flank))
		} else {
			if(distType == "percent")
				region <- c(round(geneAnno[geneIndex, "end"] + (up / 100 * widths[geneIndex]) + flank), round(geneAnno[geneIndex, "end"] - down / 100 * widths[geneIndex] - flank))
			else
				region <- c(round(geneAnno[geneIndex, "end"] + up + flank), round(geneAnno[geneIndex, "end"] - down - flank))
		}
		return(region)
	})

	# Get sampling positions in the extracted region for every gene.
	if(distType == "percent")
	{
		samplingPositions <- lapply(1:nrow(geneAnno), function(geneIndex)
		{
			round(flank + (0:(length(positions)-1) / (length(positions)-1)) * widths[geneIndex]) + 1
		})
	} else {
		samplingPositions <- round(flank + (0:(length(positions)-1) * cvgResolution)) + 1
	}

	# Get coverage sampling matrix for each IP mark.
	coverageMatrices <- lapply(readsIPs, function(readsIP)
	{
		# Scale all coverages to be as if there were 10 million reads in each GRanges object.
		# Make the coverage vectors be numeric. Much faster to access than Rle.
		coverageIP <- coverage(readsIP) * 10000000 / length(readsIP)
		
		TSScoverageList <- list()
		for(geneIndex in 1:nrow(geneAnno))
		{
			# Get smoothing region. Extra flanks so that we can smooth the edge of the region nicely.
			coverageGene <- tryCatch(coverageIP[[geneAnno[geneIndex, "chr"]]][regionBasesList[[geneIndex]][1]:regionBasesList[[geneIndex]][2]], error = function(errMsg) "Error getting all bases' coverage.")

			# If could get all bases, then smooth.
			if(class(coverageGene) == "Rle") {coverageGene <- runmean(coverageGene, smoothingWidth, endrule = "constant")} else {TSScoverageList[[geneIndex]] <- "Error"; next;}

			# Get the coverage at the sampling positions.
			if(distType == "percent")
				reprCoverage <- as.numeric(coverageGene[samplingPositions[[geneIndex]]])
			else
				reprCoverage <- as.numeric(coverageGene[samplingPositions])

			names(reprCoverage) <- positionsLabels
			TSScoverageList[[geneIndex]] <- reprCoverage	
		}
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

	rm(readsIPs, coverageIPs)
	gc()
	return(coverageMatrices)
})
