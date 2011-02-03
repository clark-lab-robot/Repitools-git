setGeneric("getTSScoverage", signature = "readsIPs", function(readsIPs, ...){standardGeneric("getTSScoverage")})

setMethod("getTSScoverage", "GRangesList", function(readsIPs, seqLen = 300, exptTypes, geneAnno, up = 50, down = 100, distType = c("percent", "base"), cvgResolution = 10, smoothingWidth = 2001)
{
	if(!all(c("chr", "start", "end", "strand", "name") %in% colnames(geneAnno)))
		stop("Some of the mandatory columns of geneAnno are not present. See documentation for requirements.")

	distType <- match.arg(distType)
	readsIPs <- endoapply(readsIPs, resize, seqLen)
	chrLengths <- seqlengths(readsIPs)

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

	# Prepare for extracting Views by sorting in chromosomal order.
	
	geneAnno$original <- 1:nrow(geneAnno)
	chrOrder <- order(geneAnno$chr)
	originalOrder <- order(chrOrder)
	geneAnno <- geneAnno[chrOrder, ]
	annoGR <- GRanges(geneAnno$chr, IRanges(geneAnno$start, geneAnno$end), geneAnno$strand, order = geneAnno$original)
	whichReverse <- which(strand(annoGR) == '-')

	positions <- seq(-up, down, cvgResolution)
	if(distType == "percent")
		positionsLabels <- paste(positions, '%', sep = '')
	else
		positionsLabels <- positions
	widths <- geneAnno$end - geneAnno$start + 1
	flank = floor(smoothingWidth / 2)
  
	cat("Getting positions of regions to extract.\n")
	# Get base positions to extract for every gene.
	regionBasesList <- mapply(function(chr, start, end, strand, width)
	{
		if(strand == '+')
		{
			if(distType == "percent")
				region <- c(round(start - up / 100 * width - flank), round(start + (down / 100 * width) + flank))
			else
				region <- c(round(start - up - flank), round(start + down + flank))
		} else {
			if(distType == "percent")
				region <- c(round(end - down / 100 * width - flank), round(end + up / 100 * width + flank))
			else
				region <- c(round(end - down - flank), round(end + up + flank))
		}

		if(region[1] < 1)
		{
			flankStart = rep(0, abs(region[1] - 1))
			region[1] = 1
		} else {
			flankStart = numeric()
		}
		if(region[2] > chrLengths[[chr]])
		{
			flankEnd = rep(0, region[2] - chrLengths[[chr]])
			region[2] = chrLengths[[chr]]
		} else {
			flankEnd = numeric()
		}
			return(list(region, flankStart, flankEnd))
	}, geneAnno$chr, geneAnno$start, geneAnno$end, geneAnno$strand, widths, SIMPLIFY = FALSE)

	# Create windows to use in Views().
	windowsGR <- GRanges(as.character(seqnames(annoGR)), IRanges(sapply(regionBasesList, function(regionData) regionData[[1]][1]), sapply(regionBasesList, function(regionData) regionData[[1]][2])))
	windowsRL <- as(windowsGR, "RangesList")

	# Get sampling positions in the extracted region for every gene.
	if(distType == "percent")
	{
		samplingPositions <- lapply(widths, function(width)
		{
			round(flank + (0:(length(positions)-1) / (length(positions)-1)) * width) + 1
		})
	} else {
		samplingPositions <- rep(list(round(flank + (0:(length(positions)-1) * cvgResolution)) + 1), length(windowsGR))
	}

	cat("Getting coverage matrices.\n")
	# Get coverage sampling matrix for each IP mark.
	coverageMatrices <- lapply(readsIPs, function(readsIP)
	{
		# Scale all coverages to be as if there were 10 million reads in each GRanges object.
		coverageIP <- coverage(readsIP) * 10000000 / length(readsIP)
		coverageIP <- coverageIP[names(windowsRL)]

		# Get smoothing region. Extra flanks so that we can smooth the edge of the region nicely.
		coverageGenes <- viewApply(unlist(Views(coverageIP, windowsRL)), function(coverage) coverage)
		coverageGenes[whichReverse] <- lapply(coverageGenes[whichReverse], rev)

		# Running mean smoothing.
		coverageGenes <- lapply(coverageGenes, function(coverageGene) runmean(coverageGene, smoothingWidth, endrule = "constant"))

		# Tack on the zeros for the regions past the ends.
		coverageGenes <- mapply(function(start, coverage, end)
		                       {
				       		if(strand == '+')
		      	               			c(start, coverage, end)
						else
							c(end, coverage, start)
		      		       }, regionBasesList[[2]], coverageGenes, regionBasesList[[3]], geneAnno$strand, SIMPLIFY = FALSE)
		
		# Get the coverage at the sampling positions.
		TSScoverageTable <- t(mapply(function(genePositions, geneCoverage)
		{
			as.numeric(geneCoverage[genePositions])
		}, samplingPositions, coverageGenes))
		TSScoverageTable <- TSScoverageTable[originalOrder, ]
		
		colnames(TSScoverageTable) <- positionsLabels	
		rownames(TSScoverageTable) <- geneAnno[, "name"]
		return(TSScoverageTable)
	})

	rm(readsIPs)
	gc()
	return(coverageMatrices)
})
