setGeneric("getTSScoverage", signature = "readsIPs", function(readsIPs, ...){standardGeneric("getTSScoverage")})

setMethod("getTSScoverage", "GRangesList", function(readsIPs, seqLen = 300, exptTypes, geneAnno, up = 50, down = 100, distType = c("percent", "base"), cvgResolution = 10, smoothingWidth = 2001)
{
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
	geneAnno <- geneAnno[chrOrder, ]
	annoGR <- GRanges(geneAnno$chr, IRanges(geneAnno$start, geneAnno$end), geneAnno$strand, order = geneAnno$original)

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

		if(region[1] < 0 || region[2] > chrLengths[[chr]]) # Region overlaps past edge of chromosome.
			return(numeric())
		else
			return(region)
	}, geneAnno$chr, geneAnno$start, geneAnno$end, geneAnno$strand, widths, SIMPLIFY = FALSE)

	# Eliminate regions past edges.
	whichOutside <- which(sapply(regionBasesList, function(regionBases) length(regionBases) == 0))
	annoGR <- annoGR[-whichOutside, ]
	widths <- widths[-whichOutside]
	regionBasesList <- regionBasesList[-whichOutside]
	whichReverse <- which(strand(annoGR) == '-')

	# Create windows to use in Views().
	windowsGR <- GRanges(as.character(seqnames(annoGR)), IRanges(sapply(regionBasesList, function(bases) bases[1]), sapply(regionBasesList, function(bases) bases[2])))
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
		coverageGenes <- viewApply(unlist(Views(coverageIP, windowsRL)), function(aView) Rle(as.numeric(aView)))
		coverageGenes[whichReverse] <- lapply(coverageGenes[whichReverse], rev)

		# Running mean smoothing.
		coverageGenes <- RleList(lapply(coverageGenes, function(coverageGene) runmean(coverageGene, smoothingWidth, endrule = "constant")))

		# Get the coverage at the sampling positions.
		TSScoverageTable=do.call(rbind, viewApply(Views(coverageGenes, RangesList(lapply(samplingPositions, function(x) IRanges(x, width = 1)))), as.numeric)@listData)
		colnames(TSScoverageTable) <- positionsLabels	
		rownames(TSScoverageTable) <- elementMetadata(annoGR)[, "order"]
		return(TSScoverageTable)
	})

	rm(readsIPs)
	gc()
	return(coverageMatrices)
})
