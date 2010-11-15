setGeneric("doSeqStats", function(reads, ...){standardGeneric("doSeqStats")})

setMethod("doSeqStats", "GenomeDataList", function(reads, seqLen, whichInputs, whichControl, whichTreat, minCount, blockSize, CNlevels, blocksTable, bpUp, bpDown)
{
	require(edgeR)
	require(BSgenome.Hsapiens.UCSC.hg18)

	if(is.null(blocksTable))
	{
		inputRegions <- as.data.frame(genomeBlocks(Hsapiens, 1:25, blockSize[1]))[, 1:3]
		enrichedRegions <- as.data.frame(genomeBlocks(Hsapiens, 1:25, blockSize[2]))[, 1:3]
		colnames(inputRegions)[1] <- "chr"
		colnames(enrichedRegions)[1] <- "chr"
	} else {
		if("strand" %in% colnames(blocksTable))
		{
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, strand = blocksTable$strand, symbol = blocksTable$symbol, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			enrichedRegions$start <- ifelse(blocksTable$strand == "+", blocksTable$start - bpUp, blocksTable$end - bpDown)
			enrichedRegions$end <- ifelse(blocksTable$strand == "+", blocksTable$start + bpDown, blocksTable$end + bpUp)	
		} else {
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			positions <- round((blocksTable$start + blocksTable$end) / 2)
			enrichedRegions$start <- positions - bpUp
			enrichedRegions$end <- positions + bpUp
		}
			intervals <- mapply(function(chr, start, end) {
			segs <- data.frame(start = seq(start - blockSize[1], end + blockSize[1] + 1, blockSize[1]), end = seq(start - 1, end + 2 * blockSize[1], blockSize[1]), stringsAsFactors = FALSE)
			segs <- cbind(chr, segs, stringsAsFactors = FALSE)
			return(segs)
			}, enrichedRegions$chr, enrichedRegions$start, enrichedRegions$end, SIMPLIFY = FALSE)

			inputRegions <- do.call(rbind, intervals)
			inputRegions <- inputRegions[!duplicated(inputRegions), ]
			rm(intervals)

	 }

	inputRegionsReads <- annotationBlocksCounts(reads[whichInputs], inputRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(inputRegionsReads) >= minCount[1]
	rownames(inputRegionsReads) <- rownames(inputRegions)
	inputRegionsReads <- inputRegionsReads[keepRegions, ]
	inputRegions <- inputRegions[keepRegions, ]

	enrichedRegionsReads <- annotationBlocksCounts(reads[c(whichControl, whichTreat)], enrichedRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(enrichedRegionsReads) >= minCount[2]
	rownames(enrichedRegionsReads) <- rownames(enrichedRegions)
	enrichedRegionsReads <- enrichedRegionsReads[keepRegions, ]
	enrichedRegions <- enrichedRegions[keepRegions, ]

	CNs <- getCN(list(inputs = inputRegions, enriched = enrichedRegions), inputRegionsReads)
	CNblocks <- cut(CNs[, 2], breaks = quantile(CNs[,2], p=(0:CNlevels)/CNlevels), include.lowest=TRUE)
	indicesByCN <- split(1:nrow(enrichedRegions), CNblocks)
	CNlevels <- sapply(indicesByCN, function(indices) {median(CNs[indices, 2])})
	tables <- lapply(indicesByCN, function(indices) enrichedRegionsReads[indices, ])

	readsTotals <- colSums(enrichedRegionsReads)

	groups <- character()
	groups[(length(whichControl)+1):(length(whichControl) + length(whichTreat))] <- "T"
	groups[1:length(whichControl)] <- "C"
	treatInSubset <- which(groups == "T")
	controlInSubset <- which(groups == "C")
	analysis <- mapply(function(tableForCN, cnLevel) {
                                                          cnSF <- numeric()
	                                                  cnSF[controlInSubset] <- 1
	                                                  cnSF[treatInSubset] <- cnLevel
   	                                                  dge <- DGEList(counts = tableForCN, group = factor(groups), lib.size = readsTotals * cnSF, genes = rownames(tableForCN))
  	                                                  dge <- estimateCommonDisp(dge) 
  	                                                  results <- exactTest(dge, pair = c("C","T"))
							  return(cbind(results$table, TreatmentCN = cnLevel)) # Relative to control being x1.
                                                         }, tables, CNlevels, SIMPLIFY = FALSE
                          )

	normTo10mFs <- 1/(laneCounts(reads[c(whichControl, whichTreat)]) / 10000000)
	ERRTo10M <- matrix(nrow = nrow(enrichedRegionsReads), ncol = ncol(enrichedRegionsReads), dimnames = dimnames(enrichedRegionsReads))
	colnames(ERRTo10M) <- paste(colnames(ERRTo10M), "Per 10 Million Reads")
	for(colIndex in 1:ncol(ERRTo10M))
		ERRTo10M[, colIndex] <- round(enrichedRegionsReads[, colIndex] *  normTo10mFs[colIndex])

	resultsTable <- data.frame()
	for(CNindex in 1:length(analysis))
	{
		analysisSmall <- analysis[[CNindex]]
		analysisExpanded <- matrix(nrow = length(indicesByCN[[CNindex]]), ncol = 4, dimnames = list(rownames(enrichedRegions[indicesByCN[[CNindex]], ]), colnames(analysisSmall)))
		map <- match(rownames(analysisSmall), rownames(analysisExpanded))
		analysisExpanded[map, ] <- as.matrix(analysisSmall)

		resultsTable <- rbind(resultsTable, cbind(enrichedRegions[indicesByCN[[CNindex]], ], analysisExpanded, ERRTo10M[indicesByCN[[CNindex]], ]))
	}
	resultsTable <- resultsTable[order(resultsTable$chr, resultsTable$start), ]
	
	resultsTable$adj.p.val <- p.adjust(resultsTable[, "p.value"], method = "BH")
	resultsTable$zScore <- sign(resultsTable$logFC) * abs(qnorm(resultsTable$p.value / 2))
	resultsTable$zeroReads <- apply(resultsTable[, grep("Reads", colnames(resultsTable))], 1, function(aRow) length(which(aRow == 0)))
	resultsTable$totalReads <- rowSums(resultsTable[, grep("Million", colnames(resultsTable))], na.rm = TRUE)

	return(resultsTable)
})

setMethod("doSeqStats", "GRangesList", function(reads, seqLen, whichInputs, whichControl, whichTreat, minCount, blockSize, CNlevels, blocksTable, bpUp, bpDown)
{
	require(edgeR)
	require(BSgenome.Hsapiens.UCSC.hg18)

	if(is.null(blocksTable))
	{
		inputRegions <- as.data.frame(genomeBlocks(Hsapiens, 1:25, blockSize[1]))[, 1:3]
		enrichedRegions <- as.data.frame(genomeBlocks(Hsapiens, 1:25, blockSize[2]))[, 1:3]
		colnames(inputRegions)[1] <- "chr"
		colnames(enrichedRegions)[1] <- "chr"
	} else {
		if("strand" %in% colnames(blocksTable))
		{
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, strand = blocksTable$strand, symbol = blocksTable$symbol, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			enrichedRegions$start <- ifelse(blocksTable$strand == "+", blocksTable$start - bpUp, blocksTable$end - bpDown)
			enrichedRegions$end <- ifelse(blocksTable$strand == "+", blocksTable$start + bpDown, blocksTable$end + bpUp)	
		} else {
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			positions <- round((blocksTable$start + blocksTable$end) / 2)
			enrichedRegions$start <- positions - bpUp
			enrichedRegions$end <- positions + bpUp
		}
			intervals <- mapply(function(chr, start, end) {
			segs <- data.frame(start = seq(start - blockSize[1], end + blockSize[1] + 1, blockSize[1]), end = seq(start - 1, end + 2 * blockSize[1], blockSize[1]), stringsAsFactors = FALSE)
			segs <- cbind(chr, segs, stringsAsFactors = FALSE)
			return(segs)
			}, enrichedRegions$chr, enrichedRegions$start, enrichedRegions$end, SIMPLIFY = FALSE)

			inputRegions <- do.call(rbind, intervals)
			inputRegions <- inputRegions[!duplicated(inputRegions), ]
			rm(intervals)

	 }

	inputRegionsReads <- annotationBlocksCounts(reads[whichInputs], inputRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(inputRegionsReads) >= minCount[1]
	rownames(inputRegionsReads) <- rownames(inputRegions)
	inputRegionsReads <- inputRegionsReads[keepRegions, ]
	inputRegions <- inputRegions[keepRegions, ]

	enrichedRegionsReads <- annotationBlocksCounts(reads[c(whichControl, whichTreat)], enrichedRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(enrichedRegionsReads) >= minCount[2]
	rownames(enrichedRegionsReads) <- rownames(enrichedRegions)
	enrichedRegionsReads <- enrichedRegionsReads[keepRegions, ]
	enrichedRegions <- enrichedRegions[keepRegions, ]

	CNs <- getCN(list(inputs = inputRegions, enriched = enrichedRegions), inputRegionsReads)
	CNblocks <- cut(CNs[, 2], breaks = quantile(CNs[,2], p=(0:CNlevels)/CNlevels), include.lowest=TRUE)
	indicesByCN <- split(1:nrow(enrichedRegions), CNblocks)
	CNlevels <- sapply(indicesByCN, function(indices) {median(CNs[indices, 2])})
	tables <- lapply(indicesByCN, function(indices) enrichedRegionsReads[indices, ])

	readsTotals <- colSums(enrichedRegionsReads)

	groups <- character()
	groups[(length(whichControl)+1):(length(whichControl) + length(whichTreat))] <- "T"
	groups[1:length(whichControl)] <- "C"
	treatInSubset <- which(groups == "T")
	controlInSubset <- which(groups == "C")
	analysis <- mapply(function(tableForCN, cnLevel) {
                                                          cnSF <- numeric()
	                                                  cnSF[controlInSubset] <- 1
	                                                  cnSF[treatInSubset] <- cnLevel
   	                                                  dge <- DGEList(counts = tableForCN, group = factor(groups), lib.size = readsTotals * cnSF, genes = rownames(tableForCN))
  	                                                  dge <- estimateCommonDisp(dge) 
  	                                                  results <- exactTest(dge, pair = c("C","T"))
							  return(cbind(results$table, TreatmentCN = cnLevel)) # Relative to control being x1.
                                                         }, tables, CNlevels, SIMPLIFY = FALSE
                          )

	normTo10mFs <- 1/(elementLengths(reads[c(whichControl, whichTreat)]) / 10000000)
	ERRTo10M <- matrix(nrow = nrow(enrichedRegionsReads), ncol = ncol(enrichedRegionsReads), dimnames = dimnames(enrichedRegionsReads))
	colnames(ERRTo10M) <- paste(colnames(ERRTo10M), "Per 10 Million Reads")
	for(colIndex in 1:ncol(ERRTo10M))
		ERRTo10M[, colIndex] <- round(enrichedRegionsReads[, colIndex] *  normTo10mFs[colIndex])

	resultsTable <- data.frame()
	for(CNindex in 1:length(analysis))
	{
		analysisSmall <- analysis[[CNindex]]
		analysisExpanded <- matrix(nrow = length(indicesByCN[[CNindex]]), ncol = 4, dimnames = list(rownames(enrichedRegions[indicesByCN[[CNindex]], ]), colnames(analysisSmall)))
		map <- match(rownames(analysisSmall), rownames(analysisExpanded))
		analysisExpanded[map, ] <- as.matrix(analysisSmall)

		resultsTable <- rbind(resultsTable, cbind(enrichedRegions[indicesByCN[[CNindex]], ], analysisExpanded, ERRTo10M[indicesByCN[[CNindex]], ]))
	}
	resultsTable <- resultsTable[order(resultsTable$chr, resultsTable$start), ]
	
	resultsTable$adj.p.val <- p.adjust(resultsTable[, "p.value"], method = "BH")
	resultsTable$zScore <- sign(resultsTable$logFC) * abs(qnorm(resultsTable$p.value / 2))
	resultsTable$zeroReads <- apply(resultsTable[, grep("Reads", colnames(resultsTable))], 1, function(aRow) length(which(aRow == 0)))
	resultsTable$totalReads <- rowSums(resultsTable[, grep("Million", colnames(resultsTable))], na.rm = TRUE)

	return(resultsTable)
})

