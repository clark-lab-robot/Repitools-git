doSeqStats <- function(reads, whichInputs, whichControl, whichTreat, minCount, blockSize, blocksTable = NULL)
{
	require(edgeR)
	require(BSgenome.Hsapiens.UCSC.hg18)

	chrNames <- paste("chr", c(1:22, "X", "Y"), sep = "")
	chrLengths <- seqlengths(Hsapiens)[chrNames]
	if(is.null(blocksTable))
	{
		intervals <- lapply(chrLengths, function(chrLength) {
		                                                    segs <- data.frame(start = seq(1, chrLength - blockSize[1] + 1, blockSize[1]), end = seq(blockSize[1], chrLength, blockSize[1]))
		                                                    segs[nrow(segs), "end"] <- chrLength; return(segs)
								    }
			   )
		for(index in 1:length(intervals))
		{
			intervals[[index]] <- cbind(chr = seqnames(Hsapiens)[index], intervals[[index]])
		}
	} else {
		intervals <- mapply(function(chr, start, end) {
							      	segs <- data.frame(start = seq(start - blockSize[1], end + blockSize[1] + 1, blockSize[1]), end = seq(start - 1, end + 2 * blockSize[1], blockSize[1]), stringsAsFactors = FALSE)
								segs <- cbind(chr, segs, stringsAsFactors = FALSE)
								return(segs)
							      },
						              blocksTable$chr, blocksTable$start, blocksTable$end, SIMPLIFY = FALSE)
	}

	inputRegions <- do.call(rbind, intervals)
	inputRegions <- inputRegions[!duplicated(inputRegions), ]
	rm(intervals)

	inputRegionsReads <- annotationBlocksCounts(reads[whichInputs], inputRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(inputRegionsReads) >= minCount[1]
	rownames(inputRegionsReads) <- rownames(inputRegions)
	inputRegionsReads <- inputRegionsReads[keepRegions, ]
	inputRegions <- inputRegions[keepRegions, ]

	if(is.null(blocksTable))
	{
		enrichedIntervals <-  lapply(chrLengths, function(chrLength) {
                                                                             segs <- data.frame(start = seq(1, chrLength - blockSize[2] + 1, blockSize[2]), end = seq(blockSize[2], chrLength, blockSize[2]))
                                                                             return(segs)
                                                                             }
                                            )
		for(index in 1:length(enrichedIntervals))
		{
			enrichedIntervals[[index]] <- cbind(chr = seqnames(Hsapiens)[index], enrichedIntervals[[index]], stringsAsFactors = FALSE)
		}	
		enrichedRegions <- do.call(rbind, enrichedIntervals)
		rm(enrichedIntervals)
	} else {
		if("strand" %in% colnames(blocksTable))
		{
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, strand = blocksTable$strand, symbol = blocksTable$symbol, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			enrichedRegions$start <- ifelse(blocksTable$strand == "+", blocksTable$start - 1000, blocksTable$end - 1000)
			enrichedRegions$end <- ifelse(blocksTable$strand == "+", blocksTable$start + 1000, blocksTable$end + 1000)
		} else {
			enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, start = NA, end = NA, featureStart = blocksTable$start, featureEnd = blocksTable$end, stringsAsFactors = FALSE)
			positions <- round((blocksTable$start + blocksTable$end) / 2)
			enrichedRegions$start <- positions - 500
			enrichedRegions$end <- positions + 500
		}
	 }

	enrichedRegionsReads <- annotationBlocksCounts(reads[c(whichControl, whichTreat)], enrichedRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(enrichedRegionsReads) >= minCount[2]
	rownames(enrichedRegionsReads) <- rownames(enrichedRegions)
	enrichedRegionsReads <- enrichedRegionsReads[keepRegions, ]
	enrichedRegions <- enrichedRegions[keepRegions, ]

	CNs <- getCN(list(inputs = inputRegions, enriched = enrichedRegions), inputRegionsReads)
	CNblocks <- cut(CNs[, 2], breaks = quantile(CNs[,2], p=(0:10)/10), include.lowest=TRUE)
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

	normTo10mFs <- 1/(laneCounts(reads[c(whichControl, whichTreat)])/10000000)
	ERRTo10M <- matrix(nrow = nrow(enrichedRegionsReads), ncol = ncol(enrichedRegionsReads), dimnames = dimnames(enrichedRegionsReads))
	colnames(ERRTo10M) <- paste(colnames(ERRTo10M), "Per 10 Million Reads")
	for(colIndex in 1:ncol(ERRTo10M))
		ERRTo10M[, colIndex] <- round(enrichedRegionsReads[, colIndex] *  normTo10mFs[colIndex])

	resultsTable <- data.frame()
	for(CNindex in 1:length(analysis))
	{
		analysisSmall <- analysis[[CNindex]]
		analysisExpanded <- matrix(nrow = nrow(enrichedRegions[indicesByCN[[CNindex]], ]), ncol = 4, dimnames = list(rownames(enrichedRegions[indicesByCN[[CNindex]], ]), colnames(analysisSmall)))
		map <- match(rownames(analysisSmall), rownames(analysisExpanded))
		analysisExpanded[map, ] <- as.matrix(analysisSmall)

		resultsTable <- rbind(resultsTable, cbind(enrichedRegions[indicesByCN[[CNindex]], ], analysisExpanded, ERRTo10M[indicesByCN[[CNindex]], ]))
	}
	resultsTable <- resultsTable[order(resultsTable$chr, resultsTable$start), ]
	
	resultsTable$adj.p.val <- p.adjust(resultsTable[, "p.value"], method = "BH")
	resultsTable$zScore <- sign(resultsTable$logFC)*abs(qnorm(resultsTable$p.value/2))
	resultsTable$zeroReads <- apply(resultsTable[, grep("Reads", colnames(resultsTable))], 1, function(aRow) length(which(aRow == 0)))
	resultsTable$totalReads <- rowSums(resultsTable[, grep("Million", colnames(resultsTable))], na.rm = TRUE)

	return(resultsTable)
}
