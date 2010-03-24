doSeqStats <- function(reads, whichInputs, whichControl, whichTreat, minCount, blockSize, blocksTable)
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
	keepRegions <- rowSums(inputRegionsReads) > minCount[1]
	rownames(inputRegionsReads) <- rownames(inputRegions)
	inputRegionsReads <- inputRegionsReads[keepRegions, ]
	inputRegions <- inputRegions[keepRegions, ]

	if(is.null(blocksTable))
	{
		enrichedIntervals <-  lapply(chrLengths, function(chrLength) {segs <- data.frame(start = seq(1, chrLength - blockSize[2] + 1, blockSize[2]), end = seq(blockSize[2], chrLength, blockSize[2])); segs[nrow(segs), "end"] <- chrLength; return(segs)})
		for(index in 1:length(enrichedIntervals))
		{
			enrichedIntervals[[index]] <- cbind(chr = seqnames(Hsapiens)[index], enrichedIntervals[[index]], stringsAsFactors = FALSE)
		}
		enrichedRegions <- do.call(rbind, enrichedIntervals)
		rm(enrichedIntervals)
	} else {
		enrichedRegions <- data.frame(name = blocksTable$name, chr = blocksTable$chr, start = NA, end = NA, stringsAsFactors = FALSE)
		enrichedRegions$start <- ifelse(blocksTable$strand == "+", blocksTable$start - 1000, blocksTable$end - 1000)
		enrichedRegions$end <- ifelse(blocksTable$strand == "+", blocksTable$start + 1000, blocksTable$end + 1000)
	 }

	enrichedRegionsReads <- annotationBlocksCounts(reads[c(whichControl, whichTreat)], enrichedRegions, seqLen = 300, verbose = FALSE)
	keepRegions <- rowSums(enrichedRegionsReads) > minCount[2]
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
	groups[whichTreat] <- "T"
	groups[whichControl] <- "C"
	groups <- groups[which(nchar(groups) == 1)]
	treatInSubset <- which(groups == "T")
	controlInSubset <- which(groups == "C")
	analysis <- mapply(function(tableForCN, cnLevel) {
                                                          cnSF <- numeric()
	                                                  cnSF[controlInSubset] <- 1
	                                                  cnSF[treatInSubset] <- cnLevel
   	                                                  dge <- DGEList(counts=tableForCN, group=factor(groups), lib.size = readsTotals * cnSF)
  	                                                  dge <- estimateCommonDisp(dge) 
  	                                                  results <- exactTest(dge, pair=c("C","T"))
							  return(cbind(results$table, TreatmentCN = cnLevel)) # Relative to control being x1.
                                                         }, tables, CNlevels, SIMPLIFY = FALSE
                          )
	resultsTable <- data.frame()
	for(CNindex in 1:length(analysis))
	{
		resultsTable <- rbind(resultsTable, cbind(enrichedRegions[indicesByCN[[CNindex]], ], analysis[[CNindex]], enrichedRegionsReads[indicesByCN[[CNindex]], ]))
	}
	resultsTable <- resultsTable[order(resultsTable$chr, resultsTable$start), ]
	resultsTable$adj.p.val <- p.adjust(resultsTable[, "p.value"], method = "BH")
	resultsTable$zScore <- sign(resultsTable$logFC)*abs(qnorm(resultsTable$adj.p.val/2))
	return(resultsTable)
}
