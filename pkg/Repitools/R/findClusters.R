.makeClusters <- function(scoresOnChr, wSize, cut, count = FALSE) {

    mergeOverlaps <- function (query, subject) {
        temp <- IRanges(slice(coverage(c(query, subject)), lower=1))
        temp[countOverlaps(temp, query)>0]
    }

    clustersByColumn <- function(scoreByColumn) {
                medianByColumn <- rollmedian(scoreByColumn, k = wSize, na.pad = TRUE)
                mbc.cut <- medianByColumn > cut
                mbc.cut[is.na(mbc.cut)] <- FALSE
                mbc.2cut <- rollmean(as.integer(mbc.cut), k = wSize, na.pad = TRUE) > (2 / wSize)
                mbc.2cut[is.na(mbc.2cut)] <- FALSE
                sbc.cut <- scoreByColumn > 0

                # Window extension.
                mediansCut <- mbc.cut & mbc.2cut
                extend <- sbc.cut
                
                extended <- mergeOverlaps(IRanges(mediansCut), IRanges(extend))
                
		# There must be at least 3 genes in an EAD / LRES region, by Saul's definition.
                extended <- unique(extended[width(extended) >= 3])
		
                if (count == TRUE)
		{
			return(length(extended))
		} else {
			as.numeric(coverage(extended, width = length(scoreByColumn)))
		}
            }
    
    clusters <- numeric()
    if(count == TRUE)
    {
    	clusters <- apply(scoresOnChr, 2, clustersByColumn)
    } else {
    	clusters <- clustersByColumn(scoresOnChr[, 1])
    }
    
    return(clusters)
}

findClusters <- function(statsTable, posCol, scoreCol, windowSize = 5, cutoff = 0.05, trend = c("down", "up"), nPermutations = 100, getFDRs = FALSE, verbose = TRUE)
{
    require(IRanges)
    require(zoo)
    trend <- match.arg(trend)
    
    # Do check in case users pass in some rows with NA scores. For biological meaningfulness.
    statsTable <- statsTable[!is.na(statsTable[, scoreCol]), ]

    scoreVect <- numeric()
    if(trend == "down")
    {
        scoreVect = -statsTable[, scoreCol]
    } else {
        scoreVect = statsTable[, scoreCol]
    }

    perms <- 1:nPermutations
    scoresByChr <- split(data.frame(scoreVect, sapply(perms, function(x) scoreVect[sample(nrow(statsTable))])), statsTable$chr)

    scoreMedCutoffs <- seq(1, 10, 0.05)
    
    clustersAtCutoffs <- lapply(scoreMedCutoffs, function(aCutoff) {
                                                                   regionCounts <- colSums(do.call(rbind, lapply(scoresByChr, .makeClusters, windowSize, aCutoff, TRUE)))
								   return(c(realCount = regionCounts[1], randCount = round(mean(regionCounts[2:length(regionCounts)]))))
                                                                   })


    allFDR <- sapply(clustersAtCutoffs, function(clustersPerCutoff) {
                  currFDR <- clustersPerCutoff[2] / clustersPerCutoff[1]
                  if((is.nan(currFDR))) # 0 / 0 case
                      currFDR = 0
                  return(currFDR)
                  })

    FDRtable <- data.frame(cutoff = scoreMedCutoffs, FDR = allFDR)
    chosenCutoff <- scoreMedCutoffs[match(TRUE, allFDR < cutoff)]
    
    if(verbose == TRUE)
        cat("Using the cutoff", chosenCutoff, "for a FDR of <", cutoff, "\n")

    clusterCol <- unlist(lapply(scoresByChr, function(scoresOnChr) .makeClusters(scoresOnChr, windowSize, chosenCutoff)))

    clusterNum = 1
    for(index in 2:length(clusterCol))
    {
    	if(clusterCol[index - 1] == 0 && clusterCol[index] > 0)
	{
		clusterCol[index] = clusterNum
	} else if(clusterCol[index - 1] > 0 && clusterCol[index] > 0)
	{
		clusterCol[index] = clusterNum
	} else if((clusterCol[index - 1] > 0 && clusterCol[index] == 0))
	{
		clusterNum = clusterNum + 1
	}
    }

    statsTable$cluster <- clusterCol

    if(getFDRs == TRUE)
    {
    	return(list(table = statsTable, FDRs = FDRtable))
    } else {
    		return(statsTable)
    }
    
}
