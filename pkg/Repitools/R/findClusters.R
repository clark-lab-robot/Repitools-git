.getMedians <- function(sTable, sCol, wSize)
{
	chrs <- sTable$chr
	scoreMed <- rep(NA, nrow(sTable))
	score <- sTable[, sCol]
	for(index in 1 : (nrow(sTable) - (wSize - 1)))
	{
		if(length(unique(chrs[index:(index + wSize - 1)])) == 1) # all in window on same chromosome.
		{
			scoreMed[index + floor(wSize / 2)] = median(score[index:(index + (wSize - 1))])
		}
	}
	return(scoreMed)
}

.makeClusters <- function(aTable, sCol, medCol, wSize, cut, cType)
{
	score <- aTable[, sCol]
	scoreMed <- aTable[, medCol]
	cluster <- rep(0, length(score))
	count = 1
	for(index in 1:(nrow(aTable) - (wSize - 1)))
	{
		if(!any(is.na(scoreMed[index:(index + wSize - 1)])) && ifelse(cType == "down", score[index + floor(wSize / 2)] < 0, score[index + floor(wSize / 2)] > 0) && ifelse(cType == "down", scoreMed[index + floor(wSize / 2)] < cut, scoreMed[index + floor(wSize / 2)] > cut) && length(which(scoreMed[index : (index + wSize - 1)] > cut)) >= 2)
		{
			indices <- index:(index + wSize - 1)
			leftIndex <- indices[ceiling(wSize / 2)]
			rightIndex <- indices[ceiling(wSize / 2)]
			for(scoresIndex in 1:(wSize / 2))
			{
				if(cType == "down")
				{
					if(score[leftIndex - 1]  < 0 && (aTable[leftIndex, "chr"] == aTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(score[rightIndex + 1] < 0 && (aTable[rightIndex, "chr"] == aTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				} else {
					if(score[leftIndex - 1] > 0 && (aTable[leftIndex, "chr"] == aTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(score[rightIndex + 1] > 0 && (aTable[rightIndex, "chr"] == aTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				}
			}
			if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
			{
				cluster[leftIndex:rightIndex] = count
				count = count + 1
			}
		}
	}
	return(cluster)
}

.makeBlocks <- function(regionCol)
{
	clusterNum = 1
	for(index in 2:length(regionCol))
	{
		if(regionCol[index - 1] == 0 && regionCol[index] > 0)
		{
			regionCol[index] = clusterNum
		} else if(regionCol[index - 1] > 0 && regionCol[index] > 0)
		{
			regionCol[index] = clusterNum
		} else if((regionCol[index - 1] > 0 && regionCol[index] == 0))
		{
			clusterNum = clusterNum + 1
		}
	}
	return(regionCol)	
}

.extendLeft <- function(cluster, sTable, sCol, cType)
{
	score <- sTable[, sCol]
	for(index in 2:(nrow(sTable) - 1))
	{
		if(cluster[index - 1] == 0 && cluster[index] > 0)
		{
			extIndex = index - 1
			clusterCode = cluster[index]
			while(extIndex > 0 && ifelse(cType == "down", score[extIndex] < 0, score[extIndex] > 0) && (sTable[index, "chr"] == sTable[extIndex, "chr"]))
			{
				cluster[extIndex] <- clusterCode
				extIndex <- extIndex - 1	
			}
		}
	}
	return(cluster)
}

.extendRight <- function(cluster, sTable, sCol, cType)
{
	score <- sTable[, sCol]
	for(index in 2:(nrow(sTable) - 1))
	{
		if(cluster[index] > 0 && cluster[index + 1] == 0)
		{
			extIndex = index + 1
			clusterCode = cluster[index]
			while(extIndex <= length(cluster) && ifelse(cType == "down", score[extIndex] < 0, score[extIndex] > 0) && (sTable[index, "chr"] == sTable[extIndex, "chr"]))
			{
				cluster[extIndex] <- clusterCode
				extIndex <- extIndex + 1	
			}
		}
	}
	return(cluster)
}

findClusters <- function(statsTable, posCol, scoreCol, windowSize = 5, cutoff = 0.05, trend = c("down", "up"), nPermutations = 10, getFDRs = FALSE, verbose = TRUE)
{
	require(multicore)
	trend <- match.arg(trend)

	# Do check in case users pass in some rows with NA scores. For biological meaningfulness.
	statsTable <- statsTable[!is.na(statsTable[, scoreCol]), ]
	
	statsTable <- statsTable[order(statsTable$chr, statsTable[, posCol]), ]
	statsTable$scoreMed <- .getMedians(statsTable, scoreCol, windowSize)

	perms <- 1:nPermutations
	if(verbose == TRUE)
		cat("Generating", nPermutations, "sets of random scores and calculating medians.")
	randTables <- mclapply(perms, function(iteration)
	{
		statsTableRand <- statsTable[, -ncol(statsTable)]
		statsTableRand[, scoreCol] <- statsTable[sample(1:nrow(statsTable)), scoreCol]
		statsTableRand$scoreMedRand <- .getMedians(statsTableRand, scoreCol, windowSize)
		return(statsTableRand)
	}, mc.cores = 5
	)
	
	if(trend == "down")
	{
		scoreMedCutoffs <- seq(-1, -10, -0.05)	
	} else {
		scoreMedCutoffs <- seq(1, 10, 0.05)
	}

	clustersAtCutoffs <- mclapply(scoreMedCutoffs, function(aCutoff)
	{
		if(verbose == TRUE)
			cat("Processing real scores at", aCutoff, '\n')
		return(.makeBlocks(.makeClusters(statsTable, scoreCol, ncol(statsTable), windowSize, aCutoff, cType = trend)))
	}, mc.cores = 6
	)
		
	randTablesClusters <- mclapply(randTables, function(rTable)
	{
		if(verbose == TRUE)
			cat("Processing a set of random scores.\n")
		clusts <- mclapply(scoreMedCutoffs, function(aCutoff)
		{
			if(verbose == TRUE)	
				cat("Processing randomised scores at", aCutoff, '\n')
			return(.makeBlocks(.makeClusters(rTable, scoreCol, ncol(rTable), windowSize, aCutoff, cType = trend)))
		}, mc.cores = 2
		)
		if(verbose == TRUE)
			cat("Finishing processing a set of random scores.\n")
		return(clusts)

	}, mc.cores = 4
	)

	FDRtables <- mclapply(randTablesClusters, function(randTableClusters)
	{
		cutoffFDR <- data.frame(cutoff = scoreMedCutoffs, FDR = NA)
		fdr <- mapply(function(real, rand) {
		                                    currFDR <- max(rand) / max(real)
		                                    if((is.nan(currFDR))) # 0 / 0 case
		                                    	currFDR = 0
		                                    return(currFDR)
					           }, clustersAtCutoffs, randTableClusters)
	        cutoffFDR$FDR = fdr
		return(cutoffFDR)
		if(verbose == TRUE)
			cat("Finished FDRs for one random table.")
	}, mc.cores = 4
	)
	
	cutoffs <- sapply(FDRtables, function(aTable)
	{
		for(index in 1:nrow(aTable))
		{
			if(aTable$FDR[index] < cutoff)
			{
				break
			}
		}
		return(scoreMedCutoffs[index])
	}
	)

	chosenCutoff <- mean(cutoffs)
	
	if(verbose == TRUE)
		cat("Using the cutoff", chosenCutoff, "for a FDR of <", cutoff, "\n")

	cluster <- .makeBlocks(.makeClusters(statsTable, scoreCol, ncol(statsTable), windowSize, chosenCutoff, cType = trend))
	
	cluster <- .extendLeft(cluster, statsTable, scoreCol, cType = trend)
	cluster <- .makeBlocks(.extendRight(cluster, statsTable, scoreCol, cType = trend))

	statsTable$cluster <- cluster

	if(getFDRs == TRUE)
	{
		return(list(table = statsTable, FDRs = FDRtables))
	} else {
		return(statsTable)
	}
}
