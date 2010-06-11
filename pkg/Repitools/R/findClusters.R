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

findClusters <- function(statsTable, posCol, scoreCol, windowSize = 5, cutoff = 0.05, trend = c("down", "up"), getFDRs = FALSE)
{
	trend <- match.arg(trend)
	
	statsTable <- statsTable[order(statsTable$chr, statsTable[, posCol]), ]

	scoreMed <- rep(NA, nrow(statsTable))
	score <- statsTable[, scoreCol]
	chrs <- statsTable$chr

	for(index in 1 : (nrow(statsTable) - (windowSize - 1)))
	{
		if(length(unique(chrs[index:(index + windowSize - 1)])) == 1) # all in window on same chromosome.
		{
			scoreMed[index + floor(windowSize / 2)] = median(score[index:(index + (windowSize - 1))])
		}
	}
	statsTable$scoreMed <- scoreMed

	statsTableRand <- statsTable
	whichNotNAs <- which(!is.na(scoreMed))
	statsTableRand[whichNotNAs, ] <- statsTable[sample(whichNotNAs), ]
	scoreMedRand = rep(NA, nrow(statsTableRand))
	scoreRand <- statsTableRand[, scoreCol]
	for(index in 1:(nrow(statsTableRand) - (windowSize - 1)))
	{
		if(length(unique(chrs[index:(index + windowSize - 1)])) == 1) # Just assume same chrs as real table.
		{
			scoreMedRand[index + floor(windowSize / 2)] = median(scoreRand[index:(index + (windowSize - 1))])
		}
	}
	
	if(trend == "down")
	{
		scoreMedCutoffs <- seq(-1, -10, -0.05)	
	} else {
		scoreMedCutoffs <- seq(1, 10, 0.05)
	}
	cutoffFDR <- data.frame(cutoff = scoreMedCutoffs, FDR = NA)
	for(cutoffIndex in 1:length(scoreMedCutoffs))
	{
		scoreMedCutoff = scoreMedCutoffs[cutoffIndex]
		cluster <- rep(0, length(score))
		clusterRand <- rep(0, length(scoreRand))
		count = 1
		countRand = 1
		for(index in 1:(nrow(statsTable) - (windowSize - 1)))
		{
			if(trend == "down")
			{
				if(!any(is.na(scoreMed[index:(index + windowSize - 1)])) && score[index + floor(windowSize / 2)] < 0 && scoreMed[index + floor(windowSize / 2)] < scoreMedCutoff && length(which(scoreMed[index : (index + windowSize - 1)] < scoreMedCutoff)) >= 2)
				{
					indices <- index:(index + windowSize - 1)
					centreIndex <- indices[ceiling(windowSize / 2)]
					leftIndex <- indices[ceiling(windowSize / 2)]
					rightIndex <- indices[ceiling(windowSize / 2)]
					for(scoresIndex in 1:(windowSize / 2))
					{
						if(score[leftIndex - 1] < 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
							leftIndex = leftIndex - 1
						if(score[rightIndex + 1] < 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
							rightIndex = rightIndex + 1
					}
					if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
					{
						cluster[leftIndex:rightIndex] = count
						count = count + 1
					}
				}

				if(!any(is.na(scoreMedRand[index:(index + windowSize - 1)])) && scoreRand[index + floor(windowSize / 2)] < 0 && scoreMedRand[index + floor(windowSize / 2)] < scoreMedCutoff && length(which(scoreMedRand[index : (index + windowSize - 1)] < scoreMedCutoff)) >= 2)
				{
					indices <- index:(index + windowSize - 1)
					centreIndex <- indices[ceiling(windowSize / 2)]
					leftIndex <- indices[ceiling(windowSize / 2)]
					rightIndex <- indices[ceiling(windowSize / 2)]
					for(scoresIndex in 1:(windowSize / 2))
					{
						if(scoreRand[leftIndex - 1] < 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
							leftIndex = leftIndex - 1
						if(scoreRand[rightIndex + 1] < 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
							rightIndex = rightIndex + 1
					}
					if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
					{
						clusterRand[leftIndex:rightIndex] = count
						countRand = countRand + 1
					}
				}
			} else {
				if(!any(is.na(scoreMed[index:(index + windowSize - 1)])) && score[index + floor(windowSize / 2)] > 0 && scoreMed[index + floor(windowSize / 2)] > scoreMedCutoff && length(which(scoreMed[index : (index + windowSize - 1)] > scoreMedCutoff)) >= 2)
				{
					indices <- index:(index + windowSize - 1)
					centreIndex <- indices[ceiling(windowSize / 2)]
					leftIndex <- indices[ceiling(windowSize / 2)]
					rightIndex <- indices[ceiling(windowSize / 2)]
					for(scoresIndex in 1:(windowSize / 2))
					{
						if(score[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
							leftIndex = leftIndex - 1
						if(score[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
							rightIndex = rightIndex + 1
					}
					if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
					{
						cluster[leftIndex:rightIndex] = count
						count = count + 1
					}
				}

				if(!any(is.na(scoreMedRand[index:(index + windowSize - 1)])) && scoreRand[index + floor(windowSize / 2)] > 0 && scoreMedRand[index + floor(windowSize / 2)] > scoreMedCutoff && length(which(scoreMedRand[index : (index + windowSize - 1)] > scoreMedCutoff)) >= 2)
				{
					indices <- index:(index + windowSize - 1)
					centreIndex <- indices[ceiling(windowSize / 2)]
					leftIndex <- indices[ceiling(windowSize / 2)]
					rightIndex <- indices[ceiling(windowSize / 2)]
					for(scoresIndex in 1:(windowSize / 2))
					{
						if(scoreRand[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
							leftIndex = leftIndex - 1
						if(scoreRand[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
							rightIndex = rightIndex + 1
					}
					if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
					{
						clusterRand[leftIndex:rightIndex] = count
						countRand = countRand + 1
					}
				}
			}
		}

		cluster <- .makeBlocks(cluster)
		clusterRand <- .makeBlocks(clusterRand)

		currFDR <- max(clusterRand) / max(cluster) # max region number equivalent to total count 
		if((is.nan(currFDR))) # 0 / 0 case
			currFDR = 0
		
		cutoffFDR[cutoffIndex, "FDR"] <- currFDR
	}
	
	cutIndex <- numeric()
	for(index in 1:length(cutoffFDR$FDR))
	{
		if(cutoffFDR$FDR[index] < cutoff)
		{
			cutIndex = index
			break
		}
	}
	
	scoreMedCutoff <- cutoffFDR[cutIndex, "cutoff"]
	cat("Using the cutoff", scoreMedCutoff, "for a FDR of", cutoffFDR[cutIndex, "FDR"], '\n')

	cluster <- rep(0, length(score))
	count = 1
	for(index in 1:(nrow(statsTable) - (windowSize - 1)))
	{
		if(trend == "down")
		{
			if(!any(is.na(scoreMed[index:(index+(windowSize - 1))])) && score[index + floor(windowSize / 2)] < 0 && scoreMed[index + floor(windowSize / 2)] < scoreMedCutoff && length(which(scoreMed[index : (index + windowSize - 1)] < scoreMedCutoff)) >= 2)
			{
				indices <- index:(index + windowSize - 1)
				centreIndex <- indices[ceiling(windowSize / 2)]
				leftIndex <- indices[ceiling(windowSize / 2)]
				rightIndex <- indices[ceiling(windowSize / 2)]
				for(tsIndex in 1:(windowSize / 2))
				{
					if(score[leftIndex - 1] < 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(score[rightIndex + 1] < 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				}
				if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
				{
					cluster[leftIndex:rightIndex] = count
					count = count + 1
				}
			}
		} else {
			if(!any(is.na(scoreMed[index:(index+(windowSize - 1))])) && score[index + floor(windowSize / 2)] > 0 && scoreMed[index + floor(windowSize / 2)] > scoreMedCutoff && length(which(scoreMed[index : (index + windowSize - 1)] > scoreMedCutoff)) >= 2)
			{
				indices <- index:(index + windowSize - 1)
				centreIndex <- indices[ceiling(windowSize / 2)]
				leftIndex <- indices[ceiling(windowSize / 2)]
				rightIndex <- indices[ceiling(windowSize / 2)]
				for(tsIndex in 1:(windowSize / 2))
				{
					if(score[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(score[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				}
				if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
				{
					cluster[leftIndex:rightIndex] = count
					count = count + 1
				}
			}
		}
	}

	cluster <- .makeBlocks(cluster)
	
	for(rowIndex in 2:(nrow(statsTable) - 1))
	{
		if(trend == "down")
		{
			if(cluster[rowIndex - 1] == 0 && cluster[rowIndex] > 0)
			{
				extIndex = rowIndex - 1
				clusterCode = cluster[rowIndex]
				while(score[extIndex] < 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
				{
					cluster[extIndex] <- clusterCode
					extIndex <- extIndex - 1
				}
			}
		
			if(cluster[rowIndex] > 0 && cluster[rowIndex + 1] == 0)
			{
				extIndex = rowIndex + 1
				clusterCode = cluster[rowIndex]
				while(score[extIndex] < 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
				{
					cluster[extIndex] <- clusterCode
					extIndex <- extIndex + 1	
				}
			}
		}
		else {
			if(cluster[rowIndex - 1] == 0 && cluster[rowIndex] > 0)
			{
				extIndex = rowIndex - 1
				clusterCode = cluster[rowIndex]
				while(score[extIndex] > 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
				{
					cluster[extIndex] <- clusterCode
					extIndex <- extIndex - 1
				}
			}
		
			if(cluster[rowIndex] > 0 && cluster[rowIndex + 1] == 0)
			{
				extIndex = rowIndex + 1
				clusterCode = cluster[rowIndex]
				while(score[extIndex] > 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
				{
					cluster[extIndex] <- clusterCode
					extIndex <- extIndex + 1	
				}
			}
		}
	}

	statsTable$cluster <- cluster

	if(getFDRs == TRUE)
	{
		return(list(table = statsTable, FDRs = cutoffFDR))
	} else {
		return(statsTable)
	}
}
