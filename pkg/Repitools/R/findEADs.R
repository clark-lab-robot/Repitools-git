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

findEADs <- function(statsTable, windowSize = 5, cutoff = 0.05, getFDRs = FALSE)
{
	statsTable <- statsTable[order(statsTable$chr, (statsTable$start + statsTable$end) / 2), ]

	tMed <- rep(NA, nrow(statsTable))
	t <- statsTable$t
	chrs <- statsTable$chr

	for(index in 1 : (nrow(statsTable) - (windowSize - 1)))
	{
		if(length(unique(chrs[index:(index + windowSize - 1)])) == 1) # all in window on same chromosome.
		{
			tMed[index + floor(windowSize / 2)] = median(t[index:(index + (windowSize - 1))])
		}
	}
	statsTable$tMed <- tMed

	statsTableRand <- statsTable
	whichNotNAs <- which(!is.na(tMed))
	statsTableRand[whichNotNAs, ] <- statsTable[sample(whichNotNAs), ]
	tMedRand = rep(NA, nrow(statsTableRand))
	tRand <- statsTableRand$t
	for(index in 1:(nrow(statsTableRand) - (windowSize - 1)))
	{
		if(length(unique(chrs[index:(index + windowSize - 1)])) == 1) # Just assume same chrs as real table.
		{
			tMedRand[index + floor(windowSize / 2)] = median(tRand[index:(index + (windowSize - 1))])
		}
	}
	
	tMedCutoffs <- seq(1, 10, 0.05)
	cutoffFDR <- data.frame(cutoff = tMedCutoffs, FDR = NA)
	for(cutoffIndex in 1:length(tMedCutoffs))
	{
		tMedCutoff = tMedCutoffs[cutoffIndex]
		EAD <- rep(0, length(t))
		EADrand <- rep(0, length(tRand))
		count = 1
		countRand = 1
		for(index in 1:(nrow(statsTable) - (windowSize - 1)))
		{	
			if(!any(is.na(tMed[index:(index + windowSize - 1)])) && t[index + floor(windowSize / 2)] > 0 && tMed[index + floor(windowSize / 2)] > tMedCutoff && length(which(tMed[index : (index + windowSize - 1)] > tMedCutoff)) >= 2)
			{
				indices <- index:(index + windowSize - 1)
				centreIndex <- indices[ceiling(windowSize / 2)]
				leftIndex <- indices[ceiling(windowSize / 2)]
				rightIndex <- indices[ceiling(windowSize / 2)]
				for(tsIndex in 1:(windowSize / 2))
				{
					if(t[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(t[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				}
				if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
				{
					EAD[leftIndex:rightIndex] = count
					count = count + 1
				}
			}

			if(!any(is.na(tMedRand[index:(index + windowSize - 1)])) && tRand[index + floor(windowSize / 2)] > 0 && tMedRand[index + floor(windowSize / 2)] > tMedCutoff && length(which(tMedRand[index : (index + windowSize - 1)] > tMedCutoff)) >= 2)
			{
				indices <- index:(index + windowSize - 1)
				centreIndex <- indices[ceiling(windowSize / 2)]
				leftIndex <- indices[ceiling(windowSize / 2)]
				rightIndex <- indices[ceiling(windowSize / 2)]
				for(tsIndex in 1:(windowSize / 2))
				{
					if(tRand[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
						leftIndex = leftIndex - 1
					if(tRand[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
						rightIndex = rightIndex + 1
				}
				if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
				{
					EADrand[leftIndex:rightIndex] = count
					countRand = countRand + 1
				}
			}
		}

		EAD <- .makeBlocks(EAD)
		EADrand <- .makeBlocks(EADrand)

		currFDR <- max(EADrand) / max(EAD) # max region number equivalent to total count 
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
	
	tMedCutoff <- cutoffFDR[cutIndex, "cutoff"]
	cat("Using the cutoff", tMedCutoff, "for a FDR of", cutoffFDR[cutIndex, "FDR"], '\n')

	EAD <- rep(0, length(t))
	count = 1
	for(index in 1:(nrow(statsTable) - (windowSize - 1)))
	{	
		if(!any(is.na(tMed[index:(index+(windowSize - 1))])) && t[index + floor(windowSize / 2)] > 0 && tMed[index + floor(windowSize / 2)] > tMedCutoff && t[index + floor(windowSize / 2)] > 0 && length(which(tMed[index : (index + windowSize - 1)] > tMedCutoff)) >= 2)
		{
			indices <- index:(index + windowSize - 1)
			centreIndex <- indices[ceiling(windowSize / 2)]
			leftIndex <- indices[ceiling(windowSize / 2)]
			rightIndex <- indices[ceiling(windowSize / 2)]
			for(tsIndex in 1:(windowSize / 2))
			{
				if(t[leftIndex - 1] > 0 && (statsTable[leftIndex, "chr"] == statsTable[leftIndex - 1, "chr"]))
					leftIndex = leftIndex - 1
				if(t[rightIndex + 1] > 0 && (statsTable[rightIndex, "chr"] == statsTable[rightIndex + 1, "chr"]))
					rightIndex = rightIndex + 1
			}
			if(leftIndex < rightIndex - 1) # Make sure at least 3 genes in region.
			{
				EAD[leftIndex:rightIndex] = count
				count = count + 1
			}
		}
	}

	EAD <- .makeBlocks(EAD)
	
	for(rowIndex in 2:(nrow(statsTable) - 1))
	{
		if(EAD[rowIndex - 1] == 0 && EAD[rowIndex] > 0)
		{
			extIndex = rowIndex - 1
			EADcode = EAD[rowIndex]
			while(t[extIndex] > 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
			{
				EAD[extIndex] <- EADcode
				extIndex <- extIndex - 1
			}
		}
	
		if(EAD[rowIndex] > 0 && EAD[rowIndex + 1] == 0)
		{
			extIndex = rowIndex + 1
			EADcode = EAD[rowIndex]
			while(t[extIndex] > 0 && (statsTable[rowIndex, "chr"] == statsTable[extIndex, "chr"]))
			{
				EAD[extIndex] <- EADcode
				extIndex <- extIndex + 1	
			}
		}
	}

	statsTable$EAD <- EAD

	if(getFDRs == TRUE)
	{
		return(list(table = statsTable, FDRs = cutoffFDR))
	} else {
		return(statsTable)
	}
}


