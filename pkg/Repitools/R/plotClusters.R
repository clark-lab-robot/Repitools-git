plotClusters <- function(clustersTable, posCol, scoreCol, scoreType, ylim)
{
	chrs <- paste("chr", c(1:22, "X", "Y"), sep = "")
	offsets <- seq(500, 75000, 1500)

	clustersTable <- clustersTable[clustersTable$cluster != 0, ]
	clustersTable <- clustersTable[order(factor(clustersTable$chr, levels = chrs)), ]
	clustersList <- split(clustersTable, factor(paste(clustersTable$chr, clustersTable$cluster), levels = unique(paste(clustersTable$chr, clustersTable$cluster))))

	invisible(lapply(clustersList, function(cluster)
	{
		cluster$modPos = cluster[, posCol]
		repeat
		{
			changed = rep(FALSE, nrow(cluster))
			for(index in 2:nrow(cluster))
			{
				if(cluster[index - 1, "modPos"] > cluster[index, "modPos"] - 250000)
				{		changed[index] = TRUE
						backIndex = index
						while(changed[backIndex - 1] == TRUE)
							backIndex = backIndex - 1
						cluster[index, "modPos"] <- cluster[index, "modPos"] + offsets[index - backIndex + 1]
				}
			}
			if(!any(changed))
				break
		}
		par(oma = c(3, 2, 6, 2), mar = c(5, 4, 5, 1))
		plot(cluster$modPos, cluster[, scoreCol], type = "h", xlab = "", ylab = scoreType, ylim = ylim, lwd = 2, xaxt = "n")
		title(cluster[1, "chr"], outer = TRUE)
		mtext("gene", side = 1, line = 5)
		axis(1, cluster$modPos, cluster$symbol, las = 2)
		axis(3, cluster$modPos, cluster[, posCol], las = 2, cex.axis = 0.8)
	}
	)
	)
	return()
}
