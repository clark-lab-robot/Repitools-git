.getAvgQScore <- function(cycTable){
                                   totalCounts <- sum(cycTable[, "Count"])
				   totalScore = 0
				   for(index in 1:nrow(cycTable))
				   	totalScore = totalScore + (cycTable[index, "Score"] * cycTable[index, "Count"])
				   avg <- totalScore / totalCounts
				   return(unlist(avg))
                                  }

# paths : absoulte paths to folders of lanes for one experiment.
# names : names of experiments.
genQC <- function(paths, exptName, laneNames)
{
	require(ShortRead)

	tablesListList <- lapply(paths, function(path){
                                        # First, the average quality score over the cycles.
				    	QAfastq <- qa(path, "fastq.gz", "fastq")
					sbCyc <- split(QAfastq[["perCycle"]][["quality"]], QAfastq[["perCycle"]][["quality"]][, "Cycle"])
					avgScores <- sapply(sbCyc, .getAvgQScore)

					# Percentage base calls over the cycles.
					aligned <- readAligned(path, "uniq.map.gz", "Bowtie")
					fastq <- readFastq(path, "fastq.gz")
					freqBasesAll <- apply(alphabetByCycle(sread(fastq), c("G", "A", "T", "C", "N")), 2, function(aCol) aCol / sum(aCol))
					freqBasesAligned <- apply(alphabetByCycle(sread(aligned), c("G", "A", "T", "C", "N")), 2, function(aCol) aCol / sum(aCol))
					colnames(freqBasesAll) <- 1:ncol(freqBasesAll)
					colnames(freqBasesAligned) <- 1:ncol(freqBasesAligned)

					# Mismatch matrix over the cycles.
					mmColumn <- aligned@alignData@data$mismatch
					
					# Make complements of negative strand sequences to get back to the original FASTQ sequence.
					whichMapNeg <- which(strand(aligned) == '-')
					mmColumn[whichMapNeg] <- gsub("G", "c", mmColumn[whichMapNeg])
					mmColumn[whichMapNeg] <- gsub("C", "g", mmColumn[whichMapNeg])
					mmColumn[whichMapNeg] <- gsub("A", "t", mmColumn[whichMapNeg])
					mmColumn[whichMapNeg] <- gsub("T", "a", mmColumn[whichMapNeg])
					mmColumn[whichMapNeg] <- toupper(mmColumn[whichMapNeg])

					mmMatrix <- do.call(rbind, strsplit(unlist(strsplit(mmColumn, ",", fixed = TRUE)), ":", fixed = TRUE))
					
					mmTable <- table(mmMatrix[,2], factor(mmMatrix[,1], levels = as.character(0:35)))
					mmTable <- apply(mmTable, 2, function(cycCol) cycCol / sum(cycCol))
					colnames(mmTable) <- as.numeric(colnames(mmTable)) + 1
					gc()

					return(list(avgQual = avgScores, freqAll = freqBasesAll, freqAligned = freqBasesAligned, mmTable = mmTable))
				    }
                             )
	names(tablesListList) <- laneNames
	
	pdf(paste("QC ", exptName, ".pdf", sep = ''), h = 8.3, w = 11.7)
	avgQs <- sapply(tablesListList, function(tList) tList[[1]])
	matplot(avgQs, type = "l", lty = 1, lwd = 2, col = 1:8, main = paste("Average Quality Over Cycles for", exptName), xlab = "Cycle", ylab = "Quality", ylim = c(15, 40))
	legend("topright", legend = names(tablesListList), lty = 1, lwd = 2, col = 1:8)
	layout(matrix(c(1:8), ncol = 4, byrow = TRUE))
	invisible(mapply(function(tList, name){
					     	matplot(y = t(tList[[2]]),  main = paste("All Bases For", name), xlab = "Cycle", ylab = "Base", type = "l", lty = 1, lwd = 2, xlim = c(1, 36), ylim = c(0.00, 0.40), new = TRUE)
						abline(h = 0.25, col = "red")
						legend("topright", legend = c("G", "A", "T", "C", "N"), lty = 1, lwd = 2, col = 1:5, cex = 0.5)
					     }, tablesListList, names(tablesListList)
			)
                 )
	layout(matrix(c(1:8), ncol = 4, byrow = TRUE))
	invisible(mapply(function(tList, name){
				     		matplot(y = t(tList[[3]]),  main = paste("Aligned Bases For", name), xlab = "Cycle", ylab = "Base", type = "l", lty = 1, lwd = 2, xlim = c(1, 36), ylim = c(0.00, 0.40))
						abline(h = 0.25, col = "red")
						legend("topright", legend = c("G", "A", "T", "C", "N"), lty = 1, lwd = 2, col = 1:5, cex = 0.5)
					     }, tablesListList, names(tablesListList)
			)
		 )
	par(oma = c(1, 1, 1, 1))
	layout(matrix(c(1:4), ncol = 2, byrow = TRUE), widths=c(5,1))
	invisible(mapply(function(tList, name){
						par(mai = c(1,1.2,0.2,0.5))
					     	matplot(y = t(tList[[4]]),  main = paste("Mismatches For", name), xlab = "Cycle", ylab = "Mismatch", type = "l", lty = rep(c(1, 2), each = 8), lwd = 2, xlim = c(1, 36), ylim = c(0.00, 1.00), col = c("black", "darkblue", "purple", "darkgreen", "cyan", "lightgreen", "red", "yellow"))
						abline(v = 1:36, col = c("lightgrey", "darkgrey"))
						par(mai = c(0,0,0.2,0))
						plot.new()
						legend("topleft", title = "Reference > Read", legend = rownames(tList[[4]]), lty = rep(c(1, 2), each = 8), lwd = 2, col = c("black", "darkblue", "purple", "darkgreen", "cyan", "lightgreen", "red", "yellow"))
					     }, tablesListList, names(tablesListList)
			)
		)
	dev.off()

	return(tablesListList)
}
