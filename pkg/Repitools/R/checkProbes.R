checkProbes <- function(regionsTable, probesTable)
{
	hits <- as.data.frame(table(probesTable$name), stringsAsFactors = FALSE)
	colnames(hits) <- c("name", "number")
	regions <- GRanges(seqnames = regionsTable$chr, ranges = IRanges(start = regionsTable$start, end = regionsTable$end))
	probes <- GRanges(seqnames = probesTable$chr, ranges = IRanges(start = probesTable$start, end = probesTable$end), name = probesTable$name)
	overlaps <- findOverlaps(regions, probes)@matchMatrix
	indsPerRegion <- split(overlaps[, 2], overlaps[, 1])

	pdf("Probe Specificity.pdf", h = 8, w = 12)
	invisible(mapply(function(region, probesInRegion) {
		par(oma=c(4, 1, 1, 1))
		hitsMatch <- match(probesTable[probesInRegion, "name"], hits[, "name"])
		plot(x = probesTable[probesInRegion, "start"], y = hits[hitsMatch, "number"], xlim = c(regionsTable[region, "start"], regionsTable[region, "end"]), yaxt="n", type = 'p', pch = 19, cex = 0.5, main = regionsTable[region, "name"], xlab = "Position", ylab = "Probe Hits", col = c("red", "green", "blue"))
		ymax = max(hits[hitsMatch, "number"])
		if(ymax > 100)
		{
			ticks = c(1, seq(25, ymax, 25))
		} else if(ymax > 10)
		{
			ticks = c(1, seq(5, ymax, 5))
		} else {
			ticks = 1:ymax
		}
		axis(2, ticks, las = 2)
		mtext("TSS", 3, 1, at = regionsTable[region, "pos"], col="red")
		abline(v = regionsTable[region, "pos"], col = "red")
		abline(h = 1, lty = 2, col = "darkgrey")
		mtext(regionsTable[region, "chr"], 1, font = 2, at = 0.52, outer = TRUE)
	}, as.numeric(names(indsPerRegion)), indsPerRegion))
	dev.off()
}
