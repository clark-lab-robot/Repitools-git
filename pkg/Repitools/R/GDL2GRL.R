setGeneric(".GDL2GRL", signature = "gdl", function(gdl, ...)
                              {standardGeneric(".GDL2GRL")})

setMethod(".GDL2GRL", "GenomeDataList", function(gdl)
{
	grl <- seqapply(gdl, function(gd) {
		gr <- do.call(c, lapply(names(gd), function(chr) {
			pos <- gd[[chr]]
			starts <- c(pos[[names(pos)[1]]], pos[[names(pos)[2]]])
			GRanges(chr, IRanges(starts, width = 1), 
			rep(names(pos), elementLengths(pos)))
		}))
	})
	
	names(grl) <- names(gdl)
	elementMetadata(grl) <- elementMetadata(gdl)
	metadata(grl) <- metadata(gdl)
	grl
})
