setGeneric("laneCounts", function(cs){standardGeneric("laneCounts")})

setMethod(laneCounts, "GenomeDataList", function(cs)
{
	lib.sizes <- IRanges::as.list(BSgenome::gdapply(cs, function(x) length(x[["+"]]) + length(x[["-"]]) ))
	lib.sizes <- lapply(lib.sizes, IRanges::as.list)
	sapply(lib.sizes, function(x) sum(unlist(x)))
})
