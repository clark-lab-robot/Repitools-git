laneCounts <- function(cs) {
	if ("GenomeData" %in% class(cs)) cs <- GenomeDataList(list(cs))
	stopifnot("GenomeDataList" %in% class(cs))
	#require(IRanges)
	#require(BSgenome)
        # modification by Mark R -- Jan 20 2010
	lib.sizes <- IRanges::as.list(BSgenome::gdapply(cs, function(x) length(x[["+"]]) + length(x[["-"]]) ))
	lib.sizes <- lapply(lib.sizes, IRanges::as.list)
	sapply(lib.sizes, function(x) sum(unlist(x)))
        #chrSizes <- gdapply(cs, FUN=function(v) length(v[["+"]])+length(v[["-"]]))
        #gdapply(cs, FUN=function(u) sum( sapply(u, FUN=function(v) length(v[["+"]])+length(v[["-"]])) ) )
}
