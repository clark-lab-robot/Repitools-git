setGeneric("mergeReplicates", signature = "reads", function(reads, types, ...)
                                         {standardGeneric("mergeReplicates")})

setMethod("mergeReplicates", "GRangesList", function(reads, types, verbose = TRUE)
{
	if(is.null(types) == TRUE)
		stop("Mandatory argument 'types' not provided.\n")
	if(length(types) != length(reads))
		stop("'types' and 'reads' lengths differ.\n")

	if(verbose == TRUE) cat("Unlisting GRangesList.\n")
	readsGR <- unlist(reads, use.names = FALSE)
	rdTypes <- Rle(types, elementLengths(reads))
	if(verbose == TRUE) cat("Splitting by types.\n")
	reads <- split(readsGR, rdTypes)
	metadata(reads) <- list(names(reads))
  	gc()
	if(verbose == TRUE) cat("Pooled GRangesList created.\n")
	reads
})

setMethod("mergeReplicates", "GenomeDataList", function(reads, types, verbose = TRUE)
{
	if(is.null(types) == TRUE)
		stop("Mandatory argument 'types' not provided.\n")
	if(length(types) != length(reads))
		stop("'types' and 'reads' lengths differ.\n")

	reads <- .GDL2GRL(reads)
	mergeReplicates(reads, types, verbose)
})
