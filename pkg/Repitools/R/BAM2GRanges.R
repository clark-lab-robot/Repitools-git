setGeneric("BAM2GRanges", signature = "path", function(path, ...)
                                {standardGeneric("BAM2GRanges")})
setGeneric("BAM2GRangesList", signature = "paths", function(paths, ...)
                                  {standardGeneric("BAM2GRangesList")})

setMethod("BAM2GRanges", "character",
    function(path, what = c("rname", "strand", "pos", "qwidth"),
             flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
	     verbose = TRUE)
{
    require(Rsamtools)
    if(length(path) > 1)
	stop("This method is only for one BAM file. See ?BAM2GRangesList.")

    if(verbose == TRUE)
	cat("Reading BAM file ", basename(path), ".\n", sep = '')
    filters <- ScanBamParam(what = what, flag = flag)
    bam <- scanBam(path, param = filters)[[1]]
    anno <- !names(bam) %in% c("rname", "strand", "pos", "qwidth")
    if(verbose == TRUE) cat("Creating GRanges from BAM datalist.\n")
    gr <- GRanges(bam$rname, IRanges(bam$pos, width = bam$qwidth),
                  bam$strand, bam[anno])
})

setMethod("BAM2GRangesList", "character",
    function(paths, what = c("rname", "strand", "pos", "qwidth"),
             flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
             verbose = TRUE)
{
    require(Rsamtools)
    GRangesList(lapply(paths, function(x) BAM2GRanges(x, what, flag)))
})
