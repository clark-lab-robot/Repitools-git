setGeneric("enrichmentCalc", function(x, ...) standardGeneric("enrichmentCalc"))

setMethod("enrichmentCalc", "GenomeDataList",
    function(x, ...)
{
    enrichmentCalc(.GDL2GRL(x), ...)
})

setMethod("enrichmentCalc", "GRangesList",
    function(x, organism, seqLen=NULL, verbose=TRUE)
{
    xNames <- if (is.null(names(x))) 1:length(x) else names(x)
    ans <- lapply(1:length(x), function(i)
        {
            if (verbose) cat(xNames[i])
            enrichmentCalc(x[[i]], organism, seqLen, verbose)
        })
    if (verbose) cat("\n")
    ans
})

setMethod("enrichmentCalc", "GRanges",
    function(x, organism, seqLen=NULL, verbose=TRUE)
{
    if (!all(levels(seqnames(x)) %in% seqnames(organism)))
        stop("Chromosome name mismatch bewteen 'x' and 'organism'")
    seqlengths(x) <- seqlengths(organism)[levels(seqnames(x))]
    if (!is.null(seqLen)) x <- resize(x, seqLen)
    covTable <- colSums(table(coverage(x)))
    covTable <- data.frame(as.numeric(names(covTable)), covTable)
    colnames(covTable) <- c("coverage", "bases")
    if (verbose) cat("; ")
    covTable
})
