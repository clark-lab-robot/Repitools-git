setGeneric("sequenceCalc", function(x, ...) standardGeneric("sequenceCalc"))

setMethod("sequenceCalc", "GRanges",
    function(x, organism, pattern, fixed=TRUE, positions=FALSE)
{
    chrs <- unique(as.character(seqnames(x)))
    names(chrs) <- chrs
    if (!all(chrs %in% seqnames(organism))) stop("Chromosome name mismatch bewteen x and organism")
    hits <- as(RangesList(lapply(chrs, function(x) IRanges(start(matchPattern(pattern, organism[[x]], fixed=TRUE)), width=1))), "GRanges")
    if (!positions) return(countOverlaps(x, hits))
    scores <- vector(mode='list', length=length(x))
    temp <- findOverlaps(x, hits)@matchMatrix
    temp[,2] <- start(hits)[temp[,2]]-start(x)[temp[,1]]
    scores[unique(temp[,1])] <- split(temp[,2], temp[,1])
    scores
})

setMethod("sequenceCalc", "data.frame",
    function(x, window=500, organism, positions=FALSE, ...)
{
    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width=1), seqlengths=seqlengths(organism)[unique(x$chr)])
    x <- resize(x, window, fix="center")
    if (positions) {
        pos <- sequenceCalc(x, organism, positions=TRUE, ...)
        return(lapply(pos, function(x) if (!is.null(x)) return(x-window/2)))
    } else sequenceCalc(x, organism, ...)
})
