setGeneric("gcContentCalc", function(x, ...){standardGeneric("gcContentCalc")})

setMethod("gcContentCalc", "GRanges", function(x, organism) {
    strand(x) <- "+"
    temp <- getSeq(organism, x, as.character=FALSE)
    tempAlphabet <- alphabetFrequency(temp, as.prob=TRUE)
    (tempAlphabet[,"C"]+tempAlphabet[,"G"])
})

setMethod("gcContentCalc", "data.frame", function(x, window=500, organism) {
    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width=1), seqlengths=seqlengths(organism)[unique(x$chr)])
    x <- resize(x, window, fix="center")
    gcContentCalc(x, organism)
})
