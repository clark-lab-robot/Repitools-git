setGeneric("cpgDensityCalc", function(x, ...) standardGeneric("cpgDensityCalc"))

setMethod("cpgDensityCalc", "GenomeDataList",
function(x, ...)
{
    cpgDensityCalc(.GDL2GRL(x), ...)
})


setMethod("cpgDensityCalc", "GRangesList",
    function(x, window=500, wFunction=c("none","linear","exp","log"), organism, verbose=TRUE, seqLen=NULL)
{
    xNames <- if (is.null(names(x))) 1:length(x) else names(x)
    ans <- lapply(1:length(x), function(i)
        {
            if (verbose) cat(xNames[i])
            cpgDensityCalc(x[[i]], window, wFunction, organism, verbose, seqLen)
        })
    if (verbose) cat("\n")
    ans
})


setMethod("cpgDensityCalc", "GRanges",
    function(x, window=500, wFunction=c("none","linear","exp","log"), organism, verbose=TRUE, seqLen=NULL)
{
    wFunction <- match.arg(wFunction)
    if (!is.null(seqLen)) x <- resize(x, seqLen)
    x <- resize(x, window, fix="center")
    if(wFunction == "none") {
        cpgDensity <- sequenceCalc(x, organism, DNAString("CG"))
    } else {
        CGfinds <- sequenceCalc(x, organism, DNAString("CG"), positions=TRUE)
        CGfinds <- lapply(CGfinds, function(u) if(!is.null(u)) abs(u-window/2) else u)
        if(wFunction == "linear") {
            cpgDensity <- sapply(CGfinds, function(d) sum(1-(d/(window/2))))
        } else if(wFunction == "log") {
            cpgDensity <- sapply(CGfinds, function(d) sum(log2(2-(d/(window/2)))))
        } else {
            cpgDensity <- sapply(CGfinds, function(d) sum(exp(-5*d/(window/2))))	
        }
        rm(CGfinds)
    }    
    if (verbose) cat("; ")
    return(cpgDensity)
})

setMethod("cpgDensityCalc", "data.frame",
    function(x, organism, ...)
{
    if (is.null(x$position)) x$position <- ifelse(x$strand == '+', x$start, x$end)
    x <- GRanges(x$chr, IRanges(x$position, width=1), seqlengths=seqlengths(organism)[unique(x$chr)])
    cpgDensityCalc(x, organism=organism, ...)
})
