setGeneric("cpgDensityCalc", function(x, ...) standardGeneric("cpgDensityCalc"))

setMethod("cpgDensityCalc", "GenomeDataList",
function(x, ...)
{
    sequenceCalc(.GDL2GRL(x), ...)
})

setMethod("cpgDensityCalc", "GRangesList",
function(x, ...)
{
    sequenceCalc(.GDL2GRL(x), ...)
})

setMethod("cpgDensityCalc", "GRanges",
    function(x, seq.len, ...)
{
    x <- resize(x, seq.len)
    positionsDF <- data.frame(chr = as.character(seqnames(x)), position = round((start(x) + end(x)) / 2))
    return(cpgDensityCalc(positionsDF, ...))
})

setMethod("cpgDensityCalc", "data.frame",
    function(x, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=TRUE)
{
    wFunction <- match.arg(wFunction)
    if(wFunction == "none") {
        cpgDensity <- sequenceCalc(x, window, organism, FALSE, DNAString("CG"))
    } else {
        CGfinds <- sequenceCalc(x, window, organism, TRUE, DNAString("CG"))
        distances <- lapply(CGfinds, function(u) if(!is.null(u)) abs(u) else u)
        if(wFunction == "linear") {
            cpgDensity <- sapply(distances, function(d) sum(1-(d/(window/2))))
        } else if(wFunction == "log") {
            cpgDensity <- sapply(distances, function(d) sum(log2(2-(d/(window/2)))))
        } else {
            cpgDensity <- sapply(distances, function(d) sum(exp(-5*d/(window/2))))	
        }
        rm(CGfinds)
    }
    gc()
    return(cpgDensity)
})
