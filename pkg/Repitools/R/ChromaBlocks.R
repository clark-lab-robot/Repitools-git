ChromaBlocks <- function(rs.ip, rs.input, organism, chrs, ipWidth=100, inputWidth=500, preset=NULL, blockWidth=NULL, minBlocks=NULL, extend=NULL, cutoff=NULL, FDR=0.01, nPermutations=5, nCutoffs=20, cutoffQuantile=0.98, verbose=TRUE, seqLen=NULL) {

    mergeOverlaps <- function (query, subject) {
        ov <- matchMatrix(findOverlaps(query, subject, select="all"))
        query.ov <- unique(ov[,1])
        subjectStarts <- tapply(ov[,2], ov[,1], function(x) min(start(subject[x])))
        subjectEnds <- tapply(ov[,2], ov[,1], function(x) max(end(subject[x])))
        start(query)[query.ov] <- ifelse(start(query[query.ov])>subjectStarts, subjectStarts, start(query[query.ov]))
        end(query)[query.ov] <- ifelse(end(query[query.ov])<subjectEnds, subjectEnds, end(query[query.ov]))
        query
    }
    
    callRegions <- function(bins, RPKM=NULL, cutoffs) {
        callCutoff <- function(dat, cutoff) {
            dat$cutoff <- Ringo::sliding.meansd(dat$position, as.integer(dat$RPKM>=cutoff), ipWidth*blockWidth/2)[,"mean"]>=minBlocks/blockWidth
            temp <- dat$mean>cutoff&dat$cutoff
            #clean the ends of the chromosomes
            temp[c(1:(blockWidth/2), (length(temp)-(blockWidth/2)):length(temp))] <- FALSE
            #Expand the regions out to the ends
            toExtend <- which(temp)
            for (i in -(blockWidth/2):(blockWidth/2)) temp[toExtend+i] <- TRUE
            tempRle <- Rle(temp)
            tempStarts <- start(tempRle)[runValue(tempRle)==TRUE]
            tempEnds <- end(tempRle)[runValue(tempRle)==TRUE]
            if (!is.null(extend)){
                extendRle <- Rle(dat$extend)
                extendRanges <- IRanges(dat$position[start(extendRle)[runValue(extendRle)==TRUE]], dat$position[end(extendRle)[runValue(extendRle)==TRUE]])
                mergeOverlaps(IRanges(dat$position[tempStarts], dat$position[tempEnds]), extendRanges)
            } else IRanges(dat$position[tempStarts], dat$position[tempEnds])
        }
        
        callChr <- function(dat) {
            if (verbose) cat(".")
            dat$mean <- Ringo::sliding.meansd(dat$position, dat$RPKM, ipWidth*blockWidth/2)[,"mean"]
            if (!is.null(extend)) dat$extend <- dat$mean<extend 
            if (length(cutoffs)>1) sapply(cutoffs, function(x) callCutoff(dat, x)) else callCutoff(dat, cutoffs)
        }
        
        RPKM.split <- GenomeData(split(data.frame(position=(start(bins)+end(bins))/2, RPKM=if (is.null(RPKM)) bins$RPKM else RPKM), space(bins)))
        regions <- lapply(RPKM.split, callChr)
        if (verbose) cat("\n")
        if (length(cutoffs)>1) rowSums(sapply(regions, function(x) sapply(x,length))) else GenomeData(regions)
    }
    
    if (preset=="small") {
        blockWidth=10
        minBlocks=5
    } else if (preset=="large") {
        blockWidth=50
        minBlocks=25
        extend=0.1
    } else stopifnot(!is.null(blockWidth), !is.null(minBlocks))
    if (verbose) cat("Creating bins\n")
    IPbins <- genomeBlocks(organism, chrs, ipWidth)
    InputBins <- genomeBlocks(organism, chrs, inputWidth, ipWidth)
    if (verbose) cat("Counting IP lanes: ")
    ipCounts <- annotationBlocksCounts(rs.ip, IPbins, seqLen=seqLen, verbose=verbose)
    #pool & normalise IP lanes & turn into RPKM - reads per kb (ipWidth/1000) per million (/lanecounts*1000000)
    ipCounts <- rowSums(ipCounts)/sum(elementLengths(rs.ip))*1000000/(ipWidth/1000)
    
    if (verbose) cat("Counting Input lanes: ")
    inputCounts <- annotationBlocksCounts(rs.input, InputBins, seqLen=seqLen, verbose=verbose)
    #pool & normalise Input lanes
    inputCounts <- rowSums(inputCounts)/sum(elementLengths(rs.input))*1000000/(inputWidth/1000)
    
    IPbins$RPKM <- ipCounts-inputCounts
    rm(ipCounts, inputCounts)
    #scale RPMK to have mean==0
    IPbins$RPKM <- IPbins$RPKM-mean(IPbins$RPKM)
    #okay find regions
    if (is.null(cutoff)) {
        cutoffs <- seq(0.1, quantile(IPbins$RPKM, cutoffQuantile), length.out=nCutoffs)
        negRegions <- sapply(1:nPermutations, function(u) {
            if (verbose) cat("Permutation",u)
            callRegions(IPbins, IPbins$RPKM[sample(nrow(IPbins))], cutoffs=cutoffs)
        })
        if (verbose) cat("Testing positive regions")
        posRegions <- callRegions(IPbins, cutoffs=cutoffs)
        FDRTable <- cbind(cutoffs=cutoffs, pos=posRegions, negRegions, FDR=rowMeans(negRegions)/posRegions)
        cutoff <- cutoffs[match(TRUE, FDRTable[,"FDR"]<FDR)]
        if (is.na(cutoff)) {
            warning("No cutoff below FDR", FDR, "found! Analysis halted! Try increasing cutoffQuantile or FDR")
            return(data=IPbins, FDRTable=FDRTable)
        }
        if (verbose) cat("Using cutoff of",cutoff,"for a FDR of",FDR,"\n")
    } else {
        if (verbose) cat("Using supplied cutoff of",cutoff,"\n")
        FDRTable <- NULL
    }
    list(data=IPbins, regions=callRegions(IPbins, cutoffs=cutoff), FDRTable=FDRTable, cutoff=cutoff)
}
