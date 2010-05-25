ChromaBlocks <- function(rs.ip, rs.input, organism, chrs, ipWidth=100, inputWidth=500, preset=NULL, blockWidth=NULL, cutoff=NULL, minBlocks=NULL, extend=NULL, verbose=TRUE, seqLen=NULL) {
	if (preset=="small") {
		blockWidth=10
		cutoff=1.5
		minBlocks=5
	} else if (preset=="large") {
		blockWidth=50
		cutoff=0.5
		minBlocks=25
		extend=0.1
	} else stopifnot(!is.null(blockWidth), !is.null(cutoff), !is.null(minBlocks))
	if (verbose) cat("Creating bins\n")
    IPbins <- genomeBlocks(organism, chrs, ipWidth)
    InputBins <- genomeBlocks(organism, chrs, inputWidth, ipWidth)
	if (verbose) cat("Counting IP lanes\n")
	ipCounts <- annotationBlocksCounts(rs.ip, IPbins, seqLen=seqLen, verbose=verbose)
	#pool & normalise IP lanes & turn into RPKM - reads per kb (ipWidth/1000) per million (/lanecounts*1000000)
	ipCounts <- rowSums(ipCounts)/sum(laneCounts(rs.ip))*1000000/(ipWidth/1000)

	if (verbose) cat("Counting Input lanes\n")
	inputCounts <- annotationBlocksCounts(rs.input, InputBins, seqLen=seqLen, verbose=verbose)
	#pool & normalise Input lanes
	inputCounts <- rowSums(inputCounts)/sum(laneCounts(rs.input))*1000000/(inputWidth/1000)

	RPMK <- ipCounts-inputCounts
	#scale RPMK to have mean==0
	RPMK <- RPMK-mean(RPMK)
	#okay find regionz
	if (verbose) cat("Smoothing and calling regions\n")
	RPMK.split <- GenomeData(split(data.frame(position=IPbins$position, RPMK=RPMK, stringsAsFactors=FALSE), IPbins$chr))
	regions <- GenomeData(vector(mode='list', length=length(chrs)))
	names(regions) <- names(RPMK.split)
	for (chr in 1:length(RPMK.split)) {
		if (verbose) cat(names(RPMK.split)[chr], " ", sep="")
		RPMK.split[[chr]]$mean <- Ringo::sliding.meansd(RPMK.split[[chr]]$position, RPMK.split[[chr]]$RPMK, ipWidth*blockWidth/2)[,"mean"]
		RPMK.split[[chr]]$cutoff <- Ringo::sliding.meansd(RPMK.split[[chr]]$position, as.integer(RPMK.split[[chr]]$RPMK>=cutoff), ipWidth*blockWidth/2)[,"mean"]>=minBlocks/blockWidth
		temp <- RPMK.split[[chr]]$mean>cutoff&RPMK.split[[chr]]$cutoff
		#clean the ends of the chromosomes
		temp[c(1:(blockWidth/2), (length(temp)-(blockWidth/2)):length(temp))] <- FALSE
		#Expand the regions out to the ends
		toExtend <- which(temp)
		for (i in -(blockWidth/2):(blockWidth/2)) temp[toExtend+i] <- TRUE
		tempRle <- Rle(temp)
		tempStarts <- start(tempRle)[runValue(tempRle)==TRUE]
		tempEnds <- end(tempRle)[runValue(tempRle)==TRUE]
		regions[[chr]] <- IRanges(start=RPMK.split[[chr]]$position[tempStarts], end=RPMK.split[[chr]]$position[tempEnds])
		}
	cat("\n")
	GenomeDataList(list(data=RPMK.split, regions=regions))
}