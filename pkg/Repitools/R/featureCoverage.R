setGeneric("featureCoverage", signature = c("x", "anno"), function(x, anno, ...)
                                           {standardGeneric("featureCoverage")})
setGeneric(".validate", signature = c("anno"), function(anno, up, down)
                                        {standardGeneric(".validate")})

setMethod(".validate", "GRanges", function(anno, up, down)
{
    str <- strand(anno)

    if(-up > down)
    {
	if('*' %in% str)
		stop("Boundaries closer to starts of chromosomes are closer to
                      the ends of the chromosomes, than the endmost boundaries.\n")
	else
		stop("Upstream boundaries are more downstream than downstream
                      boundaries.\n")
    }
	
    if(any(str == '*') && any(str %in% c('+', '-')))
	stop("Annotation contains mixed feature types.")

    # For unstranded features.
    if(any(str %in% '*') && up != down)
	stop("Different upstream and downstream window distances don't make
              sense for unstranded features.\n")
    validated <- NULL; rm(validated); # Dummy to trick R CMD check
    validated <<- TRUE
})

setMethod("featureCoverage", c("GRanges", "GRanges"), function(x, anno, up = 50,
           down = 100, dist = c("percent", "base"), freq = 10, s.width = 1000,
           verbose = TRUE)
{
    if(!"validated" %in% ls())
    {
	.validate(anno, up, down)
	dist <- match.arg(dist)
    }

    if(!"cvgSamps" %in% ls())
    {
	str <- strand(anno)
	st <- start(anno)
	en <- end(anno)	
	wd <- width(anno) 
	pos <- str == '+'

	if(verbose == TRUE) cat("Calculating sampling positions.\n")	
	if(all(str == '*'))
	{
	    refPoints <- round((st + en) / 2)
	} else {
	    refPoints <- ifelse(pos, st, en)
	}

	if(dist == "percent")
	{
	    starts = as.numeric(ifelse(pos, refPoints - wd * up/100,
                                            refPoints - wd * down/100))
	    ends = as.numeric(ifelse(pos, refPoints + wd * down/100,
                                          refPoints + wd * up/100))
	} else {
	    starts = as.numeric(ifelse(pos, refPoints - up,
                                            refPoints - down))
	    ends = as.numeric(ifelse(pos, refPoints + down,
                                          refPoints + up))
	}

	covWinds <- GRanges(seqnames = seqnames(anno),
                            IRanges(start = starts, end = ends),
                            strand = str)

	posns <- seq(-up, down, freq)
	if(dist == "percent")
	    posLabels <<- paste(posns, '%')
	else
	    posLabels <<- posns
	nPos <- length(posns)

	# Make ranges for each sample point.
	cvgSamps <<- rep(covWinds, each = nPos)
	gapSize <- ifelse(dist == "percent", width(cvgSamps) / (nPos - 1),
                                             rep(freq, length(cvgSamps)))
	ranges(cvgSamps) <- IRanges(start = as.numeric(
                              ifelse(strand(cvgSamps) %in% c('+', '*'), 
                                     start(cvgSamps) + 0:(nPos - 1) * gapSize,
                                     end(cvgSamps) - 0:(nPos - 1) * gapSize)
                                                      ),
                                     width = 1)

	# Find upper bound of how far a sampling position could be
        # outside of a chromosome.
	maxOut <<- ifelse(dist == "percent", max(width(anno)) * 
                                                 max(abs(up), abs(down)) / 100,
                                             max(abs(up), abs(down)))
	cvgSamps <- shift(cvgSamps, maxOut)

	# Find order to get back from RangesList order to original order.
	chrOrd <<- order(as.character(seqnames(anno)))
	anno <<- anno[chrOrd]
	oldOrd <<- order(chrOrd)
    }

    # Qualitatively near identical to running mean smoothing.
    if(verbose == TRUE) cat("Extending all reads to smoothing width.\n")
    seqlengths(x) <- rep(NA, length(seqlengths(x)))
    x <- resize(x, s.width)

    # Infer chromosome end positions from feature annotations.
    # This means the user doesn't have to have a BSgenome object of their samples'
    # genome installed on their computer.
    maxAnno <- sapply(split(ranges(anno), seqnames(anno)), function(y) max(end(y)))
    maxReads <- sapply(split(ranges(x), seqnames(x)), function(y) max(end(y)))
    maxReads <- maxReads[names(maxAnno)]
    chrLens <- mapply(max, maxAnno, maxReads)
    gc()

    # Remove reads on chromosomes not in annotation.
    # Shift all reads and sampling positions right by the out-of-bounds upper bound,
    # and then adjust chr lengths accordingly.
    x <- x[seqnames(x) %in% names(chrLens)]
    x <- shift(x, maxOut)
    seqlengths(x) <- chrLens + 2 * maxOut
    gc()

    # Get coverage.
    if(verbose == TRUE) cat("Calculating coverage.\n")
    # Scale all coverages for total reads.
    cvg <- coverage(x)
    sampChr <- as(cvgSamps, "RangesList")
    gc()

    # Do sampling, per chromosome.
    if(verbose == TRUE) cat("Sampling coverage.\n")
    cvgMat <- lapply(names(sampChr), function(y)
	            {
			inds <- sampChr[[y]]
			chrMat <- cvg[[y]][inds, drop = TRUE]
			matrix(chrMat, ncol = length(posLabels), byrow = TRUE)		   	
		    })
    cvgMat <- do.call(rbind, cvgMat)
    	
    colnames(cvgMat) <- posLabels
    if("name" %in% names(elementMetadata(anno))) 
	rownames(cvgMat) <- elementMetadata(anno)[, "name"]
    else
	rownames(cvgMat) <- 1:nrow(cvgMat)
    cvgMat <- cvgMat[oldOrd, ]

    rm(cvg)
    if(length(grep("featureCoverage", sys.calls())) <= 2)
	rm(cvgSamps, posLabels, maxOut, chrOrd, anno, oldOrd, validated, pos = ".GlobalEnv")
    gc()
    new("CoverageList", marks = "Undefined", cvgs = list(cvgMat), up = up, down = down,
                        dist = dist, freq = freq, s.width = s.width)
})

setMethod("featureCoverage", c("GRangesList", "GRanges"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
	dist <- match.arg(dist)
	.validate(anno, up, down)
	
	cvgs <- lapply(as.list(x), function(y) 
	           {
		   	featureCoverage(y, anno, up, down, dist, freq, s.width, verbose)
		   })
	rm(cvgSamps, posLabels, maxOut, chrOrd, anno, oldOrd, validated, pos = ".GlobalEnv")
	if(!is.null(names(x)))
		marks <- names(x)
	else
		marks <- unname(sapply(cvgs, names))
	new("CoverageList", marks = marks, cvgs = unname(sapply(cvgs, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod("featureCoverage", c("character", "GRanges"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
	dist <- match.arg(dist)
	.validate(anno, up, down)

	cvgs <- lapply(x, function(y)
	           {
		   	readsGR <- BAM2GRanges(y)
		   	featureCoverage(readsGR, anno, up, down, dist, freq,
                                        s.width, verbose)
		   })
	rm(cvgSamps, posLabels, maxOut, chrOrd, anno, oldOrd, validated, pos = ".GlobalEnv")
	if(!is.null(names(x)))
		marks <- x
	else
		marks <- unname(sapply(cvgs, names))
	new("CoverageList", marks = marks, cvgs = unname(sapply(cvgs, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod("featureCoverage", c("ANY", "data.frame"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
	colMissing <- setdiff(c("chr", "start", "end"), colnames(anno))
	if(length(colMissing > 0))
		stop("Columns ", paste(colMissing, collapse = ", "),
                     " of annotation are not present.")

	annoGR <- GRanges(anno$chr,
                          IRanges(anno$start, anno$end),
                          if("strand" %in% colnames(anno)) anno$strand else '*',
			  name = if("name" %in% colnames(anno)) anno$name)
	featureCoverage(x, annoGR, up, down, dist, freq, s.width, verbose)
})

setMethod("featureCoverage", c("GenomeDataList", "ANY"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
	dist <- match.arg(dist)
	.validate(anno, up, down)
	x <- .GDL2GRL(x)

	featureCoverage(x, anno, up, down, dist, freq, s.width, verbose)
})
