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

    if(!"cvg.samps" %in% ls())
    {
		str <- strand(anno)
		st <- start(anno)
		en <- end(anno)	
		wd <- width(anno) 
		pos <- str == '+'

		if(verbose == TRUE) cat("Calculating sampling positions.\n")	
		if(all(str == '*'))
		{
	    	ref.points <- round((st + en) / 2)
		} else {
	    	ref.points <- ifelse(pos, st, en)
		}

		if(dist == "percent")
		{
	    	starts = as.numeric(ifelse(pos, ref.points - wd * up/100,
                                            ref.points - wd * down/100))
	    	ends = as.numeric(ifelse(pos, ref.points + wd * down/100,
                                          ref.points + wd * up/100))
		} else {
	    	starts = as.numeric(ifelse(pos, ref.points - up,
                                            ref.points - down))
	    	ends = as.numeric(ifelse(pos, ref.points + down,
                                          ref.points + up))
		}

	cov.winds <- GRanges(seqnames = seqnames(anno),
                         IRanges(start = starts, end = ends),
                                 strand = str)

	posns <- seq(-up, down, freq)
	if(dist == "percent")
	    pos.labels <<- paste(posns, '%')
	else
	    pos.labels <<- posns
	n.pos <- length(posns)

	# Make ranges for each sample point.
	cvg.samps <<- rep(cov.winds, each = n.pos)
	gap.size <- ifelse(dist == "percent", width(cvg.samps) / (n.pos - 1),
                                          rep(freq, length(cvg.samps)))
	ranges(cvg.samps) <- IRanges(start = as.numeric(
                              ifelse(strand(cvg.samps) %in% c('+', '*'), 
                                     start(cvg.samps) + 0:(n.pos - 1) * gap.size,
                                     end(cvg.samps) - 0:(n.pos - 1) * gap.size)
                                                      ),
                                     width = 1)

		# Find upper bound of how far a sampling position could be
        # outside of a chromosome.
		max.out <<- ifelse(dist == "percent", max(width(anno)) * 
                                             max(abs(up), abs(down)) / 100,
                                             max(abs(up), abs(down)))
		cvg.samps <- shift(cvg.samps, max.out)

		# Find order to get back from RangesList order to original order.
		chr.ord <<- order(as.character(seqnames(anno)))
		anno.chr <<- anno[chr.ord]
		old.ord <<- order(chr.ord)
    }

    # Qualitatively near identical to running mean smoothing.
    if(verbose == TRUE) cat("Extending all reads to smoothing width.\n")
    seqlengths(x) <- rep(NA, length(seqlengths(x)))
    x <- resize(x, s.width)

    # Infer chromosome end positions from feature annotations.
    # This means the user doesn't have to have a BSgenome object of their samples'
    # genome installed on their computer.
    max.anno <- sapply(split(ranges(anno.chr), seqnames(anno.chr)), function(y) max(end(y)))
    max.reads <- sapply(split(ranges(x), seqnames(x)), function(y) max(end(y)))
    max.reads <- max.reads[names(max.anno)]
    chr.lens <- mapply(max, max.anno, max.reads)
    gc()

    # Remove reads on chromosomes not in annotation.
    # Shift all reads and sampling positions right by the out-of-bounds upper bound,
    # and then adjust chr lengths accordingly.
    x <- x[seqnames(x) %in% names(chr.lens)]
    x <- shift(x, max.out)
    seqlengths(x) <- chr.lens + 2 * max.out
    gc()

    # Get coverage.
    if(verbose == TRUE) cat("Calculating coverage.\n")
    # Scale all coverages for total reads.
    cvg <- coverage(x) / length(x)
    samp.chr <- as(cvg.samps, "RangesList")
    gc()

    # Do sampling, per chromosome.
    if(verbose == TRUE) cat("Sampling coverage.\n")
    cvg.mat <- lapply(names(samp.chr), function(y)
	           {
			     inds <- samp.chr[[y]]
                 chr.mat <- cvg[[y]][inds, drop = TRUE]
			     matrix(chr.mat, ncol = length(pos.labels), byrow = TRUE)		   	
		    })
    cvg.mat <- do.call(rbind, cvg.mat)
    	
    colnames(cvg.mat) <- pos.labels
    if("name" %in% names(elementMetadata(anno.chr))) 
	rownames(cvg.mat) <- elementMetadata(anno.chr)[, "name"]
    else
	rownames(cvg.mat) <- 1:nrow(cvg.mat)
    cvg.mat <- cvg.mat[old.ord, ]

    rm(cvg)
    if(length(grep("featureCoverage", sys.calls())) <= 2)
	rm(cvg.samps, pos.labels, max.out, chr.ord, anno.chr, old.ord, validated,
        pos = ".GlobalEnv")
    gc()

    new("CoverageList", marks = "Undefined", cvgs = list(cvg.mat), anno = anno,
         up = up, down = down, dist = dist, freq = freq, s.width = s.width)
})

setMethod("featureCoverage", c("GRangesList", "GRanges"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
    dist <- match.arg(dist)
    .validate(anno, up, down)
    
    if(length(s.width) == 1) s.width <- rep(s.width, length(x))
    cvgs <- mapply(function(y, z) 
	           {
		   	featureCoverage(y, anno, up, down, dist, freq, z, verbose)
		   }, as.list(x), s.width, SIMPLIFY = FALSE)
    rm(cvg.samps, pos.labels, max.out, chr.ord, anno.chr, old.ord, validated,
       pos = ".GlobalEnv")
    if(!is.null(names(x)))
	marks <- names(x)
    else
	marks <- unname(sapply(cvgs, names))
    new("CoverageList", marks = marks, anno = anno, cvgs = unname(sapply(cvgs, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod("featureCoverage", c("character", "GRanges"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
    dist <- match.arg(dist)
    .validate(anno, up, down)

    if(length(s.width) == 1) s.width <- rep(s.width, length(x))
    cvgs <- mapply(function(y, z)
	           {
		   	readsGR <- BAM2GRanges(y)
		   	featureCoverage(readsGR, anno, up, down, dist, freq,
                                        s.width, verbose)
		   }, x, s.width, SIMPLIFY = FALSE)
    rm(cvg.samps, pos.labels, max.out, chr.ord, anno.chr, old.ord, validated,
       pos = ".GlobalEnv")
    if(!is.null(names(x)))
	marks <- x
    else
	marks <- unname(sapply(cvgs, names))
    new("CoverageList", marks = marks, anno = anno, cvgs = unname(sapply(cvgs, tables)),
                            up = up, down = down, dist = dist, freq = freq,
                            s.width = s.width)
})

setMethod("featureCoverage", c("ANY", "data.frame"),
    function(x, anno, up = 50, down = 100, dist = c("percent", "base"), freq = 10,
             s.width = 1000, verbose = TRUE)
{
    col.missing <- setdiff(c("chr", "start", "end"), colnames(anno))
    if(length(col.missing > 0))
	stop("Columns ", paste(col.missing, collapse = ", "),
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
