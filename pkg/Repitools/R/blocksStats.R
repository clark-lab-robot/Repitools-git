.blocksStats <- function(diffs, coordinatesTable, design, upStream, downStream, verbose, robust=FALSE, minNRobust, adjustMethod, useAsRegions, annot)
{       
  means <- tstats <- matrix(NA, nr=nrow(coordinatesTable), nc=ncol(diffs), dimnames=list(NULL,colnames(design)))
  df <- rep(0,nrow(coordinatesTable))

  for(i in 1:nrow(coordinatesTable)) {
    vind <- annot[[1]][[i]]
    vind <- unique(vind[!is.na(vind)])
    if( length(vind) < 2 )
      next
        for(j in 1:ncol(diffs)) {
          if( robust & length(vind) >= minNRobust ) {
            require(MASS)
        rf <- summary(rlm(diffs[vind,j]~1))
        tstats[i,j] <- rf$coef[3]
            means[i,j] <- rf$coef[1]
        df[i] <- rf$df[2]
          } else {
        tt <- t.test(diffs[vind,j])
        tstats[i,j] <- tt$statistic
            means[i,j] <- tt$estimate
        df[i] <- tt$parameter
          }
        }
  }

  pvals <- 2*pt(-abs(tstats),df)
  
  # adjust p-values for multiple testing
  adjpvals <- pvals
  for(i in 1:ncol(pvals))
    adjpvals[,i] <- p.adjust(pvals[,i],method=adjustMethod)
  
  xDf <- data.frame(df=df, meandiff=means, tstats=tstats, pvals=pvals, adjpvals=adjpvals)
  if(useAsRegions == TRUE)
        {
                colnames(xDf)[1] <- "df"        
        }
        else
        {
                colnames(xDf)[1] <- paste("df",upStream,downStream,sep=".")
        }
        
  if( ncol(xDf)==5 )
    colnames(xDf)[2:5] <- paste(c("meandiff","tstats","pvals","adjpvals"), gsub(".[1-9]$","",colnames(xDf)[2:5]), sep=".")
   

  cbind(coordinatesTable,xDf)
}

setMethodS3("blocksStats", "AffymetrixCelSet", function(cs, coordinatesTable, annot=NULL, probePositions = NULL, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust=TRUE, useAsRegions=FALSE, ...)
{
	require(aroma.affymetrix)

  # cs - AffymetrixCelSet to read probe-level data from
  # coordinatesTable - data frame giving genome coordinates or regions of interest coordinates
  # ind - (optional) 
  # dmP - data matrix of probes
  
  if( nrow(design) != nbrOfArrays(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
  
  
	w <- which( rowSums(design != 0) > 0 )
	cs <- extract(cs, w, verbose=verbose)
	
	if(is.null(annot))
	{
		probePositions <- getProbePositionsDf( getCdf(cs), verbose=verbose )
		if(useAsRegions == TRUE)
		{
			if(!all(c("chr", "name", "start", "end")  %in% colnames(coordinatesTable)))
				stop("Incorrect column headings for coordinatesTable. Check documentation for details.")  
			annot <- annotationBlocksLookup(probePositions, coordinatesTable, verbose=verbose)
		}
		else
		{
			if(!all(c("chr", "name", "start", "end", "strand")  %in% colnames(coordinatesTable)))
				stop("Incorrect column headings for coordinatesTable. Check documentation for details.")  
			pos <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
			
			genePositions <- data.frame(chr=coordinatesTable$chr, position=pos, 
			strand=coordinatesTable$strand, row.names=coordinatesTable$name,
			stringsAsFactors=FALSE)
																											
			# run lookup twice.  first to get a list of smaller list of probes to use
			annot <- annotationLookup(probePositions, genePositions, upStream, downStream, verbose=verbose)
			pb <- unique(unlist(annot$indexes, use.names=FALSE))
			probePositions <- probePositions[pb,]
			annot <- annotationLookup(probePositions, genePositions, upStream, downStream, verbose=verbose)
		}
	}

	dmP <- extractMatrix(cs, cells=probePositions$index, verbose=verbose)
	if(log2adjust == TRUE)
	{
		diffs <- log2(dmP) %*% design[w,]
	}
	else
	{
		diffs <- dmP %*% design[w,]
	}
	
	return(.blocksStats(diffs, coordinatesTable, design, upStream, downStream, verbose, robust, minNRobust, adjustMethod, useAsRegions, annot))
  
})

setMethodS3("blocksStats", "GenomeDataList", function(cs, coordinatesTable, design, upStream=0, downStream=2000, verbose=TRUE, useAsRegions=FALSE, seqLen=NULL, libSize="lane", Acutoff=NULL, ...) {
	if(libSize == "ref" && is.null(Acutoff))
		stop("Must give value of Acutoff if using \"ref\" normalisation.\n")
	require(edgeR)
	if(!all(c("chr", "name", "start", "end")  %in% colnames(coordinatesTable)))
		stop("Incorrect column headings for coordinatesTable. Check documentation for details.")
	if (verbose) cat("Generating table of counts\n")
	if (useAsRegions) dm <- annotationBlocksCounts(cs, coordinatesTable, seqLen, verbose) else {
		coordinatesTable$position <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
		dm <- annotationCounts(cs, coordinatesTable, upStream, downStream, seqLen, verbose)
	}
	if (libSize == "lane")
		lib.sizes <- laneCounts(cs)
	if(libSize == "inRegions")
		lib.sizes <- colSums(dm)
	if(libSize == "ref")
		lib.sizes <- colSums(dm) * calcNormFactors(dm, Acutoff=Acutoff)
	
    dmRes <- cbind(coordinatesTable, dm)
    for (i in 1:ncol(design)) {
		if (verbose) cat("Processing column",i,"of design matrix\n")
		stopifnot(sum(design[,i]==1)>0, sum(design[,i]==-1)>0, all(design[,i] %in% c(-1,0,1)))
		thisCol <- design[,i]!=0
		d <- DGEList(counts=dm[,thisCol], group=as.character(design[thisCol,i]), lib.size=lib.sizes[thisCol])
		d.disp <- estimateCommonDisp(d)
		tmp <- d.disp$pseudo.alt
		colnames(tmp) <- paste(colnames(tmp), "pseudo",sep="_")
		dmRes <- merge(dmRes, tmp, all.x=TRUE, by.x="name", by.y="row.names")
		deD <- exactTest(d.disp, pair = c("-1","1"))
		de <- topTags(deD, n=nrow(deD$table))@.Data[[1]]
		colnames(de) <- paste(colnames(de), colnames(design)[i], sep="_")
		dmRes <- merge(dmRes, de, all.x = TRUE, by.x = "name", by.y = "row.names")
	}
	dmRes[match(coordinatesTable$name, dmRes$name),]
})


setMethodS3("blocksStats", "matrix", function(cs, ndf, coordinatesTable, annot=NULL, probePositions=NULL, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust=TRUE, useAsRegions=FALSE, ...)
{
	if( nrow(design) != ncol(cs) )
		stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix.")
		
	w <- which( rowSums(design != 0) > 0 )
	
	if(is.null(annot))
	{
		probePositions <- data.frame(chr = ndf$chr, position = ndf$position, index = ndf$index, strand = ndf$strand, stringsAsFactors=FALSE)
	
		if(useAsRegions == TRUE)
		{
			if(!all(c("chr", "name", "start", "end")  %in% colnames(coordinatesTable)))
				stop("Incorrect column headings for coordinatesTable. Check documentation for details.")
			annot <- annotationBlocksLookup(probePositions, coordinatesTable)
		}
		else
		{	
			if(!all(c("chr", "name", "start", "end", "strand")  %in% colnames(coordinatesTable)))
				stop("Incorrect column headings for coordinatesTable. Check documentation for details.")
			pos <- rep(NA,nrow(coordinatesTable))
			pos[coordinatesTable$strand=="+"] <- coordinatesTable$start[coordinatesTable$strand=="+"]
			pos[coordinatesTable$strand=="-"] <- coordinatesTable$end[coordinatesTable$strand=="-"]
	
			genePositions <- data.frame(chr=coordinatesTable$chr, position=pos, strand=coordinatesTable$strand, row.names=coordinatesTable$name, stringsAsFactors=FALSE)
	
			# run lookup twice.  first to get a list of smaller list of probes to use
			annot <- annotationLookup(probePositions, genePositions, upStream, downStream, verbose=verbose)
			pb <- unique(unlist(annot$indexes, use.names=FALSE))
			probePositions <- probePositions[pb,]
			annot <- annotationLookup(probePositions, genePositions, upStream, downStream, verbose=verbose)
		}
	}

	if(log2adjust == TRUE)
	{						
		diffs <- log2(cs[probePositions$index, ]) %*% design[w,]
	} else {
		diffs <- cs[probePositions$index, ] %*% design[w,]
	}

	return(.blocksStats(diffs, coordinatesTable, design, upStream, downStream, verbose, robust, minNRobust, adjustMethod, useAsRegions, annot))
  }
)
