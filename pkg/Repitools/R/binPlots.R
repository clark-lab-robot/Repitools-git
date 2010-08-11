setMethodS3("binPlots", "default", function(rs, coordinatesTable, design=NULL, upStream=7500, downStream=2500, by=100, bw=300, libSize="lane", seqLen=NULL, verbose=FALSE, Acutoff=NULL, ...) {
	coordinatesTable$position <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
	rownames(coordinatesTable) <- coordinatesTable$name
	if(libSize == "ref" && is.null(Acutoff))
		stop("Must give value of Acutoff if using \"ref\" normalisation.\n")
	blockPos <- seq.int(-upStream, downStream, by)
	if (verbose) cat("made blockPos\n")
	annoBlocks <- data.frame(chr=rep(coordinatesTable$chr, each=length(blockPos)),
                                 start=rep(coordinatesTable$position-bw, each=length(blockPos)),
                                 end=rep(coordinatesTable$position+bw, each=length(blockPos)),
                                 strand=rep(coordinatesTable$strand, each=length(blockPos)))
	annoBlocks$start[annoBlocks$strand=="+"] <- annoBlocks$start[annoBlocks$strand=="+"] + blockPos
	annoBlocks$end[annoBlocks$strand=="+"] <- annoBlocks$end[annoBlocks$strand=="+"] + blockPos
	annoBlocks$start[annoBlocks$strand=="-"] <- annoBlocks$start[annoBlocks$strand=="-"] - blockPos
	annoBlocks$end[annoBlocks$strand=="-"] <- annoBlocks$end[annoBlocks$strand=="-"] - blockPos
	if (verbose) cat("made annoBlocks\n")
	if (!is.null(design)) {
		stopifnot(all(design %in% c(-1,0,1)), nrow(design)==length(rs))
		inUse <- !apply(design==0,1,all)
		design <- design[inUse, , drop = FALSE]
	} else inUse <- rep(TRUE, length(rs))
	annoCounts <- annotationBlocksCounts(rs[inUse], annoBlocks, seqLen, verbose)
	if (libSize == "ref") {
		if (verbose) cat("normalising to reference sample\n")
		annoCounts <- t(t(annoCounts) * calcNormFactors(annoCounts, Acutoff = Acutoff))
		
	} else { # libSize = "lane"
		if (verbose) cat("normalising to total library sizes\n")
		totalReads <- if(class(rs) == "GenomeDataList") laneCounts(rs) else elementLengths(rs[inUse])
		annoCounts <- t(t(annoCounts)/totalReads)*1000000	
	}
	if (verbose) cat("made annoCounts\n")
	if (!is.null(design)) {
		if (verbose) cat("applying design matrix\n")
		design <- apply(design, 2, function(x) {
					x[x==1] <- 1/sum(x==1)
					x[x==-1] <- -1/sum(x==-1)
					return(x)
				})
		annoCounts <- annoCounts %*% design 
	}
	annoTable <- matrix(1:nrow(annoCounts), byrow=TRUE, ncol=length(blockPos), nrow=nrow(coordinatesTable), dimnames=list(NULL, blockPos))
	if (verbose) cat("made annoTable\n")
	binPlots(annoCounts, annoTable, removeZeros=FALSE, useMean=TRUE, ...)
})

setMethodS3("binPlots", "AffymetrixCelSet", function(cs, probeMap=NULL, coordinatesTable=NULL, upStream=7500, downStream=2500, by=100, bw=300, log2adjust=TRUE, verbose=FALSE, ...) {			
	if (is.null(probeMap)) {
		if (is.null(coordinatesTable)) stop("Either probeMap or coordinatesTable must be supplied!")
		probePositions <- getProbePositionsDf( getCdf(cs), verbose=verbose )
		coordinatesTable$position <- ifelse(coordinatesTable$strand=="+", coordinatesTable$start, coordinatesTable$end)
		rownames(coordinatesTable) <- coordinatesTable$name
			
		# run lookup twice.  first to get a list of smaller list of probes to use
		annot <- annotationLookup(probePositions, coordinatesTable, upStream+bw, downStream+bw, verbose=verbose)
		pb <- unique(unlist(annot$indexes, use.names=FALSE))
		probePositions <- probePositions[pb,]
		annot <- annotationLookup(probePositions, coordinatesTable, upStream+bw, downStream+bw, verbose=verbose)
		lookupT <- makeWindowLookupTable(annot$indexes, annot$offsets,
				starts = seq(-upStream-bw, downStream-bw, by), ends = seq(-upStream+bw, downStream+bw, by))
	} else {
		if (verbose) cat("Using supplied probeMap\n")
		probePositions <- probeMap$probePositions
		lookupT <- probeMap$lookupT
	}
	
	dmM <- extractMatrix(cs, cells = probePositions$index, verbose = verbose)
	if (log2adjust) dmM <- log2(dmM)

	binPlots(dmM, lookupT, ...)
	invisible(list(lookupT=lookupT, probePositions=probePositions))
})


setMethodS3("binPlots", "matrix", function(dataMatrix, lookupTable, ordering, plotType=c("line","heatmap","terrain","boxplot"), nbins=10, cols=NULL, lwd=3, lty=1, sameScale=TRUE, symmScale=FALSE, verbose=FALSE, removeZeros=TRUE, useMean=FALSE, ...) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  plotType <- match.arg(plotType)
  if(!ncol(ordering) == ncol(dataMatrix)) {
    if (!ncol(ordering) == 1)
      stop("ordering must have either 1 column or the same number of columns as dataMatrix.")
      orderingIndex <- rep(1, ncol(dataMatrix))
  } else orderingIndex <- 1:ncol(dataMatrix)

  if( is.null(cols) ) {
    require(gplots)
	if(plotType=="line") {
	  cols <- colorpanel(nbins,"blue","green","red")
	} else {
	  cols <- colorpanel(64,"blue","white","red")
	}
  }
  
  label <- vector("list", ncol(ordering))
  .makeBins <- function(u) {
    if(class(u)=="numeric") {
      br <- quantile(u,p=(0:nbins)/nbins)
      list(breakpoints=br, intervals=cut(u,breaks=br))
	} else if(class(u)=="factor") {
	  nbins <- length(levels(u))
	  list(breakpoints=u, intervals=u)
	}
  }
  
  for(i in 1:ncol(ordering))
  {
	if(class(ordering[,i])=="numeric")
		label[[i]] <- " Order:"
	else if (class(ordering[,i])=="factor")
		label[[i]] <- " Factor:"
  }
  
  breaks <- apply(ordering, 2, .makeBins)
  if( plotType %in% c("line","heatmap","terrain")) {
    intensScores <- array(NA,dim=c(ncol(dataMatrix), ncol(lookupTable), nbins),
	                      dimnames=list(colnames(dataMatrix),colnames(lookupTable),NULL))
  } else {
    intensScores <- vector("list",ncol(dataMatrix))
	for(i in 1:length(intensScores))
		intensScores[[i]] <- vector("list", length(levels(breaks[[orderingIndex[i]]][["intervals"]])))
  }
  
  xval <- as.numeric(colnames(lookupTable))

  for(i in 1:ncol(dataMatrix)) {
	if (verbose) cat(colnames(ordering)[orderingIndex[i]],": ",sep="")
	cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	
	for(j in 1:length(cutLevels)){
		level <- cutLevels[j]
	    lookupTableSubset <- lookupTable[breaks[[orderingIndex[i]]][["intervals"]]==level, ]
	  if( plotType %in% c("line","heatmap","terrain")) {
	    intensScores[i,,j] <- .scoreIntensity(lookupTableSubset, intensities=dataMatrix[,i], minProbes=2, removeZeros=removeZeros, useMean=useMean)	
      } else {
		d <- .scoreIntensity(lookupTableSubset, intensities=dataMatrix[,i], minProbes=2, returnMatrix=TRUE, removeZeros=removeZeros, useMean=useMean)
		intensScores[[i]][[j]] <- boxplot(as.data.frame(d), plot=FALSE)
	  }
	}
  }

    if ( plotType %in% c("line","heatmap","terrain")) if (sameScale) {
      rng <- range(intensScores, na.rm=TRUE)
      if (symmScale) rng <- c(-max(abs(rng)),max(abs(rng)))
    }

  for(i in 1:ncol(dataMatrix)) {
  cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	  
  if(plotType=="boxplot") {
	iS <- intensScores[[i]]
	n <- length(iS)
	df <- diff(xval)[1]
	for(j in 1:length(iS)) {
	  xvals <- as.numeric(iS[[j]]$names)
	  bxp(iS[[j]], at=xval+(j-1)*df/n,pars=list(boxwex=.7*df/n,medcol=cols[j],boxcol=cols[j],whiskcol=cols[j],outcol=cols[j]),
		  add=(j>1),show.names=(j==1),xlim=c(xvals[1]-df/2,xvals[length(xvals)]+df/2), ...)
	}
  } else {
	dm <- intensScores[i,,]
        if (!sameScale) {
          rng <- range(dm, na.rm=TRUE)
          if (symmScale) rng <- c(-max(abs(rng)),max(abs(rng)))
        }

	titName <- paste("Signal:", colnames(dataMatrix)[i], label[orderingIndex[i]], colnames(ordering)[orderingIndex[i]], sep="")
	if(plotType=="line")
	{
		  layout(rbind(c(1, 2)), widths=c(3,2))
		  par(oma = c(0, 0, 2, 0))
		  par(mai=c(1.02,0.90,0.82,0))
		  matplot(xval,dm,type="l",col=cols,lty=lty,lwd=lwd,xlab="Position relative to TSS",ylab="Signal",ylim=rng)
		  par(mai=c(1.02,0.05,0.82,0))
		  plot.new()
		  legend(x="top", title ="Line Colours", col=cols, lty = 1, legend=cutLevels)
		  if (verbose) print(cols)
		  intervals <- breaks[[orderingIndex[i]]][["intervals"]]
		  if (verbose) print(intervals)
		  mtext(titName, line = 0.5, outer = TRUE)
	} else if(plotType=="heatmap") {
		  layout(rbind(c(1,2,3)), widths=c(1,3,1))
		  par(mai=c(1.02,0.50,0.82,0.05))
		  par(oma = c(0, 0, 0, 0))
		  image(rbind(1:nbins), col=cols,axes=F, xlab="Signal Intensity")
		  axis(2, at=(0:nbins)/nbins, labels=format(seq(rng[1], rng[2], length.out=nbins+1), digits=1))
		  par(mai=c(1.02,0.05,0.82,0.05))
		  image(xval,1:nbins,dm,xlab="Position relative to TSS", yaxt="n", ylab="Bin",col=cols,zlim=rng)
		  par(mai=c(1.02,0.05,0.82,0.50))
		  breakpoints <- breaks[[orderingIndex[i]]][["breakpoints"]]
		  plot(x=breakpoints,y=0:nbins, type="l", yaxt="n", lwd=3,xlab="log2 Expression", yaxs="i")
		  par(oma = c(0, 0, 2, 0))
		  mtext(titName, line = 0, outer = TRUE)
		} else if(plotType=="terrain") {
		  layout(1)
		  par(oma = c(0, 0, 2, 0))
  		  dm.avg <- (dm[-1, -1] + dm[-1, -(ncol(dm) - 1)] +
             		dm[-(nrow(dm) -1), -1] + dm[-(nrow(dm) -1), -(ncol(dm) - 1)]) / 4

  		  this.cols = cols[cut(dm.avg, breaks = seq(rng[1], rng[2], length.out=length(cols)), include.lowest = T)] 
		  persp(xval, 1:nbins, dm, xlab="Position relative to TSS", yaxt="n", ylab="Bin", col=this.cols, zlim=rng, theta=-25, phi=20, d=1.5, border=NA, ticktype="detailed", zlab="Signal")
		  mtext(titName, line = 0, outer = TRUE)
		}

	  }
	  par(def.par)#- reset to default
    }
})


