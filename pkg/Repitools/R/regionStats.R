#regionStats <- function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, fdrProbes = FALSE, ...) UseMethod("regionStats")

#regionStats <- function(...) UseMethod("regionStats")

.regionStats <- function(diffs, design, ch, sp, fdrLevel, nPermutations, probeWindow, meanTrim, nProbes, maxGap, twoSides, verbose, returnTrimmedMeans) {


  getBed <- function(score, ch, sp, cut=NULL, nProbes=10, maxGap, twoSides) {
    if( is.null(cut) )
      stop("Need to specify 'cut'.")
    posInd <- getRegions(score, ch, sp, nProbes, maxGap, cut, twoSides, doJoin=TRUE)
    if( is.null(posInd) )
      return(list())
    posReg <- data.frame(chr=paste("chr",ch[posInd$start],sep=""),
                         start=sp[posInd$start],
	                 end=sp[posInd$end], score=0, 
                         startInd=posInd$start, 
                         endInd=posInd$end, 
                         stringsAsFactors=F)
    for(i in 1:nrow(posInd))
      posReg$score[i] <- round(median(score[ (posInd$start[i]:posInd$end[i]) ]),3)
    posReg
  }

  
  getRegions <- function(score, ch, sp, nProbes, maxGap, cutoff, twoSides, doJoin) {
    getRegionsChr <- function(ind, score, sp, nProbes, maxGap, cutoff, doJoin) {
#pad the beginning & end
      probes <- c(FALSE, score[ind] > cutoff, FALSE)

      #insert FALSEs in to break up regions with gaps>maxGap
      probeGaps <- which(diff(sp[ind])>maxGap)
      num.gaps <- length(probeGaps)

      ind.2 <- rep(NA, length(ind)+num.gaps)
      ind.gaps <- probeGaps+1:num.gaps
      ind.nogaps <- (1:length(ind.2))[-ind.gaps]
      ind.2[ind.nogaps] <- ind

      probes.2 <- rep(FALSE, length(probes)+num.gaps)
      probes.gaps <- probeGaps+1:num.gaps+1
      probes.nogaps <- (1:length(probes.2))[-probes.gaps]
      probes.2[probes.nogaps] <- probes
 
      df <- diff(probes.2)
      st <- ind.2[which(df==1)]
      en <- ind.2[which(df==-1)-1]

      #sort starts & ends again
      st <- sort(st)
      en <- sort(en)

      #join regions with < maxGap basepairs between positive probes
      if (doJoin) {
        gap.w <- which(sp[st[-1]]-sp[rev(rev(en)[-1])] < maxGap)
        if (length(gap.w)>0) {
          st <- st[-(gap.w+1)]
          en <- en[-gap.w]
        }
      }

      w <- (en-st+1) >= nProbes
      if (sum(w)==0)
        return(data.frame(start=NULL, end=NULL))
      else
        data.frame(start=st,end=en)[w,]
    }

    chrInds <- split(1:length(score), ch)
    regTable <- data.frame(start=NULL, end=NULL)
    for (i in 1:length(chrInds)) regTable <- rbind(regTable, getRegionsChr(chrInds[[i]], score, sp, nProbes, maxGap, cutoff, doJoin))
    if (twoSides) for (i in 1:length(chrInds)) regTable <- rbind(regTable, getRegionsChr(chrInds[[i]], -score, sp,  nProbes, maxGap, cutoff, doJoin))
    return(regTable)
  }

  
  fdrTable <- function(realScore, permScore, ch, sp, cutsLength, nProbes, maxGap, twoSides, minCutoff = .5, maxCutoff=max( abs(permScore), na.rm=TRUE ), verbose) {
    require(gsmoothr)
    cuts <- seq(minCutoff,maxCutoff,length=cutsLength)

    fdr <- matrix(,nr=length(cuts),nc=4)
    colnames(fdr) <- c("cutoff","neg","pos","fdr")
    for(i in 1:length(cuts)) {
      pos <- nrow(getRegions(realScore, ch, sp, nProbes, maxGap, cuts[i], twoSides, doJoin=FALSE))
      neg <- nrow(getRegions(permScore, ch, sp, nProbes, maxGap, cuts[i], twoSides, doJoin=FALSE))
      fdr[i,] <- c(cuts[i],neg,pos,min(neg/pos,1))
      if (verbose) cat(".")
    }
    if (verbose) cat("\n")
    as.data.frame(fdr)
  }

  
  uch <- unique(ch)

  tmeanReal <- matrix(,nr=nrow(diffs),nc=ncol(diffs))
  tmeanPerms <- lapply( as.list(colnames(design)), FUN=function(u) {
    matrix(NA,nr=nrow(diffs),nc=nPermutations)
  })
  regions <- fdrTabs <- vector("list", length(tmeanPerms))
  names(regions) <- names(fdrTabs) <- colnames(design)
  
  ifelse( verbose, print(gc()), gc())
  
  
  # calculate smoothed statistics
  for(col in 1:ncol(diffs)) {
    if( verbose )
	  cat("Calculating trimmed means for column", col, "of design matrix:\n")
    for(ii in 1:length(uch)) {
      if( verbose )
	    cat(" ", uch[ii], "-", sep="")
	  w <- which(ch == uch[ii])
	  
	  tmeanReal[w,col] <- gsmoothr::tmeanC(sp[w], diffs[w,col], probeWindow=probeWindow, trim=meanTrim, nProbes=nProbes)
  	  if( verbose )
	      cat("R")
	  for(j in 1:ncol(tmeanPerms[[col]])) {
	    s <- sample(1:nrow(tmeanReal))
	    tmeanPerms[[col]][w,j] <- gsmoothr::tmeanC(sp[w], diffs[s,col][w], probeWindow=probeWindow, trim=meanTrim, nProbes=nProbes)
		if( verbose )
	      cat(".")

	  }
	}


    if( verbose )
      cat("\nCalculating FDR table.\n")
    # calculate FDR table
	mx <- max(abs(tmeanPerms[[col]]),na.rm=TRUE)


	z <- apply(tmeanPerms[[col]], 2, FUN=function(u) fdrTable(tmeanReal[,col], u, ch, sp, 40, nProbes, maxGap, twoSides, maxCut=mx, verbose=verbose))
	fdrTabs[[col]] <- z[[1]]
	for(i in 2:length(nPermutations))
	  fdrTabs[[col]][,2:3] <- fdrTabs[[col]][,2:3] + z[[i]][,2:3]
	fdrTabs[[col]]$fdr <- pmin(fdrTabs[[col]]$neg/fdrTabs[[col]]$pos,1)  # re-adjust FDR calculation over all permutations
	
	# select lowest cutoff such that FDR is achieved
	w <- which(fdrTabs[[col]]$fdr < fdrLevel )
	cut <- min( fdrTabs[[col]]$cut[w], na.rm=TRUE )
	
    if( verbose )
      cat("Using cutoff of", cut, "for FDR of", fdrLevel,"\n")
	  
	regions[[col]] <- getBed(tmeanReal[,col], ch, sp, cut, nProbes, maxGap, twoSides)
	
	
  }
  
  if(returnTrimmedMeans == TRUE)
  {
	return(list(regions=regions,tmeanReal=tmeanReal,tmeanPerms=tmeanPerms,fdrTables=fdrTabs))
  } else {
	return(list(regions=regions, fdrTables=fdrTabs))
  }

}


setMethodS3("regionStats","AffymetrixCelSet",function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ind=NULL, returnTrimmedMeans = FALSE, ...) {

  require(aroma.affymetrix)

  # cs - AffymetrixCelSet to read probe-level data from
  # design - design matrix
  # ind - (optional) 

  
  cdf <- getCdf(cs)
    
  if( is.null(ind) )
    ind <- getCellIndices( cdf, useNames=FALSE, unlist=TRUE)

  if( nrow(design) != nbrOfArrays(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
	
  acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
  ch <- acp[ind,1,drop=TRUE]
  sp <- acp[ind,2,drop=TRUE]
  
  # cut down on the amount of data read, if some rows of the design matrix are all zeros
  w <- which( rowSums(design != 0) > 0 )
  cs <- extract(cs,w, verbose=verbose)
  dmP <- log2(extractMatrix(cs,cells=ind,verbose=verbose))
  
  # compute probe-level score of some contrast
  diffs <- dmP %*% design[w,]

  w <- rowSums( is.na(diffs) )==0
  if( verbose )
    cat("Removing", sum(!w), "rows, due to NAs.\n")
	
  diffs <- diffs[w,,drop=FALSE]
  ch <- ch[w]
  sp <- sp[w]
  
  rm(dmP)
  ifelse( verbose, print(gc()), gc())

  return(.regionStats(diffs, design, ch, sp, fdrLevel, nPermutations, probeWindow, meanTrim, nProbes, maxGap, twoSides, verbose, returnTrimmedMeans))
})


setMethodS3("regionStats","default",function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, maxGap=500, twoSides=TRUE, verbose=TRUE, ndf, returnTrimmedMeans = FALSE, ...) {
  #nimblegen data

  # cut down on the amount of data read, if some rows of the design matrix are all zeros
  w <- which( rowSums(design != 0) > 0 )
  diffs = cs %*% design

  w <- rowSums( is.na(diffs) )==0
  if( verbose )
    cat("Removing", sum(!w), "rows, due to NAs.\n")


  return(.regionStats(diffs, design, gsub("chr","",ndf$chr), ndf$position, fdrLevel, nPermutations, probeWindow, meanTrim, nProbes, maxGap, twoSides, verbose, returnTrimmedMeans))
})

