setGeneric("cpgBoxplots", function(this, ...){standardGeneric("cpgBoxplots")})

.cpgBoxplots <- function(dm, bins, gcContent, nBins, calcDiff, pdfFile, mfrow, col, ylim, gcCount, cb, sampleNames)
{
  if(calcDiff){
    title1 <- paste( col, paste(sampleNames,collapse="-"), sep="=" )
  }else{
    title1 <- paste( paste(col,sampleNames,sep="="), collapse="," )
	}
	
  if( !is.null(pdfFile) ) {
    pdf(pdfFile,width=10,height=10)
	par(mfrow=mfrow)
  }
	
  actualNBins <- length( levels(bins) )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Loop through G+C contents and 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  bp <- NULL
  count <- 0
  for(j in gcContent) {
    w <- which( gcCount==j)
	if( length(w)==0 )
	  next
	  
    title2 <- paste("[Probe G+C =",j,"]","Percentage of Probes:",round(length(w)/length(cb)*100,2))
	main <- paste(title1,title2,sep="\n")
	
	if (calcDiff) {
    count <- count+1
	  bp[[count]] <- boxplot( dm[w,1]-dm[w,2] ~ bins[w], xlim=c(0,actualNBins),ylim=ylim,col=col,main=main,las=2,cex.axis=.8,cex.main=.9)
    names(bp) <- paste(colnames(dm), collapse="-")
	} else {
    count <- count+1
    bp[[count]] <- boxplot( dm[w,1] ~ bins[w], at=(1:actualNBins)/2, boxwex=.4,xlim=c(0,actualNBins),ylim=ylim,col=col[1],main=main,las=2,cex.axis=.8,cex.main=.9)
    count <- count+1
    bp[[count]] <- boxplot( dm[w,2] ~ bins[w], at=((actualNBins+1):(actualNBins*2))/2,boxwex=.4, main="",col=col[2],add=TRUE,las=2,cex.axis=.8)
    names(bp) <- colnames(dm)
	}
  }
  
  if( !is.null(pdfFile) )
    dev.off()
	
  invisible(bp)
}

.createBins <- function(u, nBins) {
	q<-quantile(u,prob=(0:nBins)/nBins)
	q[1] <- q[1]-.000000001
	n <- length(q)
	q[n] <- q[n]+.000000001
	cut(u,breaks=q)
}

setMethod("cpgBoxplots", "AffymetrixCelSet", function(this, samples=c(1,2), subsetChrs="chr[1-5]", gcContent=7:18, 
                                                     calcDiff=FALSE, verbose=FALSE, nBins=40, pdfFile=NULL,
													 ylim=if (calcDiff) c(-5,6) else c(4,15), 
													 col=if (calcDiff) "salmon" else c("lightgreen","lightblue"),
													 mfrow=if (!is.null(pdfFile)) c(2,2) else c(1,1)) {
													 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  
  if( length(samples) != 2 )
    stop("Can only do boxplots on a pair of samples.")
	
  if(calcDiff && length(col) != 1)
    stop("calcDiff=TRUE, but length(col) != 1.")
	
  if(!calcDiff && length(col) != 2)
    stop("calcDiff=FALSE, but length(col) != 2.")
	
  if( max(samples) > nbrOfArrays(this) )
    stop("'samples' is out of range.")
  	
  cdf <- getCdf(this)
  mainCdf <- getMainCdf(cdf)
  
  if (is.null(subsetChrs))
    units <- seq_len(nbrOfUnits(cdf))
  else
    units <- indexOf(cdf,subsetChrs)
	
  if( length(units) == 0 )
    stop("'units' is length 0.  Specify an appropriate 'subsetChrs' argument.")
	
	 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, sprintf("Reading indices for %d units (unique CDF) ",length(units)));
  indices <- getCellIndices(cdf,units=units,stratifyBy="pm",verbose=verbose)
  indices <- unlist(indices,use.names=FALSE)
  verbose && exit(verbose);
  
  verbose && enter(verbose, sprintf("Reading indices for %d units (main CDF) ",length(units)));
  mainIndices <- getCellIndices(mainCdf,units=units,stratifyBy="pm",verbose=verbose)
  mainIndices <- unlist(mainIndices,use.names=FALSE)
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Counting bases
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, sprintf("Counting bases for %d probes",length(mainIndices)));
  acs <- AromaCellSequenceFile$byChipType(getChipType(mainCdf))
  cb <- countBases(acs,cells=mainIndices)
  gcCount <- rowSums( cb[,c("C","G")] )
  verbose && exit(verbose);
  
  cs <- extract(this,samples)
  sampleNames <- getNames(cs)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reading data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Reading intensity data");
  dm <- extractMatrix(cs,cells=indices,verbose=verbose)
  dm <- log2(dm)
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reading CpG density data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Reading and binning Cpg data");
  acc <- AromaCellCpgFile$byChipType(getChipType(cdf))
  cpgDens <- acc[indices,1,drop=TRUE]
  bins <- .createBins(cpgDens, nBins)
  verbose && exit(verbose);

  .cpgBoxplots(dm, bins, gcContent, nBins, calcDiff, pdfFile, mfrow, col, ylim, gcCount, cb, sampleNames)												 
} 
)

setMethod("cpgBoxplots", "matrix", function(this, ndfTable, organism, samples=c(1,2), subsetChrs="chr[1-5]", gcContent=7:18, 
                                                     calcDiff=FALSE, verbose=FALSE, nBins=40, pdfFile=NULL,
													 ylim=if (calcDiff) c(-5,6) else c(4,15), 
													 col=if (calcDiff) "salmon" else c("lightgreen","lightblue"),
													 mfrow=if (!is.null(pdfFile)) c(2,2) else c(1,1)) {
													 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  
  if( length(samples) != 2 )
    stop("Can only do boxplots on a pair of samples.")
	
  if(calcDiff && length(col) != 1)
    stop("calcDiff=TRUE, but length(col) != 1.")
	
  if(!calcDiff && length(col) != 2)
    stop("calcDiff=FALSE, but length(col) != 2.")
	
  if( max(samples) > ncol(this) )
    stop("'samples' is out of range.")
  
  if (is.null(subsetChrs))
    usefulProbeIndices <- 1:nrow(ndfTable)
  else
    usefulProbeIndices <- grep(subsetChrs, ndfTable$chr)
	
  if( length(usefulProbeIndices) == 0 )
    stop("'units' is length 0.  Specify an appropriate 'subsetChrs' argument.")
	
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Counting bases
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, sprintf("Counting bases for %d probes",length(usefulProbeIndices)))
  ndfTable <- ndfTable[usefulProbeIndices, ]
  gcCount <- sapply(gregexpr("[CG]", ndfTable$sequence), length)
  cb <- sapply(ndfTable$sequence, length)
  verbose && exit(verbose)

  densities <- cpgDensityCalc(ndfTable, 300, organism = organism)
  bins <- .createBins(densities, nBins) 	
  
  sampleNames <- colnames(this)[samples]
  
  .cpgBoxplots(this[usefulProbeIndices, samples], bins, gcContent, nBins, calcDiff, pdfFile, mfrow, col, ylim, gcCount, cb, sampleNames)
}
)
