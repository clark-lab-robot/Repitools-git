processNDF <- function(filename, ncols=768) {
	ndfTemp <- read.table(filename, sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="", colClasses=c('NULL','NULL','NULL','NULL','NULL','character','NULL','NULL','NULL','NULL','NULL','NULL','character','integer','NULL','integer','integer'))
	
	#determine if reverse strand probe are marked by _RS or RS	
	if (length(grep("_RS", ndfTemp$PROBE_ID))>0) rsSymbol = "_RS" else rsSymbol="RS"

	chr <- gsub(paste("FS.*|",rsSymbol,".*",sep=""),"",ndfTemp$PROBE_ID)
	chr <- gsub("CHR", "chr", chr)
	chr <- gsub("chr0", "chr", chr)

	#return the centre of the probes sequence
	position <- ndfTemp$POSITION+trunc(nchar(ndfTemp$PROBE_SEQUENCE)/2)
	strand <- rep("+", nrow(ndfTemp))
	strand[grep(rsSymbol,ndfTemp$PROBE_ID)] = "-"

	#return the GC content of each probe
	GC <- sapply(ndfTemp$PROBE_SEQUENCE, function(x) {sum(strsplit(x,split="")[[1]] %in% c("C","G"))/nchar(x)}, USE.NAMES=FALSE)
	ndf <- data.frame(chr, position, strand, index=ndfTemp$Y*ncols+ndfTemp$X, sequence=ndfTemp$PROBE_SEQUENCE, GC, stringsAsFactors=FALSE)
	ndf <- ndf[order(ndf$chr, ndf$position, ndf$strand),]
	return(ndf)
}

loadPairFile <- function(filename, ndf, ncols=768) {
	pairTemp <- read.table(filename, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	pairTemp$PM[pairTemp$PM==0] = 1
	return(log2(pairTemp$PM[match(ndf$index,pairTemp$Y*ncols+pairTemp$X)]))
}

loadSampleDirectory <- function(path, ndf, what="Cy3", ncols=768) {
	#what="Cy3" or "Cy5" or "Cy3/Cy5" or "Cy5/Cy3"
	if (!what %in% c("Cy3","Cy5","Cy3/Cy5","Cy5/Cy3","Cy3andCy5","Cy5andCy3")) stop('parameter "what" must be "Cy3", "Cy5", "Cy3/Cy5", "Cy5/Cy3", "Cy3andCy5" or "Cy5andCy3"')	
	files <- dir(path, pattern=".pair$")
	samples <- unique(gsub("_532.pair|_635.pair","",files))
        if (what=="Cy3andCy5") 
          tempData <- matrix(NA, nrow=nrow(ndf), ncol=2*length(samples), dimnames=list(NULL, paste(rep(samples,2), rep(c("Cy3","Cy5"), each=length(samples)), sep="_")))
        else
          tempData <- matrix(NA, nrow=nrow(ndf), ncol=length(samples), dimnames=list(NULL, samples))
	for (i in 1:length(samples)) {
		if (what=="Cy3")
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_532.pair",sep=""), ndf, ncols)
		else if (what=="Cy5")
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_635.pair",sep=""), ndf, ncols)
		else if (what=="Cy3/Cy5")
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_532.pair",sep=""), ndf, ncols)-loadPairFile(paste(path,"/",samples[i],"_635.pair",sep=""), ndf, ncols)
		else if (what=="Cy5/Cy3")
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_635.pair",sep=""), ndf, ncols)-loadPairFile(paste(path,"/",samples[i],"_532.pair",sep=""), ndf, ncols)
		else if (what=="Cy3andCy5") {
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_532.pair",sep=""), ndf, ncols) 
                  tempData[,i+length(samples)] <- loadPairFile(paste(path,"/",samples[i],"_635.pair",sep=""), ndf, ncols)
                } else if (what=="Cy5andCy3") {
		  tempData[,i] <- loadPairFile(paste(path,"/",samples[i],"_635.pair",sep=""), ndf, ncols) 
                  tempData[,i+length(samples)] <- loadPairFile(paste(path,"/",samples[i],"_532.pair",sep=""), ndf, ncols)
                }
	}
	return(tempData)
}