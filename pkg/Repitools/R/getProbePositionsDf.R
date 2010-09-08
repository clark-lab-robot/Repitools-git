setMethodS3("getProbePositionsDf", "AffymetrixCdfFile", function(cdf, ..., verbose=-20) {

    require(aroma.affymetrix)

	chrs <- paste( "chr", c(1:22,"X","Y","M"), sep="" )
    names(chrs) <- 1:25
	
    ind <- getCellIndices(cdf,...,useNames=FALSE,unlist=TRUE,verbose=verbose)
    acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
    ch <- acp[ind,1,drop=TRUE]
    sp <- acp[ind,2,drop=TRUE]

    data.frame(chr=chrs[as.character(ch)],position=sp,index=ind,stringsAsFactors=FALSE)
})
