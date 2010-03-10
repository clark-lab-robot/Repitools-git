writeWig <- function(rs.cov, filename, sample=20, dropZero=TRUE, normalise=1, verbose=TRUE) {
	if (class(rs.cov)!="GenomeData" | class(rs.cov[[1]])!="Rle") stop("rs.cov must be a coverage GenomeData object!")
	f1 <- if (grepl("gz$",filename)) gzfile(filename, open="wt") else file(filename, open="wt")
	scipen <- getOption("scipen")
	options(scipen=100)
	for (i in 1:length(rs.cov)) {
		cat("writing out",names(rs.cov)[i],"\n")
		writeLines(paste("variableStep chrom=", names(rs.cov)[i]," span=", sample, sep=""), f1)
		bp <- 1:(length(rs.cov[[i]]) %/% sample)*sample
		bp.score <- as.integer(rs.cov[[i]][bp])/normalise
		temp <- if (dropZero) paste(bp[bp.score!=0], bp.score[bp.score!=0], sep="\t")
			else paste(bp, bp.score, sep="\t")
		writeLines(temp, f1)
	}
	options(scipen=scipen)
	close(f1)
}

