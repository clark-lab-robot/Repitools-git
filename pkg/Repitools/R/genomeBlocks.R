genomeBlocks <- function(organism, chrs, width, spacing=width) {
    chr.length <- seqlengths(organism)[chrs]
    windows <- do.call(c, lapply(names(chr.length), function(x) GRanges(seqnames = x, ranges = IRanges(start=seq.int(spacing/2, chr.length[x], spacing)-width/2+1, end=seq.int(spacing/2, chr.length[x], spacing)+width/2))))
    return(windows)
}
