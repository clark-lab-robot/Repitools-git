genomeBlocks <- function(organism, chrs, width, spacing=width) {
    chr.length <- seqlengths(organism)[chrs]
    temp <- lapply(names(chr.length), function(x) IRanges(start=seq.int(spacing/2, chr.length[x], spacing)-width/2+1, end=seq.int(spacing/2, chr.length[x], spacing)+width/2))
    names(temp) <- names(chr.length)
    RangedData(RangesList(temp))        
}
