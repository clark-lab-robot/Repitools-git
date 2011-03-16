setGeneric(".annoDF2GR", signature = "anno", function(anno, ...)
                                {standardGeneric(".annoDF2GR")})

setMethod(".annoDF2GR", "data.frame", function(anno)
{
    GRanges(anno$chr,
	    IRanges(anno$start, anno$end),
            if("strand" %in% colnames(anno)) anno$strand else '*',
            name = if("name" %in% colnames(anno)) anno$name)
})
