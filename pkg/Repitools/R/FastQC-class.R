setClass("FastQC",
    representation(
        Basic_Statistics="data.frame", 
        Per_base_sequence_quality="data.frame",
        Per_sequence_quality_scores="data.frame",
        Per_base_sequence_content="data.frame",
        Per_base_GC_content="data.frame",
        Per_sequence_GC_content="data.frame",
        Per_base_N_content="data.frame",
        Sequence_Length_Distribution="data.frame",
        Sequence_Duplication_Levels="data.frame",
        Overrepresented_sequences="data.frame"
    )
)

setMethod("show", "FastQC", 
    function(object) {
        cat("FastQC object of ", object@Basic_Statistics$Value[3], " ", object@Basic_Statistics$Value[4], "bp sequences, containing:\n", sep = "")
        cat(slotNames(object), sep="\n")
    }
)


setClass("SequenceQC",
    representation(
        Unaligned="FastQC",
        Aligned="FastQC",
        Mismatches="data.frame"
    )
)

setMethod("show", "SequenceQC", function(object) {
    cat("SequenceQC object of ", object@Unaligned@Basic_Statistics$Value[4], "bp sequences\n", sep="")
    cat(object@Unaligned@Basic_Statistics$Value[3], " before alignment\n", sep="")
    cat(object@Aligned@Basic_Statistics$Value[3], " after alignment\n", sep="")
})

setClass("SequenceQCSet",
    representation(NAMES="character"),
    contains="list",
    prototype=list()
)

#Constructor
SequenceQCSet <- function(x) {
    new("SequenceQCSet", x, NAMES=names(x))
}

setValidity("SequenceQCSet", function(object) {
    if (!all(sapply(object@.Data, class)=="SequenceQC"))
        return("All elements must be 'SequenceQC' objects")
    TRUE
})

setMethod("show", "SequenceQCSet", function(object) {
    cat("SequenceQCSet object containing ", length(object@.Data), " element(s)\n", sep="")
    cat(paste(object@NAMES, collapse=" "), "\n")
})

setMethod("names", "SequenceQCSet", function(x) x@NAMES)
setMethod("[", "SequenceQCSet", function(x, i, j, ..., drop) {
    if (is.character(i)) i <- match(i, x@NAMES)
    new("SequenceQCSet", x@.Data[i], NAMES=x@NAMES[i])
})

setMethod("[[", "SequenceQCSet", function(x, i, j, ...) {
    if (is.character(i)) i <- match(i, x@NAMES)
    x@.Data[[i]]
})

readFastQC <- function(filename) {
    if (length(filename)!=1) stop("Length of filename must be 1")
    if (!file.exists(filename)) stop(filename, "does not exist")
    temp <- readLines(filename)
    #remove #'s
    temp <- gsub("#", "", temp)
    #remove END_MODULE lines
    temp <- temp[!grepl(">>END_MODULE", temp)]
    #split into modules
    temp <- split(temp, cumsum(grepl("^>>", temp)))[-1]
    names(temp) <- sapply(temp, function(x) {
        gsub("^>>", "", gsub("\t.*", "", gsub(" ", "_", x[1])))
    })
    #Parse each module
    temp <- lapply(temp, function(x) {
        if (length(x)==1) return(data.frame())
        x <- strsplit(x[-1], split="\t")
        tab <- as.data.frame(do.call(rbind, x[-1]), stringsAsFactors=FALSE)
        for (i in 1:ncol(tab)) if (!any(is.na(suppressWarnings(as.numeric(tab[,i]))))) tab[,i] <- as.numeric(tab[,i]) 
        colnames(tab) <- x[[1]]
        tab
    })
    do.call(new, c(list(Class="FastQC"), temp))
}
